#pragma once
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include <map>
#include <cstring>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include "core.hpp"
#include "functions.hpp"
#include "RgswSample.hpp"
#include "LweSample.hpp"

// How the database is kept on SSD when write_to_file is set.
//   kSsdLegacyNtt (0):     the NTT-form database (one uint64 per coefficient,
//                          DEGREE1 coefficients per polynomial) is written and
//                          read back verbatim. 8 bytes per 1-byte plaintext
//                          coefficient, and capped at kMaxDbPolys polynomials.
//   kSsdCompactStaged (1): the raw plaintext (one byte per coefficient, P=256)
//                          of ALL ROW*COL polynomials is written. Load = one
//                          sequential read into a staging buffer, then widen
//                          to uint64 and forward-NTT into the working buffer.
//   kSsdCompactStream (2): same file as (1), but loaded in fixed-size chunks
//                          with POSIX_FADV_WILLNEED prefetch of the next chunk
//                          so the kernel's readahead overlaps the (single
//                          threaded) widen+NTT work.
enum VIASsdFormat : int { kSsdLegacyNtt = 0, kSsdCompactStaged = 1, kSsdCompactStream = 2 };

class VIAServer
{
private:
    uint64_t* contiguous_db = nullptr;
    bool write_to_file_ = false;
    int ssd_format_ = kSsdCompactStaged;
    uint64_t db_checksum_ = 0;
    int64_t last_load_us_ = 0;
    std::vector<uint64_t> retain_polys_;
    std::map<uint64_t, std::vector<uint8_t>> retained_;
    static constexpr const char* kDbPath = "/tmp/via_db.bin";
    // RAM cap on the number of NTT-form polynomials held in memory (32 GiB).
    static constexpr uint64_t kMaxDbPolys = UINT64_C(1) << 21;
    // Polynomials processed per widen/NTT batch (64 * 16 KiB = 1 MiB).
    static constexpr uint64_t kBatchPolys = 64;
    // Chunk size for kSsdCompactStream (16 MiB of plaintext = 8192 polys).
    static constexpr uint64_t kStreamChunkPolys = 8192;

    static uint64_t numPolys() { return ROW * COL; }
    static uint64_t numSlots() { return (ROW * COL > kMaxDbPolys) ? kMaxDbPolys : (ROW * COL); }

    static bool writeAll(int fd, const char* ptr, uint64_t n)
    {
        while (n > 0) {
            ssize_t w = write(fd, ptr, n);
            if (w <= 0) return false;
            ptr += w;
            n -= static_cast<uint64_t>(w);
        }
        return true;
    }

    static uint64_t readAll(int fd, char* ptr, uint64_t n)
    {
        uint64_t total = 0;
        while (n > 0) {
            ssize_t r = read(fd, ptr, n);
            if (r <= 0) break;
            ptr += r;
            n -= static_cast<uint64_t>(r);
            total += static_cast<uint64_t>(r);
        }
        return total;
    }

    // Widen one plaintext polynomial (bytes, values in [0, P)) to uint64 words.
    static inline void widenPoly(const uint8_t* __restrict src, uint64_t* __restrict dst)
    {
        for (uint64_t c = 0; c < DEGREE1; c++) dst[c] = src[c];
    }

    // Widen + forward NTT `count` plaintext polys starting at global poly index
    // `first`. Slot = poly index mod numSlots() (same mapping the First
    // Dimension uses via `% fdDbSize`).
    void encodeBatch(const uint8_t* src, uint64_t first, uint64_t count, int64_t& parse_ns, int64_t& derive_ns)
    {
        const uint64_t slots = numSlots();
        auto t0 = std::chrono::steady_clock::now();
        for (uint64_t q = 0; q < count; q++)
            widenPoly(src + q * DEGREE1, database[(first + q) % slots]);
        auto t1 = std::chrono::steady_clock::now();
        for (uint64_t q = 0; q < count; q++)
        {
            uint64_t* p = database[(first + q) % slots];
            nttQ2Deg1.ComputeForward(p, p, 1, 1);
        }
        auto t2 = std::chrono::steady_clock::now();
        parse_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
        derive_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    }

    uint64_t checksumDb() const
    {
        const uint64_t n = numSlots() * DEGREE1;
        uint64_t a0 = 0, a1 = 0, a2 = 0, a3 = 0;
        for (uint64_t i = 0; i + 3 < n; i += 4)
        {
            a0 = a0 * 0x9E3779B97F4A7C15ULL + contiguous_db[i];
            a1 = a1 * 0xC2B2AE3D27D4EB4FULL + contiguous_db[i + 1];
            a2 = a2 * 0x165667B19E3779F9ULL + contiguous_db[i + 2];
            a3 = a3 * 0x27D4EB2F165667C5ULL + contiguous_db[i + 3];
        }
        return a0 ^ (a1 << 1) ^ (a2 << 2) ^ (a3 << 3);
    }

    static uint64_t nttTableBytes(const intel::hexl::NTT& t)
    {
        return 8 * (t.GetRootOfUnityPowers().size() + t.GetPrecon32RootOfUnityPowers().size() +
                    t.GetPrecon64RootOfUnityPowers().size() + t.GetAVX512RootOfUnityPowers().size() +
                    t.GetAVX512Precon32RootOfUnityPowers().size() + t.GetAVX512Precon52RootOfUnityPowers().size() +
                    t.GetAVX512Precon64RootOfUnityPowers().size() + t.GetInvRootOfUnityPowers().size() +
                    t.GetPrecon32InvRootOfUnityPowers().size() + t.GetPrecon52InvRootOfUnityPowers().size() +
                    t.GetPrecon64InvRootOfUnityPowers().size());
    }

public:
    void getExpectedPolynomial(uint64_t idx_row_col, uint64_t* out) {
        // HEXL signature is ComputeInverse(result, operand, ...); the original
        // call had them swapped, which zeroed the DB entry and returned zeros.
        nttQ2Deg1.ComputeInverse(out, database[idx_row_col % numSlots()], 1, 1);
    }

    // Ground truth for verification: plaintext bytes of selected polynomials,
    // captured while the DB is generated (before any test corruption). With
    // the RAM cap, the last polynomial mapped to a slot is the one kept.
    void setRetainPolys(const std::vector<uint64_t>& polys) { retain_polys_ = polys; }
    bool getRetainedPlaintext(uint64_t poly, std::vector<uint8_t>& out) const
    {
        auto it = retained_.find(poly % numSlots());
        if (it == retained_.end()) return false;
        out = it->second;
        return true;
    }

    // Negative-test hook: perturb coefficient `coeff` of polynomial `poly` by
    // +1 in the stored DB (RAM working buffer in RAM mode, the SSD file in SSD
    // mode). Call after SetupDatabase(); not timed.
    void corruptDatabaseForTest(uint64_t poly, uint64_t coeff)
    {
        const uint64_t slot = poly % numSlots();
        std::vector<uint64_t> tmp(DEGREE1);
        auto bump = [&](uint64_t* ntt_poly) {
            nttQ2Deg1.ComputeInverse(tmp.data(), ntt_poly, 1, 1);
            tmp[coeff] = (tmp[coeff] + 1) % MODULUS_Q2;
            nttQ2Deg1.ComputeForward(ntt_poly, tmp.data(), 1, 1);
        };
        if (!write_to_file_) {
            bump(database[slot]);
        } else {
            int fd = open(kDbPath, O_RDWR);
            if (fd == -1) { std::cerr << "corrupt: open failed" << std::endl; return; }
            if (ssd_format_ == kSsdLegacyNtt) {
                std::vector<uint64_t> w(DEGREE1);
                off_t off = static_cast<off_t>(slot * DEGREE1 * sizeof(uint64_t));
                if (pread(fd, w.data(), DEGREE1 * sizeof(uint64_t), off) != static_cast<ssize_t>(DEGREE1 * sizeof(uint64_t))) std::cerr << "corrupt: pread failed" << std::endl;
                bump(w.data());
                if (pwrite(fd, w.data(), DEGREE1 * sizeof(uint64_t), off) != static_cast<ssize_t>(DEGREE1 * sizeof(uint64_t))) std::cerr << "corrupt: pwrite failed" << std::endl;
            } else {
                uint8_t b = 0;
                off_t off = static_cast<off_t>(poly * DEGREE1 + coeff);
                if (pread(fd, &b, 1, off) != 1) std::cerr << "corrupt: pread failed" << std::endl;
                b = static_cast<uint8_t>(b + 1);
                if (pwrite(fd, &b, 1, off) != 1) std::cerr << "corrupt: pwrite failed" << std::endl;
            }
            fdatasync(fd);
            posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
            close(fd);
        }
        std::cout << "│   ► TEST: corrupted DB poly " << poly << " coeff " << coeff << (write_to_file_ ? " (on SSD)" : " (in RAM)") << std::endl;
    }
    void setWriteToFile(bool flag) { write_to_file_ = flag; }
    void setSsdFormat(int f) { ssd_format_ = f; }

    // Compares the NTT-form working buffer (after the timed load) with the
    // checksum taken at setup time, before the buffer was wiped. Not timed.
    bool VerifyLoadedDatabase()
    {
        if (!write_to_file_) return true;
        bool ok = (checksumDb() == db_checksum_);
        std::cout << "│   ► Server Load Integrity Check: " << (ok ? "OK" : "FAILED") << std::endl;
        return ok;
    }

	// constructor
	VIAServer() : nttQ1_1(DEGREE1, MODULUS_Q1_1), nttQ1_2(DEGREE1, MODULUS_Q1_2), nttQ2(DEGREE2, MODULUS_Q2), nttQ2Deg1(DEGREE1, MODULUS_Q2), nttQ3(DEGREE2, MODULUS_Q3)
    {
        FirstDimVec.resize(ROW);
        FirstDimResMasks.resize(COL);
        FirstDimResVals.resize(COL);
        for (uint64_t i = 0; i < COL; i++)
        {
            FirstDimResMasks[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
            FirstDimResVals[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        }

        SecondDimVec.resize(COL);
        temp_message1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        temp_message1_1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        temp_message1_2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        decomp_scratch_ = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        temp_prod1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        gadget_poly1.resize(GADGET_L_MAX);
        gadget_poly2.resize(GADGET_L_MAX);
        for (int i = 0; i < GADGET_L_MAX; i++)
        {
            gadget_poly1[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
            gadget_poly2[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        }

        // NTT-form working buffer (server RAM). Contents are produced by
        // SetupDatabase() and, in SSD mode, re-produced by LoadDatabase().
        uint64_t dbSize = numSlots();
        database.resize(dbSize);
        contiguous_db = static_cast<uint64_t*>(std::aligned_alloc(4096, dbSize * DEGREE1 * sizeof(uint64_t)));
        for (uint64_t i = 0; i < dbSize; i++)
        {
            database[i] = contiguous_db + (i * DEGREE1);
        }
    }

	// destructor
	~VIAServer()
    {
        std::free(temp_message1);
        std::free(temp_message1_1);
        std::free(temp_message1_2);
        std::free(temp_prod1);
        std::free(decomp_scratch_);
        for (uint64_t i = 0; i < COL; i++)
        {
            std::free(FirstDimResMasks[i]);
            std::free(FirstDimResVals[i]);
        }
        for (int i = 0; i < GADGET_L_MAX; i++)
        {
            std::free(gadget_poly1[i]);
            std::free(gadget_poly2[i]);
        }
        if (contiguous_db != nullptr) {
            std::free(contiguous_db);
            contiguous_db = nullptr;
        }
    }

	// init_Database
	// Generates a random plaintext database (ROW*COL polynomials of DEGREE1
	// coefficients mod P = 256, i.e. 8 bits per coefficient as reported in
	// "Database size") and preprocesses it into NTT form mod Q2. In SSD mode,
	// also writes the on-disk representation, syncs it, evicts it from the page
	// cache and wipes the in-RAM NTT buffer so that the answer can only use what
	// LoadDatabase() reads back.
	void SetupDatabase()
    {
        std::cout << "\u2502   \u25BA" << " Setup Database: ";
        const uint64_t nPolys = numPolys();
        const uint64_t slots = numSlots();
        const bool compact = write_to_file_ && ssd_format_ != kSsdLegacyNtt;

        int fd = -1;
        if (compact) {
            fd = open(kDbPath, O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (fd == -1) std::cerr << "Failed to open " << kDbPath << " for writing!" << std::endl;
        }

        std::mt19937_64 rng(0x5649415044ULL);
        std::vector<uint8_t> chunk(kStreamChunkPolys * DEGREE1);
        int64_t parse_ns = 0, derive_ns = 0;
        bool write_ok = true;
        for (uint64_t p = 0; p < nPolys; p += kStreamChunkPolys)
        {
            uint64_t cnt = std::min<uint64_t>(kStreamChunkPolys, nPolys - p);
            uint64_t nbytes = cnt * DEGREE1;
            for (uint64_t b = 0; b < nbytes; b += 8)
            {
                uint64_t r = rng();
                std::memcpy(chunk.data() + b, &r, 8);
            }
            if (fd != -1) write_ok = write_ok && writeAll(fd, reinterpret_cast<const char*>(chunk.data()), nbytes);
            if (!retain_polys_.empty())
                for (uint64_t q = 0; q < cnt; q++)
                    for (uint64_t r : retain_polys_)
                        if ((p + q) % slots == r % slots)
                            retained_[r % slots].assign(chunk.data() + q * DEGREE1, chunk.data() + (q + 1) * DEGREE1);
            for (uint64_t q = 0; q < cnt; q += kBatchPolys)
                encodeBatch(chunk.data() + q * DEGREE1, p + q, std::min<uint64_t>(kBatchPolys, cnt - q), parse_ns, derive_ns);
        }
        std::cout << static_cast<double>(parse_ns + derive_ns) / 1e9 << " s" << std::endl;
        std::cout << "\u2502 " << std::endl;

        if (!write_to_file_) return;

        uint64_t total_size = 0;
        if (ssd_format_ == kSsdLegacyNtt) {
            // Legacy: write the NTT-form buffer (only numSlots() polys).
            fd = open(kDbPath, O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (fd == -1) {
                std::cerr << "Failed to open " << kDbPath << " for writing!" << std::endl;
                return;
            }
            total_size = slots * DEGREE1 * sizeof(uint64_t);
            write_ok = writeAll(fd, reinterpret_cast<const char*>(contiguous_db), total_size);
        } else {
            total_size = nPolys * DEGREE1;
        }
        if (fd != -1) {
            fdatasync(fd);
            posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
            close(fd);
        }
        if (!write_ok) std::cerr << "Short write to " << kDbPath << "!" << std::endl;

        const char* fmt_name = (ssd_format_ == kSsdLegacyNtt) ? "legacy-ntt-u64"
                             : (ssd_format_ == kSsdCompactStaged) ? "compact-u8-plaintext-staged"
                             : "compact-u8-plaintext-stream";
        std::cout << "│   ► Server SSD Format: " << fmt_name << std::endl;
        std::cout << "│   ► Server Long term Storage: " << total_size << " bytes" << std::endl;
        // Only DB-independent, parameter-derived state stays in RAM: HEXL NTT
        // twiddle tables (no seeds/pads in VIA's server).
        uint64_t seed_bytes = nttTableBytes(nttQ1_1) + nttTableBytes(nttQ1_2) + nttTableBytes(nttQ2) +
                              nttTableBytes(nttQ2Deg1) + nttTableBytes(nttQ3);
        std::cout << "│   ► Server Seed-Derived State In Memory: " << seed_bytes << " bytes" << std::endl;

        // Remember what the working buffer must look like after a load, then
        // wipe it (pages stay mapped) so a stale copy can't be used.
        db_checksum_ = checksumDb();
        std::memset(contiguous_db, 0, slots * DEGREE1 * sizeof(uint64_t));
    }

    void LoadDatabase()
    {
        last_load_us_ = 0;
        if (!write_to_file_) return;

        const uint64_t nPolys = numPolys();
        const uint64_t slots = numSlots();

        // Make sure the load is cold even on repeated rounds (not timed).
        {
            int efd = open(kDbPath, O_RDONLY);
            if (efd != -1) { posix_fadvise(efd, 0, 0, POSIX_FADV_DONTNEED); close(efd); }
        }

        int64_t read_ns = 0, parse_ns = 0, derive_ns = 0;
        uint64_t got = 0, want = 0;
        auto start = std::chrono::steady_clock::now();
        int fd = open(kDbPath, O_RDONLY);
        if (fd == -1) {
            std::cerr << "Failed to open " << kDbPath << " for reading!" << std::endl;
        } else if (ssd_format_ == kSsdLegacyNtt) {
            want = slots * DEGREE1 * sizeof(uint64_t);
            got = readAll(fd, reinterpret_cast<char*>(contiguous_db), want);
            close(fd);
            read_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count();
        } else if (ssd_format_ == kSsdCompactStaged) {
            want = nPolys * DEGREE1;
            uint8_t* staging = static_cast<uint8_t*>(std::aligned_alloc(4096, want));
            got = readAll(fd, reinterpret_cast<char*>(staging), want);
            close(fd);
            read_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count();
            for (uint64_t p = 0; p < nPolys; p += kBatchPolys)
                encodeBatch(staging + p * DEGREE1, p, std::min<uint64_t>(kBatchPolys, nPolys - p), parse_ns, derive_ns);
            std::free(staging);
        } else {
            want = nPolys * DEGREE1;
            const uint64_t chunk_bytes = kStreamChunkPolys * DEGREE1;
            uint8_t* staging = static_cast<uint8_t*>(std::aligned_alloc(4096, chunk_bytes));
            posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
            posix_fadvise(fd, 0, std::min<uint64_t>(chunk_bytes, want), POSIX_FADV_WILLNEED);
            for (uint64_t p = 0; p < nPolys; p += kStreamChunkPolys)
            {
                uint64_t cnt = std::min<uint64_t>(kStreamChunkPolys, nPolys - p);
                uint64_t off = p * DEGREE1;
                uint64_t next = off + cnt * DEGREE1;
                if (next < want)
                    posix_fadvise(fd, next, std::min<uint64_t>(chunk_bytes, want - next), POSIX_FADV_WILLNEED);
                auto r0 = std::chrono::steady_clock::now();
                got += readAll(fd, reinterpret_cast<char*>(staging), cnt * DEGREE1);
                read_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - r0).count();
                for (uint64_t q = 0; q < cnt; q += kBatchPolys)
                    encodeBatch(staging + q * DEGREE1, p + q, std::min<uint64_t>(kBatchPolys, cnt - q), parse_ns, derive_ns);
            }
            close(fd);
            std::free(staging);
        }
        auto end = std::chrono::steady_clock::now();
        last_load_us_ = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        if (got != want) std::cerr << "Short read from " << kDbPath << ": " << got << " of " << want << " bytes" << std::endl;

        std::cout << "│   ► Server Load Read Time: " << read_ns / 1000 << " us" << std::endl;
        std::cout << "│   ► Server Load Parse Time: " << parse_ns / 1000 << " us" << std::endl;
        std::cout << "│   ► Server Load Derive Time: " << derive_ns / 1000 << " us" << std::endl;
        std::cout << "│   ► Server Total Load Time: " << last_load_us_ << " us" << std::endl;
        // Backward-compatible name used by parse_via.py.
        std::cout << "│   ► Server Data Load Time: " << last_load_us_ << " us" << std::endl;
    }

	// Answer
	void Answer(RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel, RlweSampleQ2 result[8])
    {
        auto latency_start = std::chrono::high_resolution_clock::now();
        LoadDatabase();
        std::cout << "\u2502   \u25BA" << " Answer Time: " << std::endl;
        int64_t d[5];
        answerCore(c_rot, rsk, C_ctrl, C_sel, d);

        // ModSwitch
        std::cout << "\u2502     \u25C7" << " Final ModSwitch: ";
        auto start_FinalModSwitch = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < 8; i++)
        {
            nttQ2.ComputeInverse(SecondDimVec[i].rlweMask, SecondDimVec[i].rlweMask, 1, 1);
            modSwitch(MODULUS_Q2, MODULUS_Q3, SecondDimVec[i].rlweMask, result[i].rlweMask, DEGREE2);
            nttQ3.ComputeForward(result[i].rlweMask, result[i].rlweMask, 1, 1);
            nttQ2.ComputeInverse(SecondDimVec[i].rlweVal, SecondDimVec[i].rlweVal, 1, 1);
            modSwitch(MODULUS_Q2, MODULUS_Q4, SecondDimVec[i].rlweVal, result[i].rlweVal, DEGREE2);
        }
        auto end_FinalModSwitch = std::chrono::high_resolution_clock::now();
        auto duration_FinalModSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_FinalModSwitch - start_FinalModSwitch).count();
        std::cout << static_cast<double>(duration_FinalModSwitch) / 1000.0 << " ms" << std::endl;

        std::cout << "\u2502   \u25BA" << " Total Answer Time: " << static_cast<double>(d[0] + d[1] + d[2] + d[3] + d[4] + duration_FinalModSwitch) / 1000.0 << " ms" << std::endl;

        auto latency_end = std::chrono::high_resolution_clock::now();
        auto latency_duration = std::chrono::duration_cast<std::chrono::microseconds>(latency_end - latency_start).count();
        std::cout << "\u2502   \u25BA Server Response Latency: " << latency_duration << " us" << std::endl;
        std::cout << "\u2502   \u25BA Server Compute Time: " << (latency_duration - last_load_us_) << " us" << std::endl;
        std::cout << "\u2502 " << std::endl;
    }

	void Answer_Extraction(RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel, LweSampleQ2 result[8])
    {
        auto latency_start = std::chrono::high_resolution_clock::now();
        LoadDatabase();
        std::cout << "\u2502   \u25BA" << " Answer Time with blind extraction:: " << std::endl;
        int64_t d[5];
        answerCore(c_rot, rsk, C_ctrl, C_sel, d);

        // ModSwitch and extraction of the constant coefficient as an LWE sample.
        std::cout << "\u2502     \u25C7" << " Final ModSwitch and Extraction: ";
        auto start_ModSwitch_Extraction = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < 8; i++)
        {
            nttQ2.ComputeInverse(SecondDimVec[i].rlweMask, SecondDimVec[i].rlweMask, 1, 1);
            modSwitch(MODULUS_Q2, MODULUS_Q3, SecondDimVec[i].rlweMask, SecondDimVec[i].rlweMask, DEGREE2);
            result[i].lweMask[0] = SecondDimVec[i].rlweMask[0];
            for (int j = 1; j < DEGREE2; j++)
            {
                const uint64_t a = SecondDimVec[i].rlweMask[DEGREE2 - j];
                result[i].lweMask[j] = (a == 0) ? 0 : (MODULUS_Q3 - a);
            }
            nttQ2.ComputeInverse(SecondDimVec[i].rlweVal, SecondDimVec[i].rlweVal, 1, 1);
            modSwitch(MODULUS_Q2, MODULUS_Q4, SecondDimVec[i].rlweVal, &result[i].lweVal, 1);
        }
        auto end_ModSwitch_Extraction = std::chrono::high_resolution_clock::now();
        auto duration_ModSwitch_Extraction = std::chrono::duration_cast<std::chrono::microseconds>(end_ModSwitch_Extraction - start_ModSwitch_Extraction).count();
        std::cout << static_cast<double>(duration_ModSwitch_Extraction) / 1000.0 << " ms" << std::endl;

        std::cout << "\u2502   \u25BA" << " Total Answer Time with blind extraction: " << static_cast<double>(d[0] + d[1] + d[2] + d[3] + d[4] + duration_ModSwitch_Extraction) / 1000.0 << " ms" << std::endl;

        auto latency_end = std::chrono::high_resolution_clock::now();
        auto latency_duration = std::chrono::duration_cast<std::chrono::microseconds>(latency_end - latency_start).count();
        std::cout << "\u2502   \u25BA Server Response Latency: " << latency_duration << " us" << std::endl;
        std::cout << "\u2502   \u25BA Server Compute Time: " << (latency_duration - last_load_us_) << " us" << std::endl;
        std::cout << "\u2502 " << std::endl;
    }

    // DMux -> ModSwitch -> First Dimension -> Ring Switch -> CMux, shared by
    // both answer variants. Leaves the 8 surviving RLWE(Q2, 512) ciphertexts in
    // SecondDimVec[0..7] (NTT form). dur[] = per-phase microseconds.
    void answerCore(RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1>& C_ctrl, std::vector<RgswSampleQ2>& C_sel, int64_t dur[5])
    {
        // Dmux
        std::cout << "\u2502     \u25C7" << " Dmux: ";
        auto start_DMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < 8; i++)
        {
            copy_Deg1(c_rot[i].rlweMask1, FirstDimVec[i].rlweMask1);
            copy_Deg1(c_rot[i].rlweVal1, FirstDimVec[i].rlweVal1);
            copy_Deg1(c_rot[i].rlweMask2, FirstDimVec[i].rlweMask2);
            copy_Deg1(c_rot[i].rlweVal2, FirstDimVec[i].rlweVal2);
        }
        for (uint64_t i = 3; i < LOG_ROW; i++)
        {
            for (int64_t j = (1 << i) - 1; j >= 0; j--)
            {
                DMuxQ1(CEIL_LOG_Q1, GADGET1_L, GADGET1_L, LOG_GADGET1_BASE, LOG_GADGET1_BASE, C_ctrl[i - 3], FirstDimVec[j], FirstDimVec[j << 1], FirstDimVec[(j << 1) + 1]);
            }
        }
        auto end_DMux = std::chrono::high_resolution_clock::now();
        dur[0] = std::chrono::duration_cast<std::chrono::microseconds>(end_DMux - start_DMux).count();
        std::cout << static_cast<double>(dur[0]) / 1000.0 << " ms" << std::endl;

        // ModSwitch Q1 -> Q2 (both CRT limbs -> one Q2 limb, NTT form, degree 2048)
        std::cout << "\u2502     \u25C7" << " ModSwitch: ";
        auto start_ModSwitch = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < ROW; i++)
        {
            nttQ1_1.ComputeInverse(FirstDimVec[i].rlweMask1, FirstDimVec[i].rlweMask1, 1, 1);
            nttQ1_2.ComputeInverse(FirstDimVec[i].rlweMask2, FirstDimVec[i].rlweMask2, 1, 1);
            CRT(temp_message1, FirstDimVec[i].rlweMask1, FirstDimVec[i].rlweMask2, DEGREE1);
            modSwitch(MODULUS_Q1, MODULUS_Q2, temp_message1, FirstDimVec[i].rlweMask1, DEGREE1);
            nttQ2Deg1.ComputeForward(FirstDimVec[i].rlweMask1, FirstDimVec[i].rlweMask1, 1, 1);
            nttQ1_1.ComputeInverse(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_2.ComputeInverse(FirstDimVec[i].rlweVal2, FirstDimVec[i].rlweVal2, 1, 1);
            CRT(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1, MODULUS_Q2, temp_message1, FirstDimVec[i].rlweVal1, DEGREE1);
            nttQ2Deg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
        }
        auto end_ModSwitch = std::chrono::high_resolution_clock::now();
        dur[1] = std::chrono::duration_cast<std::chrono::microseconds>(end_ModSwitch - start_ModSwitch).count();
        std::cout << static_cast<double>(dur[1]) / 1000.0 << " ms" << std::endl;

        // First Dimension
        std::cout << "\u2502     \u25C7" << " First Dimension: ";
        auto start_FirstDim = std::chrono::high_resolution_clock::now();
        uint64_t fdDbSize = numSlots();
        for (uint64_t j = 0; j < COL; j++)
        {
            intel::hexl::EltwiseMultMod(FirstDimResMasks[j], FirstDimVec[0].rlweMask1, database[j % fdDbSize], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseMultMod(FirstDimResVals[j], FirstDimVec[0].rlweVal1, database[j % fdDbSize], DEGREE1, MODULUS_Q2, 1);
            for (uint64_t i = 1; i < ROW; i++)
            {
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[i].rlweMask1, database[(i * COL + j) % fdDbSize], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseAddMod(FirstDimResMasks[j], FirstDimResMasks[j], temp_prod1, DEGREE1, MODULUS_Q2);
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[i].rlweVal1, database[(i * COL + j) % fdDbSize], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseAddMod(FirstDimResVals[j], FirstDimResVals[j], temp_prod1, DEGREE1, MODULUS_Q2);
            }
        }
        auto end_FirstDim = std::chrono::high_resolution_clock::now();
        dur[2] = std::chrono::duration_cast<std::chrono::microseconds>(end_FirstDim - start_FirstDim).count();
        std::cout << static_cast<double>(dur[2]) / 1000.0 << " ms" << std::endl;

        // Ring switch: key-switch (a, b) from s1 to s2(X^4) over (2048, Q2),
        // then project onto coefficients 0 mod 4 (degree 512).
        std::cout << "\u2502     \u25C7" << " Ring Switch: ";
        auto start_RingSwitch = std::chrono::high_resolution_clock::now();
        std::fill(temp_message1, temp_message1 + DEGREE1, 0);
        for (uint64_t i = 0; i < COL; i++)
        {
            nttQ2Deg1.ComputeInverse(FirstDimResMasks[i], FirstDimResMasks[i], 1, 1);
            signedDecompose(FirstDimResMasks[i], DEGREE1, MODULUS_Q2, CEIL_LOG_Q2, GADGET_RSK_L, LOG_GADGET_RSK_BASE, gadget_poly1);
            for (int j = 0; j < GADGET_RSK_L; j++)
            {
                digitsToMod(gadget_poly1[j], gadget_poly1[j], DEGREE1, MODULUS_Q2);
                nttQ2Deg1.ComputeForward(gadget_poly1[j], gadget_poly1[j], 1, 1);
            }
            intel::hexl::EltwiseMultMod(FirstDimResMasks[i], gadget_poly1[0], rsk.rlweMasks[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[0], rsk.rlweVals[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseSubMod(FirstDimResVals[i], FirstDimResVals[i], temp_prod1, DEGREE1, MODULUS_Q2);

            for (int j = 1; j < GADGET_RSK_L; j++)
            {
                intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[j], rsk.rlweMasks[j], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseAddMod(FirstDimResMasks[i], FirstDimResMasks[i], temp_prod1, DEGREE1, MODULUS_Q2);
                intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[j], rsk.rlweVals[j], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseSubMod(FirstDimResVals[i], FirstDimResVals[i], temp_prod1, DEGREE1, MODULUS_Q2);
            }
            intel::hexl::EltwiseSubMod(FirstDimResMasks[i], temp_message1, FirstDimResMasks[i], DEGREE1, MODULUS_Q2);
            nttQ2Deg1.ComputeInverse(FirstDimResMasks[i], FirstDimResMasks[i], 1, 1);
            ring_project(0, FirstDimResMasks[i], SecondDimVec[i].rlweMask);
            nttQ2.ComputeForward(SecondDimVec[i].rlweMask, SecondDimVec[i].rlweMask, 1, 1);
            nttQ2Deg1.ComputeInverse(FirstDimResVals[i], FirstDimResVals[i], 1, 1);
            ring_project(0, FirstDimResVals[i], SecondDimVec[i].rlweVal);
            nttQ2.ComputeForward(SecondDimVec[i].rlweVal, SecondDimVec[i].rlweVal, 1, 1);
        }
        auto end_RingSwitch = std::chrono::high_resolution_clock::now();
        dur[3] = std::chrono::duration_cast<std::chrono::microseconds>(end_RingSwitch - start_RingSwitch).count();
        std::cout << static_cast<double>(dur[3]) / 1000.0 << " ms" << std::endl;

        // CMux
        std::cout << "\u2502     \u25C7" << " CMux: ";
        auto start_CMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = LOG_COL; i > 3; i--)
        {
            for (uint64_t j = 0; j < (UINT64_C(1) << (i - 1)); j++)
            {
                CMuxQ2(CEIL_LOG_Q2, GADGET2_L_1, GADGET2_L_2, LOG_GADGET2_BASE_1, LOG_GADGET2_BASE_2, C_sel[LOG_COL - i], SecondDimVec[j << 1], SecondDimVec[(j << 1) + 1], SecondDimVec[j]);
            }
        }
        auto end_CMux = std::chrono::high_resolution_clock::now();
        dur[4] = std::chrono::duration_cast<std::chrono::microseconds>(end_CMux - start_CMux).count();
        std::cout << static_cast<double>(dur[4]) / 1000.0 << " ms" << std::endl;
    }

    // Balanced signed gadget decomposition of `poly` (coefficient form,
    // values in [0, q), q < 2^logq): the centred value x is approximated by
    // sum_i out[i] * 2^(logq - L*logB + i*logB). Digits are int64 stored in
    // the uint64 buffers (two's complement). (The original power-of-two
    // decomposition did not centre mod q, which is wrong for q != 2^logq.)
    void signedDecompose(const uint64_t* poly, uint64_t n, uint64_t q, uint64_t logq, uint64_t L, uint64_t logB, std::vector<uint64_t*>& out)
    {
        const int64_t half_q = static_cast<int64_t>(q >> 1);
        const int64_t sq = static_cast<int64_t>(q);
        const uint64_t shift = logq - L * logB;
        const int64_t round = INT64_C(1) << (shift - 1);
        const int64_t halfB = INT64_C(1) << (logB - 1);
        const int64_t maskB = (INT64_C(1) << logB) - 1;
        int64_t* y = reinterpret_cast<int64_t*>(decomp_scratch_);
        for (uint64_t j = 0; j < n; j++)
        {
            int64_t x = static_cast<int64_t>(poly[j]);
            x = (x > half_q) ? x - sq : x;
            y[j] = (x + round) >> shift;
        }
        for (uint64_t i = 0; i + 1 < L; i++)
        {
            int64_t* o = reinterpret_cast<int64_t*>(out[i]);
            for (uint64_t j = 0; j < n; j++)
            {
                int64_t dgt = ((y[j] + halfB) & maskB) - halfB;
                o[j] = dgt;
                y[j] = (y[j] - dgt) >> logB;
            }
        }
        std::memcpy(out[L - 1], y, n * sizeof(uint64_t));
    }

    // Map signed digits to [0, p).
    static void digitsToMod(const uint64_t* in, uint64_t* out, uint64_t n, uint64_t p)
    {
        const int64_t* s = reinterpret_cast<const int64_t*>(in);
        for (uint64_t j = 0; j < n; j++)
            out[j] = static_cast<uint64_t>(s[j]) + (p & static_cast<uint64_t>(s[j] >> 63));
    }

    // half extetnal product for Q1
    void halfExternalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ1& rlevSample, uint64_t* poly, RlweSampleQ1& result)
    {
        signedDecompose(poly, DEGREE1, MODULUS_Q1, log_modulus, gadget_l, log_gadget_base, gadget_poly1);
        for (int i = 0; i < gadget_l; i++)
        {
            digitsToMod(gadget_poly1[i], gadget_poly2[i], DEGREE1, MODULUS_Q1_2);
            digitsToMod(gadget_poly1[i], gadget_poly1[i], DEGREE1, MODULUS_Q1_1);
            nttQ1_2.ComputeForward(gadget_poly2[i], gadget_poly2[i], 1, 1);
            nttQ1_1.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
        }
        intel::hexl::EltwiseMultMod(result.rlweMask1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweMask1, DEGREE1, MODULUS_Q1_1, 1);
        intel::hexl::EltwiseMultMod(result.rlweVal1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweVal1, DEGREE1, MODULUS_Q1_1, 1);

        intel::hexl::EltwiseMultMod(result.rlweMask2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweMask2, DEGREE1, MODULUS_Q1_2, 1);
        intel::hexl::EltwiseMultMod(result.rlweVal2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweVal2, DEGREE1, MODULUS_Q1_2, 1);

        for (int i = 1; i < gadget_l; i++)
        {
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweMask1, DEGREE1, MODULUS_Q1_1, 1);
            intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_prod1, DEGREE1, MODULUS_Q1_1);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweVal1, DEGREE1, MODULUS_Q1_1, 1);
            intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_prod1, DEGREE1, MODULUS_Q1_1);

            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweMask2, DEGREE1, MODULUS_Q1_2, 1);
            intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_prod1, DEGREE1, MODULUS_Q1_2);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweVal2, DEGREE1, MODULUS_Q1_2, 1);
            intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_prod1, DEGREE1, MODULUS_Q1_2);
        }
    }

    // external product for Q1 (result must not alias temp_rlweSampleQ1_b)
    void externalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result)
    {
        nttQ1_1.ComputeInverse(temp_message1_1, rlweSample.rlweMask1, 1, 1);
        nttQ1_2.ComputeInverse(temp_message1_2, rlweSample.rlweMask2, 1, 1);
        CRT(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
        halfExternalProdQ1(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1, result);
        nttQ1_1.ComputeInverse(temp_message1_1, rlweSample.rlweVal1, 1, 1);
        nttQ1_2.ComputeInverse(temp_message1_2, rlweSample.rlweVal2, 1, 1);
        CRT(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
        halfExternalProdQ1(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1, temp_rlweSampleQ1_b);
        intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_rlweSampleQ1_b.rlweMask1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_rlweSampleQ1_b.rlweVal1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_rlweSampleQ1_b.rlweMask2, DEGREE1, MODULUS_Q1_2);
        intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_rlweSampleQ1_b.rlweVal2, DEGREE1, MODULUS_Q1_2);
    }

    // DMux for Q1: result1 = C (x) ct = b*ct, result0 = ct - result1.
    // result0 may alias rlweSample; result1 must not.
    void DMuxQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result0, RlweSampleQ1& result1)
    {
        externalProdQ1(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, rlweSample, temp_rlweSampleQ1);
        copy_Deg1(temp_rlweSampleQ1.rlweMask1, result1.rlweMask1);
        copy_Deg1(temp_rlweSampleQ1.rlweVal1, result1.rlweVal1);
        copy_Deg1(temp_rlweSampleQ1.rlweMask2, result1.rlweMask2);
        copy_Deg1(temp_rlweSampleQ1.rlweVal2, result1.rlweVal2);
        intel::hexl::EltwiseSubMod(result0.rlweMask1, rlweSample.rlweMask1, temp_rlweSampleQ1.rlweMask1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseSubMod(result0.rlweVal1, rlweSample.rlweVal1, temp_rlweSampleQ1.rlweVal1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseSubMod(result0.rlweMask2, rlweSample.rlweMask2, temp_rlweSampleQ1.rlweMask2, DEGREE1, MODULUS_Q1_2);
        intel::hexl::EltwiseSubMod(result0.rlweVal2, rlweSample.rlweVal2, temp_rlweSampleQ1.rlweVal2, DEGREE1, MODULUS_Q1_2);
    }

    // half product for Q2
    void halfExternalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ2& rlevSample, uint64_t* poly, RlweSampleQ2& result)
    {
        signedDecompose(poly, DEGREE2, MODULUS_Q2, log_modulus, gadget_l, log_gadget_base, gadget_poly1);
        for (int i = 0; i < gadget_l; i++)
        {
            digitsToMod(gadget_poly1[i], gadget_poly1[i], DEGREE2, MODULUS_Q2);
            nttQ2.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
        }

        intel::hexl::EltwiseMultMod(result.rlweMask, gadget_poly1[0], rlevSample.rlweSamples[0].rlweMask, DEGREE2, MODULUS_Q2, 1);
        intel::hexl::EltwiseMultMod(result.rlweVal, gadget_poly1[0], rlevSample.rlweSamples[0].rlweVal, DEGREE2, MODULUS_Q2, 1);

        for (int i = 1; i < gadget_l; i++)
        {
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweMask, DEGREE2, MODULUS_Q2, 1);
            intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, temp_prod1, DEGREE2, MODULUS_Q2);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweVal, DEGREE2, MODULUS_Q2, 1);
            intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, temp_prod1, DEGREE2, MODULUS_Q2);
        }
    }

    // external product for Q2 (result must not alias temp_rlweSampleQ2_b)
    void externalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample, RlweSampleQ2& result)
    {
        nttQ2.ComputeInverse(temp_message1_1, rlweSample.rlweMask, 1, 1);
        halfExternalProdQ2(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1_1, result);
        nttQ2.ComputeInverse(temp_message1_1, rlweSample.rlweVal, 1, 1);
        halfExternalProdQ2(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1_1, temp_rlweSampleQ2_b);
        intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, temp_rlweSampleQ2_b.rlweMask, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, temp_rlweSampleQ2_b.rlweVal, DEGREE2, MODULUS_Q2);
    }

    // CMux for Q2: result = ct0 + C (x) (ct1 - ct0). result may alias ct0.
    void CMuxQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample0, RlweSampleQ2& rlweSample1, RlweSampleQ2& result)
    {
        intel::hexl::EltwiseSubMod(temp_rlweSampleQ2.rlweMask, rlweSample1.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseSubMod(temp_rlweSampleQ2.rlweVal, rlweSample1.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2);
        externalProdQ2(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, temp_rlweSampleQ2, temp_rlweSampleQ2_c);
        intel::hexl::EltwiseAddMod(result.rlweMask, temp_rlweSampleQ2_c.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseAddMod(result.rlweVal, temp_rlweSampleQ2_c.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2);
    }

	std::vector<uint64_t*> database;
    std::vector<RlweSampleQ1> FirstDimVec;
    std::vector<uint64_t*> FirstDimResMasks;
    std::vector<uint64_t*> FirstDimResVals;
    std::vector<RlweSampleQ2> SecondDimVec;
    intel::hexl::NTT nttQ1_1;
    intel::hexl::NTT nttQ1_2;
    intel::hexl::NTT nttQ2;
    intel::hexl::NTT nttQ2Deg1;
    intel::hexl::NTT nttQ3;
    uint64_t* temp_prod1;
    uint64_t* temp_message1;
    uint64_t* temp_message1_1;
    uint64_t* temp_message1_2;
    uint64_t* decomp_scratch_;
    std::vector<uint64_t*> gadget_poly1;
    std::vector<uint64_t*> gadget_poly2;
    RlweSampleQ1 temp_rlweSampleQ1;
    RlweSampleQ1 temp_rlweSampleQ1_b;
    RlweSampleQ2 temp_rlweSampleQ2;
    RlweSampleQ2 temp_rlweSampleQ2_b;
    RlweSampleQ2 temp_rlweSampleQ2_c;
};

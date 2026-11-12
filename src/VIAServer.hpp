#pragma once
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include <cstring>
#include <cstdio>
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
        nttQ2Deg1.ComputeInverse(database[idx_row_col], out, 1, 1);
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

    // Untimed (called before the response-latency timer starts): evicts the
    // server-data file from the page cache, then reports how many of its pages
    // are still resident (mmap + mincore), proving the following load is cold.
    // LoadDatabase() keeps its own eviction, so its timed path is unchanged.
    void EvictAndReportResidencyBeforeLoad()
    {
        if (!write_to_file_) return;
        int fd = open(kDbPath, O_RDONLY);
        if (fd == -1) {
            std::cerr << "Failed to open " << kDbPath << " for the residency check!" << std::endl;
            return;
        }
        posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
        const uint64_t page = static_cast<uint64_t>(sysconf(_SC_PAGESIZE));
        const uint64_t len = static_cast<uint64_t>(lseek(fd, 0, SEEK_END));
        uint64_t resident = 0, total = 0;
        if (len > 0) {
            total = (len + page - 1) / page;
            void* addr = mmap(nullptr, len, PROT_READ, MAP_SHARED, fd, 0);
            if (addr == MAP_FAILED) {
                std::cerr << "mmap failed for " << kDbPath << std::endl;
            } else {
                std::vector<unsigned char> vec(total);
                if (mincore(addr, len, vec.data()) != 0)
                    std::cerr << "mincore failed for " << kDbPath << std::endl;
                for (unsigned char b : vec) resident += (b & 1);
                munmap(addr, len);
            }
        }
        close(fd);
        char pct[32];
        std::snprintf(pct, sizeof(pct), "%.3f", total ? 100.0 * resident / total : 0.0);
        std::cout << "Server Page Cache Residency Before Load: " << resident << "/" << total
                  << " pages (" << pct << "%)" << std::endl;
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
	void Answer(RlweSampleQ1 c_rot[4], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel, RlweSampleQ2 result[8])
    {
        EvictAndReportResidencyBeforeLoad();
        auto latency_start = std::chrono::high_resolution_clock::now();
        LoadDatabase();
        std::cout << "\u2502   \u25BA" << " Answer Time: " << std::endl;
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
        auto duration_DMux = std::chrono::duration_cast<std::chrono::microseconds>(end_DMux - start_DMux).count();
        std::cout << static_cast<double>(duration_DMux) / 1000.0 << " ms" << std::endl;


        // ModSwitch
        std::cout << "\u2502     \u25C7" << " ModSwitch: ";
        auto start_ModSwitch = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < ROW; i++)
        {
            nttQ1_1.ComputeInverse(FirstDimVec[i].rlweMask1, FirstDimVec[i].rlweMask1, 1, 1);
            nttQ1_2.ComputeInverse(FirstDimVec[i].rlweMask2, FirstDimVec[i].rlweMask2, 1, 1);
            CRT(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1, MODULUS_Q2, temp_message1, FirstDimVec[i].rlweMask1, DEGREE1);
            nttQ2Deg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_1.ComputeInverse(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_2.ComputeInverse(FirstDimVec[i].rlweVal2, FirstDimVec[i].rlweVal2, 1, 1);
            CRT(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1, MODULUS_Q2, temp_message1, FirstDimVec[i].rlweVal1, DEGREE1);
            nttQ2Deg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
        }
        auto end_ModSwitch = std::chrono::high_resolution_clock::now();
        auto duration_ModSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_ModSwitch - start_ModSwitch).count();
        std::cout << static_cast<double>(duration_ModSwitch) / 1000.0 << " ms" << std::endl;

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
        auto duration_FirstDim = std::chrono::duration_cast<std::chrono::microseconds>(end_FirstDim - start_FirstDim).count();
        std::cout << static_cast<double>(duration_FirstDim) / 1000.0 << " ms" << std::endl;

        // ring swith
        std::cout << "\u2502     \u25C7" << " Ring Switch: ";
        auto start_RingSwitch = std::chrono::high_resolution_clock::now();
        std::fill(temp_message1, temp_message1 + DEGREE1, 0);
        for (uint64_t i = 0; i < COL; i++)
        {
            nttQ2Deg1.ComputeInverse(FirstDimResMasks[i], SecondDimVec[i].rlweMask, 1, 1);
            GadgetDecompositionDeg1(FirstDimResMasks[i], gadget_poly1, CEIL_LOG_Q2, GADGET_RSK_L, GADGET_RSK_BASE);
            for (int j = 0; j < GADGET_RSK_L; j++)
            {
                nttQ2.ComputeForward(gadget_poly1[j], gadget_poly1[j], 1, 1);
            }
            intel::hexl::EltwiseMultMod(FirstDimResMasks[i], gadget_poly1[0], rsk.rlweMasks[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[0], rsk.rlweVals[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseSubMod(FirstDimResVals[i], SecondDimVec[i].rlweVal, temp_prod1, DEGREE1, MODULUS_Q2);

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
            ring_project(0, SecondDimVec[i].rlweVal, SecondDimVec[i].rlweVal);
            nttQ2.ComputeForward(SecondDimVec[i].rlweVal, SecondDimVec[i].rlweVal, 1, 1);
        }
        auto end_RingSwitch = std::chrono::high_resolution_clock::now();
        auto duration_RingSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_RingSwitch - start_RingSwitch).count();
        std::cout << static_cast<double>(duration_RingSwitch) / 1000.0 << " ms" << std::endl;

        // CMux
        std::cout << "\u2502     \u25C7" << " CMux: ";
        auto start_CMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = LOG_COL; i > 3; i--)
        {
            for (uint64_t j = 0; j < (1 << (i - 1)); j++)
            {
                CMuxQ2(CEIL_LOG_Q2, GADGET2_L_1, GADGET2_L_2, LOG_GADGET2_BASE_1, LOG_GADGET2_BASE_2, C_sel[LOG_COL - i], SecondDimVec[j << 1], SecondDimVec[(j << 1) + 1], SecondDimVec[j]);
            }
        }
        auto end_CMux = std::chrono::high_resolution_clock::now();
        auto duration_CMux = std::chrono::duration_cast<std::chrono::microseconds>(end_CMux - start_CMux).count();
        std::cout << static_cast<double>(duration_CMux) / 1000.0 << " ms" << std::endl;



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

        std::cout << "\u2502   \u25BA" << " Total Answer Time: " << static_cast<double>(duration_DMux + duration_ModSwitch + duration_FirstDim + duration_RingSwitch + duration_CMux + duration_FinalModSwitch) / 1000.0 << " ms" << std::endl;
        
        auto latency_end = std::chrono::high_resolution_clock::now();
        auto latency_duration = std::chrono::duration_cast<std::chrono::microseconds>(latency_end - latency_start).count();
        std::cout << "\u2502   \u25BA Server Response Latency: " << latency_duration << " us" << std::endl;
        std::cout << "\u2502   \u25BA Server Compute Time: " << (latency_duration - last_load_us_) << " us" << std::endl;
        std::cout << "\u2502 " << std::endl;
    }

	void Answer_Extraction(RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel, LweSampleQ2 result[8])
    {
        EvictAndReportResidencyBeforeLoad();
        auto latency_start = std::chrono::high_resolution_clock::now();
        LoadDatabase();
        std::cout << "\u2502   \u25BA" << " Answer Time with blind extraction:: " << std::endl;
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
        auto duration_DMux = std::chrono::duration_cast<std::chrono::microseconds>(end_DMux - start_DMux).count();
        std::cout << static_cast<double>(duration_DMux) / 1000.0 << " ms" << std::endl;


        // ModSwitch
        std::cout << "\u2502     \u25C7" << " ModSwitch: ";
        auto start_ModSwitch = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < ROW; i++)
        {
            nttQ1_1.ComputeInverse(FirstDimVec[i].rlweMask1, FirstDimVec[i].rlweMask1, 1, 1);
            nttQ1_2.ComputeInverse(FirstDimVec[i].rlweMask2, FirstDimVec[i].rlweMask2, 1, 1);
            CRT(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1, MODULUS_Q2, temp_message1, FirstDimVec[i].rlweMask1, DEGREE1);
            nttQ2Deg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_1.ComputeInverse(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_2.ComputeInverse(FirstDimVec[i].rlweVal2, FirstDimVec[i].rlweVal2, 1, 1);
            CRT(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1, MODULUS_Q2, temp_message1, FirstDimVec[i].rlweVal1, DEGREE1);
            nttQ2Deg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
        }
        auto end_ModSwitch = std::chrono::high_resolution_clock::now();
        auto duration_ModSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_ModSwitch - start_ModSwitch).count();
        std::cout << static_cast<double>(duration_ModSwitch) / 1000.0 << " ms" << std::endl;


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
        auto duration_FirstDim = std::chrono::duration_cast<std::chrono::microseconds>(end_FirstDim - start_FirstDim).count();
        std::cout << static_cast<double>(duration_FirstDim) / 1000.0 << " ms" << std::endl;

        // ring swith
        std::cout << "\u2502     \u25C7" << " Ring Switch: ";
        auto start_RingSwitch = std::chrono::high_resolution_clock::now();
        std::fill(temp_message1, temp_message1 + DEGREE1, 0);
        for (uint64_t i = 0; i < COL; i++)
        {
            nttQ2Deg1.ComputeInverse(FirstDimResMasks[i], SecondDimVec[i].rlweMask, 1, 1);
            GadgetDecompositionDeg1(FirstDimResMasks[i], gadget_poly1, CEIL_LOG_Q2, GADGET_RSK_L, GADGET_RSK_BASE);
            for (int j = 0; j < GADGET_RSK_L; j++)
            {
                nttQ2.ComputeForward(gadget_poly1[j], gadget_poly1[j], 1, 1);
            }
            intel::hexl::EltwiseMultMod(FirstDimResMasks[i], gadget_poly1[0], rsk.rlweMasks[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[0], rsk.rlweVals[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseSubMod(FirstDimResVals[i], SecondDimVec[i].rlweVal, temp_prod1, DEGREE1, MODULUS_Q2);

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
        auto duration_RingSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_RingSwitch - start_RingSwitch).count();
        std::cout << static_cast<double>(duration_RingSwitch) / 1000.0 << " ms" << std::endl;

        // CMux
        std::cout << "\u2502     \u25C7" << " CMux: ";
        auto start_CMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = LOG_COL; i > 3; i--)
        {
            for (uint64_t j = 0; j < (1 << (i - 1)); j++)
            {
                CMuxQ2(CEIL_LOG_Q2, GADGET2_L_1, GADGET2_L_2, LOG_GADGET2_BASE_1, LOG_GADGET2_BASE_2, C_sel[LOG_COL - i], SecondDimVec[j << 1], SecondDimVec[(j << 1) + 1], SecondDimVec[j]);
            }
        }
        auto end_CMux = std::chrono::high_resolution_clock::now();
        auto duration_CMux = std::chrono::duration_cast<std::chrono::microseconds>(end_CMux - start_CMux).count();
        std::cout << static_cast<double>(duration_CMux) / 1000.0 << " ms" << std::endl;

        // ModSwitch and extration
        std::cout << "\u2502     \u25C7" << " Final ModSwitch and Extraction: ";
        auto start_ModSwitch_Extraction = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < 8; i++)
        {
            nttQ2.ComputeInverse(SecondDimVec[i].rlweMask, SecondDimVec[i].rlweMask, 1, 1);
            modSwitch(MODULUS_Q2, MODULUS_Q2, SecondDimVec[i].rlweMask, SecondDimVec[i].rlweMask, DEGREE2);
            result[i].lweMask[0] = SecondDimVec[i].rlweMask[0];
            for (int j = 1; j < DEGREE2; j++)
            {
                result[i].lweMask[j] = -SecondDimVec[i].rlweMask[DEGREE2-j];
            }
            nttQ2.ComputeInverse(SecondDimVec[i].rlweVal, SecondDimVec[i].rlweVal, 1, 1);
            modSwitch(MODULUS_Q2, MODULUS_Q4, SecondDimVec[i].rlweVal, &result[i].lweVal, 1);
        }
        auto end_ModSwitch_Extraction = std::chrono::high_resolution_clock::now();
        auto duration_ModSwitch_Extraction = std::chrono::duration_cast<std::chrono::microseconds>(end_ModSwitch_Extraction - start_ModSwitch_Extraction).count();
        std::cout << static_cast<double>(duration_ModSwitch_Extraction) / 1000.0 << " ms" << std::endl;

        std::cout << "\u2502   \u25BA" << " Total Answer Time with blind extraction: " << static_cast<double>(duration_DMux + duration_ModSwitch + duration_FirstDim + duration_RingSwitch + duration_CMux + duration_ModSwitch_Extraction) / 1000.0 << " ms" << std::endl;
        
        auto latency_end = std::chrono::high_resolution_clock::now();
        auto latency_duration = std::chrono::duration_cast<std::chrono::microseconds>(latency_end - latency_start).count();
        std::cout << "\u2502   \u25BA Server Response Latency: " << latency_duration << " us" << std::endl;
        std::cout << "\u2502   \u25BA Server Compute Time: " << (latency_duration - last_load_us_) << " us" << std::endl;
        std::cout << "\u2502 " << std::endl;
    }

    // half extetnal product for Q1
    void halfExternalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ1& rlevSample, uint64_t* poly, RlweSampleQ1& result)
    {
        GadgetDecompositionDeg1(poly, gadget_poly1, log_modulus, gadget_l, log_gadget_base);
        for (int i = 0; i < gadget_l; i++)
        {
            nttQ1_2.ComputeForward(gadget_poly2[i], gadget_poly1[i], 1, 1);
            nttQ1_1.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
        }
        intel::hexl::EltwiseMultMod(result.rlweMask1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweMask1, DEGREE1, MODULUS_Q1_1, 1);
        intel::hexl::EltwiseMultMod(result.rlweVal1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweVal1, DEGREE1, MODULUS_Q1_1, 1);

        intel::hexl::EltwiseMultMod(result.rlweMask2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweMask2, DEGREE1, MODULUS_Q1_1, 2);
        intel::hexl::EltwiseMultMod(result.rlweVal2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweVal2, DEGREE1, MODULUS_Q1_1, 2);

        for (int i = 1; i < gadget_l; i++)
        {
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweMask1, DEGREE1, MODULUS_Q1_1, 1);
            intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_prod1, DEGREE1, MODULUS_Q1_1);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweVal1, DEGREE1, MODULUS_Q1_1, 1);
            intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_prod1, DEGREE1, MODULUS_Q1_1);

            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweMask2, DEGREE1, MODULUS_Q1_1, 2);
            intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_prod1, DEGREE1, MODULUS_Q1_2);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweVal2, DEGREE1, MODULUS_Q1_1, 2);
            intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_prod1, DEGREE1, MODULUS_Q1_2);
        }
    }

    // external product for Q1
    void externalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result)
    {
        nttQ1_1.ComputeInverse(temp_message1_1, rlweSample.rlweMask1, 1, 1);
        nttQ1_2.ComputeInverse(temp_message1_2, rlweSample.rlweMask2, 1, 1);
        CRT(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
        halfExternalProdQ1(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1, result);
        nttQ1_1.ComputeInverse(temp_message1_1, rlweSample.rlweVal1, 1, 1);
        nttQ1_2.ComputeInverse(temp_message1_2, rlweSample.rlweVal2, 1, 1);
        CRT(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
        halfExternalProdQ1(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1, temp_rlweSampleQ1);
        intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_rlweSampleQ1.rlweMask1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_rlweSampleQ1.rlweVal1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_rlweSampleQ1.rlweMask2, DEGREE1, MODULUS_Q1_2);
        intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_rlweSampleQ1.rlweVal2, DEGREE1, MODULUS_Q1_2);
    }

    // DMux for Q1
    void DMuxQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result0, RlweSampleQ1& result1)
    {
        externalProdQ1(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, rlweSample, temp_rlweSampleQ1);
        intel::hexl::EltwiseSubMod(result0.rlweMask1, rlweSample.rlweMask1, temp_rlweSampleQ1.rlweMask1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseSubMod(result0.rlweVal1, rlweSample.rlweVal1, temp_rlweSampleQ1.rlweVal1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseSubMod(result0.rlweMask2, rlweSample.rlweMask2, temp_rlweSampleQ1.rlweMask2, DEGREE1, MODULUS_Q1_2);
        intel::hexl::EltwiseSubMod(result0.rlweVal2, rlweSample.rlweVal2, temp_rlweSampleQ1.rlweVal2, DEGREE1, MODULUS_Q1_2);

        intel::hexl::EltwiseSubMod(result1.rlweMask1, result0.rlweMask1, rlweSample.rlweMask1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseSubMod(result1.rlweVal1, result0.rlweVal1, rlweSample.rlweVal1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseSubMod(result1.rlweMask2, result0.rlweMask2, rlweSample.rlweMask2, DEGREE1, MODULUS_Q1_2);    
        intel::hexl::EltwiseSubMod(result1.rlweVal2, result0.rlweVal2, rlweSample.rlweVal2, DEGREE1, MODULUS_Q1_2);
    }

    // half product for Q2
    void halfExternalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ2& rlevSample, uint64_t* poly, RlweSampleQ2& result)
    {
        GadgetDecompositionDeg2(poly, gadget_poly1, log_modulus, gadget_l, log_gadget_base);
        for (int i = 0; i < gadget_l; i++)
        {
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

    // external product for Q2
    void externalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample, RlweSampleQ2& result)
    {
        nttQ2.ComputeInverse(temp_message1_1, rlweSample.rlweMask, 1, 1);
        halfExternalProdQ2(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1_1, result);
        nttQ2.ComputeInverse(temp_message1_1, rlweSample.rlweVal, 1, 1);
        halfExternalProdQ2(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1_1, temp_rlweSampleQ2);
        intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, temp_rlweSampleQ2.rlweMask, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, temp_rlweSampleQ2.rlweVal, DEGREE2, MODULUS_Q2);
    }

    // CMux for Q2
    void CMuxQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample0, RlweSampleQ2& rlweSample1, RlweSampleQ2& result)
    {
        intel::hexl::EltwiseSubMod(temp_rlweSampleQ2.rlweMask, rlweSample1.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseSubMod(temp_rlweSampleQ2.rlweVal, rlweSample1.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2);
        externalProdQ2(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, temp_rlweSampleQ2, result);
        intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2);
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
    std::vector<uint64_t*> gadget_poly1;
    std::vector<uint64_t*> gadget_poly2;
    RlweSampleQ1 temp_rlweSampleQ1;
    RlweSampleQ2 temp_rlweSampleQ2;
};


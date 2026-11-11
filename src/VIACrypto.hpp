#pragma once
// Self-contained, functionally correct client-side crypto for VIA (non-C).
// Replaces the use of FHE_LWE by VIAClient, whose RLWE/RGSW/rsk generation
// had several bugs (see VIA_VERIFY_FIX notes). FHE_LWE is left untouched for
// VIA-C. Conventions (shared with VIAServer):
//   * RLWE: (a, b) with b = a*s + e + m; decryption m = b - a*s.
//   * Q1 = Q1_1 * Q1_2 ciphertexts are kept as two NTT-form CRT limbs.
//   * RLev_g(mu)[i] = RLWE(mu * 2^(logq - L*logB + i*logB)).
//   * RGSW(mu) = (keyRlev = RLev(-mu*s), msgRlev = RLev(mu)), so that
//     ext(C, (a, b)) = half(keyRlev, a) + half(msgRlev, b) = RLWE(mu*(b - a*s)).
//   * Ring-switching key: rsk[j] = RLWE_{s2(X^4)}(s1 * g_j) over (2048, Q2).
#include <random>
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include "core.hpp"
#include "functions.hpp"
#include "RgswSample.hpp"
#include "LweSample.hpp"

class VIACrypto
{
public:
    static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t q)
    {
        return static_cast<uint64_t>((static_cast<unsigned __int128>(a) * b) % q);
    }
    static inline uint64_t toMod(int64_t x, uint64_t q)
    {
        int64_t r = x % static_cast<int64_t>(q);
        return static_cast<uint64_t>(r < 0 ? r + static_cast<int64_t>(q) : r);
    }

    VIACrypto()
        : n1q11(DEGREE1, MODULUS_Q1_1), n1q12(DEGREE1, MODULUS_Q1_2), n1q2(DEGREE1, MODULUS_Q2),
          n2q2(DEGREE2, MODULUS_Q2), n2q3(DEGREE2, MODULUS_Q3), rng(0x56494143ULL),
          s1(DEGREE1), s2(DEGREE2), s1_q11_ntt(DEGREE1), s1_q12_ntt(DEGREE1), s1_q2_ntt(DEGREE1),
          s2_q2_ntt(DEGREE2), s2_q3_ntt(DEGREE2), s2_q3(DEGREE2), s2emb_q2_ntt(DEGREE1),
          t1(DEGREE1), t2(DEGREE1), t3(DEGREE1)
    {
        // s1: ternary over R_{2048} (KGen_1 range [-1, 1]).
        // Distribution knobs (experiments only; defaults = values below):
        //   VIA_S1_HALF=h      s1 uniform on [-h, h]          (default 1, ternary)
        //   VIA_S2_HALF=h      s2 uniform on [-h, h]          (default 1, ternary)
        //   VIA_S2_SIGMA=x     s2 rounded Gaussian, sigma x   (overrides VIA_S2_HALF)
        //   VIA_E1_SIGMA=x     ring-1 (n1, Q1) encryption error sigma (default STDDEV1)
        //   VIA_E2_SIGMA=x     ring-2 error sigma: RLWE/RGSW over (n2, Q2) and the
        //                      ring-switching key (encrypted under s2)  (default STDDEV2)
        e1_sigma = envDouble("VIA_E1_SIGMA", STDDEV1);
        e2_sigma = envDouble("VIA_E2_SIGMA", STDDEV2);
        const int64_t h1 = static_cast<int64_t>(envDouble("VIA_S1_HALF", static_cast<double>(HALF_MODULUS_KGEN1)));
        std::uniform_int_distribution<int64_t> k1(-h1, h1);
        for (auto& v : s1) v = k1(rng);
        // s2 over R_{512}. Upstream KGen_2 nominally draws from
        // [-HALF_MODULUS_KGEN2, HALF_MODULUS_KGEN2] = [-379, 379] (and mapped
        // negatives with "+2", so the key was garbage). With the code's CMux
        // gadget (key part: L = 4 digits of 2^5, i.e. 15 truncated bits of
        // q2 ~ 2^35) a secret of that size gives ~2^27 noise per CMux level,
        // above Delta2/2 = 2^26, so decryption fails. We therefore use a ternary
        // s2 (like s1). The secret distribution only affects client key
        // generation; server work and message sizes are identical.
        // VIA_S2_HALF=<h> overrides the range for experiments.
        const int64_t h2 = static_cast<int64_t>(envDouble("VIA_S2_HALF", static_cast<double>(HALF_MODULUS_KGEN1)));
        const double s2_sigma = envDouble("VIA_S2_SIGMA", 0.0);
        if (s2_sigma > 0) {
            gaussian(s2.data(), DEGREE2, s2_sigma);
        } else {
            std::uniform_int_distribution<int64_t> k2(-h2, h2);
            for (auto& v : s2) v = k2(rng);
        }

        for (uint64_t i = 0; i < DEGREE1; i++)
        {
            s1_q11_ntt[i] = toMod(s1[i], MODULUS_Q1_1);
            s1_q12_ntt[i] = toMod(s1[i], MODULUS_Q1_2);
            s1_q2_ntt[i] = toMod(s1[i], MODULUS_Q2);
        }
        n1q11.ComputeForward(s1_q11_ntt.data(), s1_q11_ntt.data(), 1, 1);
        n1q12.ComputeForward(s1_q12_ntt.data(), s1_q12_ntt.data(), 1, 1);
        n1q2.ComputeForward(s1_q2_ntt.data(), s1_q2_ntt.data(), 1, 1);
        std::fill(s2emb_q2_ntt.begin(), s2emb_q2_ntt.end(), 0);
        for (uint64_t i = 0; i < DEGREE2; i++)
        {
            s2_q2_ntt[i] = toMod(s2[i], MODULUS_Q2);
            s2_q3[i] = toMod(s2[i], MODULUS_Q3);
            s2_q3_ntt[i] = s2_q3[i];
            s2emb_q2_ntt[i * (DEGREE1 / DEGREE2)] = s2_q2_ntt[i];
        }
        n2q2.ComputeForward(s2_q2_ntt.data(), s2_q2_ntt.data(), 1, 1);
        n2q3.ComputeForward(s2_q3_ntt.data(), s2_q3_ntt.data(), 1, 1);
        n1q2.ComputeForward(s2emb_q2_ntt.data(), s2emb_q2_ntt.data(), 1, 1);
    }

    // ---- sampling -------------------------------------------------------
    void uniform(uint64_t* out, uint64_t n, uint64_t q)
    {
        std::uniform_int_distribution<uint64_t> d(0, q - 1);
        for (uint64_t i = 0; i < n; i++) out[i] = d(rng);
    }
    void gaussian(int64_t* out, uint64_t n, double stddev)
    {
        std::normal_distribution<double> d(0.0, stddev);
        for (uint64_t i = 0; i < n; i++) out[i] = static_cast<int64_t>(std::llround(d(rng)));
    }

    // ---- RLWE over (2048, Q1 = Q1_1*Q1_2), two NTT-form limbs ----------
    // m: coefficient form, values in [0, Q1).
    void encRlweQ1(const uint64_t* m, RlweSampleQ1& ct)
    {
        std::vector<int64_t> e(DEGREE1);
        gaussian(e.data(), DEGREE1, e1_sigma);
        encLimb(m, e.data(), MODULUS_Q1_1, n1q11, s1_q11_ntt.data(), ct.rlweMask1, ct.rlweVal1);
        encLimb(m, e.data(), MODULUS_Q1_2, n1q12, s1_q12_ntt.data(), ct.rlweMask2, ct.rlweVal2);
    }
    // Debug: returns m + e in [0, Q1).
    void decRlweQ1(const RlweSampleQ1& ct, uint64_t* out)
    {
        decLimb(ct.rlweMask1, ct.rlweVal1, MODULUS_Q1_1, n1q11, s1_q11_ntt.data(), t1.data());
        decLimb(ct.rlweMask2, ct.rlweVal2, MODULUS_Q1_2, n1q12, s1_q12_ntt.data(), t2.data());
        for (uint64_t i = 0; i < DEGREE1; i++)
            out[i] = (mulmod(t1[i], CRT_COEFFICIENT1, MODULUS_Q1) + mulmod(t2[i], CRT_COEFFICIENT2, MODULUS_Q1)) % MODULUS_Q1;
    }

    // RGSW over Q1 of the constant polynomial mu (0/1).
    void encRgswQ1Const(uint64_t mu, uint64_t L, uint64_t logB, RgswSampleQ1& C)
    {
        std::vector<uint64_t> m(DEGREE1, 0);
        for (uint64_t i = 0; i < L; i++)
        {
            uint64_t g = gadgetQ1(L, logB, i);
            std::fill(m.begin(), m.end(), 0);
            m[0] = mulmod(mu, g, MODULUS_Q1);
            encRlweQ1(m.data(), C.msgRlev.rlweSamples[i]);
            for (uint64_t k = 0; k < DEGREE1; k++)
                m[k] = mulmod(toMod(-static_cast<int64_t>(mu) * s1[k], MODULUS_Q1), g, MODULUS_Q1);
            encRlweQ1(m.data(), C.keyRlev.rlweSamples[i]);
        }
    }

    // ---- RLWE over (512, Q2), NTT form ----------------------------------
    void encRlweQ2(const uint64_t* m, RlweSampleQ2& ct)
    {
        std::vector<int64_t> e(DEGREE2);
        gaussian(e.data(), DEGREE2, e2_sigma);
        for (uint64_t i = 0; i < DEGREE2; i++) t1[i] = toMod(e[i], MODULUS_Q2);
        n2q2.ComputeForward(t1.data(), t1.data(), 1, 1);
        for (uint64_t i = 0; i < DEGREE2; i++) t2[i] = m[i] % MODULUS_Q2;
        n2q2.ComputeForward(t2.data(), t2.data(), 1, 1);
        uniform(ct.rlweMask, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseMultMod(ct.rlweVal, ct.rlweMask, s2_q2_ntt.data(), DEGREE2, MODULUS_Q2, 1);
        intel::hexl::EltwiseAddMod(ct.rlweVal, ct.rlweVal, t1.data(), DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseAddMod(ct.rlweVal, ct.rlweVal, t2.data(), DEGREE2, MODULUS_Q2);
    }
    // Debug: m + e in [0, Q2), coefficient form.
    void decRlweQ2(const RlweSampleQ2& ct, uint64_t* out)
    {
        intel::hexl::EltwiseMultMod(t1.data(), ct.rlweMask, s2_q2_ntt.data(), DEGREE2, MODULUS_Q2, 1);
        intel::hexl::EltwiseSubMod(t1.data(), ct.rlweVal, t1.data(), DEGREE2, MODULUS_Q2);
        n2q2.ComputeInverse(out, t1.data(), 1, 1);
    }
    // Debug: (2048, Q2) NTT-form RLWE under s1.
    void decRlweQ2Deg1(const uint64_t* a, const uint64_t* b, uint64_t* out)
    {
        intel::hexl::EltwiseMultMod(t1.data(), a, s1_q2_ntt.data(), DEGREE1, MODULUS_Q2, 1);
        intel::hexl::EltwiseSubMod(t1.data(), b, t1.data(), DEGREE1, MODULUS_Q2);
        n1q2.ComputeInverse(out, t1.data(), 1, 1);
    }

    void encRgswQ2Const(uint64_t mu, uint64_t L_key, uint64_t logB_key, uint64_t L_msg, uint64_t logB_msg, RgswSampleQ2& C)
    {
        std::vector<uint64_t> m(DEGREE2, 0);
        for (uint64_t i = 0; i < L_msg; i++)
        {
            std::fill(m.begin(), m.end(), 0);
            m[0] = mulmod(mu, gadgetPow2(CEIL_LOG_Q2, L_msg, logB_msg, i, MODULUS_Q2), MODULUS_Q2);
            encRlweQ2(m.data(), C.msgRlev.rlweSamples[i]);
        }
        for (uint64_t i = 0; i < L_key; i++)
        {
            uint64_t g = gadgetPow2(CEIL_LOG_Q2, L_key, logB_key, i, MODULUS_Q2);
            for (uint64_t k = 0; k < DEGREE2; k++)
                m[k] = mulmod(toMod(-static_cast<int64_t>(mu) * s2[k], MODULUS_Q2), g, MODULUS_Q2);
            encRlweQ2(m.data(), C.keyRlev.rlweSamples[i]);
        }
    }

    // Ring-switching key from s1 (2048) to s2(X^4), over Q2.
    void genRsk(uint64_t L, uint64_t logB, RingSwitchingKey& rsk)
    {
        std::vector<int64_t> e(DEGREE1);
        for (uint64_t j = 0; j < L; j++)
        {
            uint64_t g = gadgetPow2(CEIL_LOG_Q2, L, logB, j, MODULUS_Q2);
            gaussian(e.data(), DEGREE1, e2_sigma);
            for (uint64_t i = 0; i < DEGREE1; i++)
            {
                t1[i] = toMod(e[i], MODULUS_Q2);
                t2[i] = mulmod(toMod(s1[i], MODULUS_Q2), g, MODULUS_Q2);
            }
            n1q2.ComputeForward(t1.data(), t1.data(), 1, 1);
            n1q2.ComputeForward(t2.data(), t2.data(), 1, 1);
            uniform(rsk.rlweMasks[j], DEGREE1, MODULUS_Q2);
            intel::hexl::EltwiseMultMod(rsk.rlweVals[j], rsk.rlweMasks[j], s2emb_q2_ntt.data(), DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseAddMod(rsk.rlweVals[j], rsk.rlweVals[j], t1.data(), DEGREE1, MODULUS_Q2);
            intel::hexl::EltwiseAddMod(rsk.rlweVals[j], rsk.rlweVals[j], t2.data(), DEGREE1, MODULUS_Q2);
        }
    }

    // ---- decoding of the final responses --------------------------------
    static int64_t decodeQ4(uint64_t m)
    {
        int64_t c = (m >= static_cast<uint64_t>(HALF_MODULUS_Q4)) ? static_cast<int64_t>(m) - static_cast<int64_t>(MODULUS_Q4) : static_cast<int64_t>(m);
        int64_t v = static_cast<int64_t>(std::llround(static_cast<double>(c) / static_cast<double>(DELTA_FINAL)));
        return ((v % static_cast<int64_t>(MODULUS_P)) + MODULUS_P) % MODULUS_P;
    }
    static uint64_t switchQ3toQ4(uint64_t x)
    {
        uint64_t r = static_cast<uint64_t>(static_cast<double>(x) * (static_cast<double>(MODULUS_Q4) / static_cast<double>(MODULUS_Q3)) + 0.5);
        return r >= MODULUS_Q4 ? r - MODULUS_Q4 : r;
    }
    // RLWE response: mask NTT-form mod Q3 (degree 512), val coefficient-form mod Q4.
    void decryptResponse(const RlweSampleQ2& ans, int64_t* out)
    {
        intel::hexl::EltwiseMultMod(t1.data(), ans.rlweMask, s2_q3_ntt.data(), DEGREE2, MODULUS_Q3, 1);
        n2q3.ComputeInverse(t1.data(), t1.data(), 1, 1);
        for (uint64_t j = 0; j < DEGREE2; j++)
        {
            uint64_t as = switchQ3toQ4(t1[j]);
            out[j] = decodeQ4((ans.rlweVal[j] + MODULUS_Q4 - as) % MODULUS_Q4);
        }
    }
    // LWE response: mask mod Q3 (n = 512), val mod Q4.
    int64_t decryptResponseLwe(const LweSampleQ2& ans)
    {
        unsigned __int128 acc = 0;
        for (uint64_t j = 0; j < DEGREE2; j++) acc += static_cast<unsigned __int128>(ans.lweMask[j] % MODULUS_Q3) * s2_q3[j];
        uint64_t as = switchQ3toQ4(static_cast<uint64_t>(acc % MODULUS_Q3));
        return decodeQ4((ans.lweVal % MODULUS_Q4 + MODULUS_Q4 - as) % MODULUS_Q4);
    }

    static uint64_t gadgetPow2(uint64_t logq, uint64_t L, uint64_t logB, uint64_t i, uint64_t q)
    {
        return (UINT64_C(1) << (logq - L * logB + i * logB)) % q;
    }
    static uint64_t gadgetQ1(uint64_t L, uint64_t logB, uint64_t i) { return gadgetPow2(CEIL_LOG_Q1, L, logB, i, MODULUS_Q1); }

    static double envDouble(const char* name, double dflt)
    {
        const char* e = std::getenv(name);
        return (e && *e) ? std::atof(e) : dflt;
    }

    double e1_sigma = STDDEV1, e2_sigma = STDDEV2;
    intel::hexl::NTT n1q11, n1q12, n1q2, n2q2, n2q3;
    std::mt19937_64 rng;
    std::vector<int64_t> s1, s2;
    std::vector<uint64_t> s1_q11_ntt, s1_q12_ntt, s1_q2_ntt, s2_q2_ntt, s2_q3_ntt, s2_q3, s2emb_q2_ntt;

private:
    std::vector<uint64_t> t1, t2, t3;

    void encLimb(const uint64_t* m, const int64_t* e, uint64_t q, intel::hexl::NTT& ntt, const uint64_t* s_ntt, uint64_t* a, uint64_t* b)
    {
        for (uint64_t i = 0; i < DEGREE1; i++) { t1[i] = m[i] % q; t2[i] = toMod(e[i], q); }
        ntt.ComputeForward(t1.data(), t1.data(), 1, 1);
        ntt.ComputeForward(t2.data(), t2.data(), 1, 1);
        uniform(a, DEGREE1, q);
        intel::hexl::EltwiseMultMod(b, a, s_ntt, DEGREE1, q, 1);
        intel::hexl::EltwiseAddMod(b, b, t2.data(), DEGREE1, q);
        intel::hexl::EltwiseAddMod(b, b, t1.data(), DEGREE1, q);
    }
    void decLimb(const uint64_t* a, const uint64_t* b, uint64_t q, intel::hexl::NTT& ntt, const uint64_t* s_ntt, uint64_t* out)
    {
        intel::hexl::EltwiseMultMod(t3.data(), a, s_ntt, DEGREE1, q, 1);
        intel::hexl::EltwiseSubMod(t3.data(), b, t3.data(), DEGREE1, q);
        ntt.ComputeInverse(out, t3.data(), 1, 1);
    }
};

#pragma once
#include "core.hpp"
#include "FHE_LWE.hpp"
#include "MlevKSKSample.hpp"

class VIA_CClient
{
public:
    VIA_CClient() : nttQ1_1C(DEGREE1, MODULUS_Q1_1C), nttQ1_2C(DEGREE1, MODULUS_Q1_2C), nttQ3C(DEGREE2, MODULUS_Q3C)
    {
        temp_poly = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        for (uint64_t i = 0; i < 10; i++)
        {
            sk[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        }
    }
	
    ~VIA_CClient()
    {
        std::free(temp_poly);
        for (uint64_t i = 0; i < 10; i++)
        {
            std::free(sk[i]);
        }
    }

    // Setup 
    void Setup(std::vector<RlevSampleQ1> ksk1, std::vector<RlevSampleQ1> ksk2, RingSwitchingKey& rsk, RlevSampleQ1& toRGSW)
    {
        std::cout << "\u2502   \u25BA" << " Client Setup Time: ";
        auto start = std::chrono::high_resolution_clock::now();
        fhe.GenRsk(CEIL_LOG_Q1C, GADGET_RSK_LC, LOG_GADGET_RSK_BASEC, rsk);

        ring_embed(1024, 2048, 0, fhe.secretKey_LWE1, temp_poly);
        fhe.encryptRlevQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, temp_poly, ksk1[0]);
        ring_embed(1024, 2048, 0, fhe.secretKey_LWE1 + 1024, temp_poly);
        fhe.encryptRlevQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, temp_poly, ksk2[0]);
        for (int i = 0; i < 10; i++)
        {
            KGen_1C(sk[i], DEGREE1);
            ring_embed(1024, 2048, 0, sk[i], temp_poly);
            fhe.encryptRlevQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, temp_poly, ksk1[i+1]);
            ring_embed(1024, 2048, 0, sk[i] + 1024, temp_poly);
            fhe.encryptRlevQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, temp_poly, ksk2[i+1]);
        }
        intel::hexl::EltwiseMultMod(fhe.temp_message1_1, fhe.secretKey1_1, fhe.secretKey1_1, DEGREE1, MODULUS_Q1_1C, 1);
        nttQ1_1C.ComputeInverse(fhe.temp_message1_1, fhe.temp_message1_1, 1, 1);
        intel::hexl::EltwiseMultMod(fhe.temp_message1_2, fhe.secretKey1_2, fhe.secretKey1_2, DEGREE1, MODULUS_Q1_2C, 1);
        nttQ1_2C.ComputeInverse(fhe.temp_message1_2, fhe.temp_message1_2, 1, 1);
        CRT(temp_poly, fhe.temp_message1_1, fhe.temp_message1_2, DEGREE1);
        fhe.encryptRlevQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, temp_poly, toRGSW);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << static_cast<double>(duration.count())/1000.0 << " ms" << std::endl;
        std::cout << "\u2502" << std::endl;
    }

    // query
    void Query(const uint64_t index, LweSampleQ1* qu)
    {
        std::cout << "\u2502   \u25BA" << " Query Time: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t idx_i = (index % (ROW * COL)) / COL;
        uint64_t idx_j = index % COL;
        uint64_t idx_k = index / (ROW * COL);
        std::vector<uint64_t> ControlBits(LOG_ROW);
        for (int i = 0; i < LOG_ROW; i++)
        {
            ControlBits[LOG_ROW - i - 1] = (idx_i >> i) & 1;
        }
        std::vector<uint64_t> SelectBits(LOG_COL);
        for (int i = 0; i < LOG_COL; i++)
        {
            SelectBits[i] = (idx_j >> i) & 1;
        }
        uint64_t msg;
        for (int i = 0; i < LOG_ROW; i++)
        {
            for (int j = 0; j < GADGET1_LC; j++)
            {
                msg = ControlBits[i] << (CEIL_LOG_Q1C - GADGET1_LC * LOG_GADGET1_BASEC + j);
                fhe.encryptLweQ1(msg, qu[i * GADGET1_LC + j]);
            }
        }
        for (int i = 0; i < LOG_COL; i++)
        {
            for (int j = 0; j < GADGET2_LC; j++)
            {
                fhe.encryptLweQ1(SelectBits[i] << (CEIL_LOG_Q1C - GADGET2_LC * LOG_GADGET2_BASEC + j), qu[LOG_ROW * GADGET1_LC + i * GADGET2_LC + j]);
            }
        }
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < GADGET2_LC; j++)
            {
                fhe.encryptLweQ1((idx_k & 1) << (CEIL_LOG_Q1C - GADGET2_LC * LOG_GADGET2_BASEC + j), qu[LOG_ROW * GADGET1_LC + LOG_COL * GADGET2_LC + i * GADGET2_LC + j]);
                idx_k >>= 1;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << static_cast<double>(duration.count())/1000.0 << " ms" << std::endl;
    }

    // query for blinded extraction
    void Query_Extraction(const uint64_t index[2], LweSampleQ1* qu)
    {
        std::cout << "\u2502   \u25BA" << " Query Time with Blinded Extraction: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t idx_i = (index[0] % (ROW * COL)) / COL;
        uint64_t idx_j = index[0] % COL;
        std::vector<uint64_t> ControlBits(LOG_ROW);
        for (int i = 0; i < LOG_ROW; i++)
        {
            ControlBits[LOG_ROW - i - 1] = (idx_i >> i) & 1;
        }
        std::vector<uint8_t> SelectBits(LOG_COL);
        for (int i = 0; i < LOG_COL; i++)
        {
            SelectBits[i] = (idx_j >> i) & 1;
        }
        for (int i = 0; i < LOG_ROW; i++)
        {
            for (int j = 0; j < GADGET1_LC; j++)
            {
                fhe.encryptLweQ1(ControlBits[i] << (CEIL_LOG_Q1C - GADGET1_LC * LOG_GADGET1_BASEC + j), qu[i * GADGET1_LC + j]);
            }
        }
        for (int i = 0; i < LOG_COL; i++)
        {
            for (int j = 0; j < GADGET2_LC; j++)
            {
                fhe.encryptLweQ1(SelectBits[i] << (CEIL_LOG_Q1C - GADGET2_LC * LOG_GADGET2_BASEC + j), qu[LOG_ROW * GADGET1_LC + i * GADGET2_LC + j]);
            }
        }
        uint64_t idx_k = index[1];
        for (int i = 0; i < 11; i++)
        {
            for (int j = 0; j < GADGET2_LC; j++)
            {
                fhe.encryptLweQ1((idx_k & 1) << (CEIL_LOG_Q1C - GADGET2_LC * LOG_GADGET2_BASEC + j), qu[LOG_ROW * GADGET1_LC + LOG_COL * GADGET2_LC + i * GADGET2_LC + j]);
                idx_k >>= 1;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << static_cast<double>(duration.count())/1000.0 << " ms" << std::endl;
    }

    // Recover RLWE
	void Recover(RlweSampleQ2 ans, int64_t* result)
    {
        std::cout << "\u2502   \u25BA" << " Recover Time: ";
        auto start = std::chrono::high_resolution_clock::now();
        intel::hexl::EltwiseMultMod(fhe.temp_prod1, ans.rlweMask, fhe.secretKey2_ntt, DEGREE2, MODULUS_Q3C, 1);
        nttQ3C.ComputeInverse(fhe.temp_prod1, fhe.temp_prod1, 1, 1);
        modSwitch(MODULUS_Q3C, MODULUS_Q4C, fhe.temp_prod1, fhe.temp_prod1, DEGREE2);
        intel::hexl::EltwiseSubMod(fhe.temp_message1, ans.rlweVal, fhe.temp_prod1, DEGREE2, MODULUS_Q4C);
        for (int j = 0; j < DEGREE2; j++)
        {
            if (fhe.temp_message1[j] < HALF_MODULUS_Q4C)
            {
                result[j] = static_cast<int64_t>(round(static_cast<double>(fhe.temp_message1[j]/static_cast<double>(DELTAC_FINAL))));
            }
            else
            {
                result[j] = static_cast<int64_t>(round(static_cast<double>(fhe.temp_message1[j] - MODULUS_Q4C)/static_cast<double>(DELTAC_FINAL)));
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << duration.count() << " \u03BCs" << std::endl;
        std::cout << "\u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;

    }

	// Recover LWE
	void Recover_Extraction(LweSampleQ2 ans, int64_t result)
    {
        std::cout << "\u2502   \u25BA" << " Recover Time with Blinded Extraction: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t temp;
        intel::hexl::EltwiseMultMod(fhe.temp_prod1, ans.lweMask, fhe.secretKey2_ntt, DEGREE2, MODULUS_Q3C, 1);
        temp = sum(fhe.temp_prod1, DEGREE2) % MODULUS_Q3C < ans.lweVal;
        temp = (temp > ans.lweVal) ? (ans.lweVal + MODULUS_Q3C - temp) : (ans.lweVal - temp);
        modSwitch(MODULUS_Q3C, MODULUS_Q4C, &temp, &temp, 1);
        result = (temp < HALF_MODULUS_Q4C) ? (static_cast<int64_t>(round(static_cast<double>(temp)/static_cast<double>(DELTAC_FINAL)))) : (static_cast<int64_t>(round(static_cast<double>(temp - MODULUS_Q4C)/static_cast<double>(DELTAC_FINAL))));
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << duration.count() << " \u03BCs" << std::endl;
        std::cout << "\u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;
    }

	FHE_LWE fhe;
    uint64_t* temp_poly;
    uint64_t* sk[10];
    intel::hexl::NTT nttQ1_1C;
    intel::hexl::NTT nttQ1_2C;
    intel::hexl::NTT nttQ3C;
    LweSampleQ1 temp_lweSample;
};
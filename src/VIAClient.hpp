#pragma once
#include "FHE_LWE.hpp"

class VIAClient
{
public:
	// constructor
	VIAClient() : nttQ3(DEGREE2, MODULUS_Q3)
    {
        temp_poly = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    }
	// destructor
	~VIAClient()
    {
        std::free(temp_poly);
    }

	// Query
	void Query(uint64_t& index, RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel)
    {
        std::cout << "\u2502   \u25BA" << " Query Time: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t idx_i = (index % (ROW * COL)) / COL;
        uint64_t idx_j = index % COL;
        uint64_t idx_k = index / (ROW * COL);
        std::vector<uint8_t> ControlBits(LOG_ROW);
        for (int i = 0; i < LOG_ROW; i++)
        {
            ControlBits[LOG_ROW - i - 1] = (idx_i >> i) & 1;
        }
        std::vector<uint8_t> SelectBits(LOG_COL - 3);
        for (int i = 0; i < LOG_COL - 3; i++)
        {
            SelectBits[i] = (idx_j >> i) & 1;
        }

        std::fill(fhe.temp_message1, fhe.temp_message1 + DEGREE1, 0);
        for (int i = 0; i < 8; i++)
        {
            if (i != (ControlBits[0] + 2 * ControlBits[1] + 4 * ControlBits[2]))
            {
                fhe.encryptRlweQ1(fhe.temp_message1, c_rot[i]);
            }
        }
        fhe.temp_message1[idx_k] = DELTA1;
        fhe.encryptRlweQ1(fhe.temp_message1, c_rot[ControlBits[0] + 2 * ControlBits[1] + 4 * ControlBits[2]]);

        std::fill(fhe.temp_message1, fhe.temp_message1 + DEGREE1, 0);
        for (int i = 3; i < LOG_ROW; i++)
        {
            fhe.temp_message1[0] = static_cast<uint64_t>(ControlBits[i]);
            fhe.encryptRgswQ1(CEIL_LOG_Q1, GADGET1_L, GADGET1_L, LOG_GADGET1_BASE, LOG_GADGET1_BASE, fhe.temp_message1, C_ctrl[i - 3]);
        }
        for (int i = 0; i < LOG_COL - 3; i++)
        {
            fhe.temp_message1[0] = static_cast<uint64_t>(SelectBits[i]);
            fhe.encryptRgswQ2(CEIL_LOG_Q2, GADGET2_L_1, GADGET2_L_2, LOG_GADGET2_BASE_1, LOG_GADGET2_BASE_2, fhe.temp_message1, C_sel[i]);
        }
        fhe.GenRsk(CEIL_LOG_Q2, GADGET_RSK_L, LOG_GADGET_RSK_BASE, rsk);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << static_cast<double>(duration.count())/1000.0 << " ms" << std::endl;
    }

	// Query for blinded extraction
	void Query_Extraction(uint64_t index[2], RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel)
    {
        std::cout << "\u2502   \u25BA" << " Query Time with Blinded Extraction: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t idx_i = (index[0] % (ROW * COL)) / COL;
        uint64_t idx_j = index[0] % COL;
        std::vector<uint8_t> ControlBits(LOG_ROW);
        for (int i = 0; i < LOG_ROW; i++)
        {
            ControlBits[LOG_ROW - i - 1] = (idx_i >> i) & 1;
        }
        std::vector<uint8_t> SelectBits(LOG_COL - 3);
        for (int i = 0; i < LOG_COL; i++)
        {
            SelectBits[i] = (idx_j >> i) & 1;
        }

        std::fill(fhe.temp_message1, fhe.temp_message1 + DEGREE1, 0);
        for (int i = 0; i < 8; i++)
        {
            if (i != (ControlBits[0] + 2 * ControlBits[1] + 4 * ControlBits[2]))
            {
                fhe.encryptRlweQ1(fhe.temp_message1, c_rot[i]);
            }
        }
        fhe.temp_message1[index[1]] = DELTA1;
        fhe.encryptRlweQ1(fhe.temp_message1, c_rot[ControlBits[0] + 2 * ControlBits[1] + 4 * ControlBits[2]]);

        std::fill(fhe.temp_message1, fhe.temp_message1 + DEGREE1, 0);
        for (int i = 3; i < LOG_ROW; i++)
        {
            fhe.temp_message1[0] = static_cast<uint64_t>(ControlBits[i]);
            fhe.encryptRgswQ1(CEIL_LOG_Q1, GADGET1_L, GADGET1_L, LOG_GADGET1_BASE, LOG_GADGET1_BASE, fhe.temp_message1, C_ctrl[i - 3]);
        }
        for (int i = 0; i < LOG_COL - 3; i++)
        {
            fhe.temp_message1[0] = static_cast<uint64_t>(SelectBits[i]);
            fhe.encryptRgswQ2(CEIL_LOG_Q2, GADGET2_L_1, GADGET2_L_2, LOG_GADGET2_BASE_1, LOG_GADGET2_BASE_2, fhe.temp_message1, C_sel[i]);
        }
        fhe.GenRsk(CEIL_LOG_Q2, GADGET_RSK_L, LOG_GADGET_RSK_BASE, rsk);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << static_cast<double>(duration.count())/1000.0 << " ms" << std::endl;
    }

	// Recover RLWE
	void Recover(RlweSampleQ2 ans[8], int64_t* result[8])
    {
        std::cout << "\u2502   \u25BA" << " Recover Time: ";
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 8; i++)
        {
            intel::hexl::EltwiseMultMod(fhe.temp_prod1, ans[i].rlweMask, fhe.secretKey2_ntt, DEGREE2, MODULUS_Q3, 1);
            nttQ3.ComputeInverse(fhe.temp_prod1, fhe.temp_prod1, 1, 1);
            modSwitch(MODULUS_Q3, MODULUS_Q4, fhe.temp_prod1, fhe.temp_prod1, DEGREE2);
            intel::hexl::EltwiseSubMod(fhe.temp_message1, ans[i].rlweVal, fhe.temp_prod1, DEGREE2, MODULUS_Q4);
            for (int j = 0; j < DEGREE2; j++)
            {
                if (fhe.temp_message1[j] < HALF_MODULUS_Q4)
                {
                    result[i][j] = static_cast<int64_t>(round(static_cast<double>(fhe.temp_message1[j]/static_cast<double>(DELTA_FINAL))));
                }
                else
                {
                    result[i][j] = static_cast<int64_t>(round(static_cast<double>(fhe.temp_message1[j] - MODULUS_Q4)/static_cast<double>(DELTA_FINAL)));
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << duration.count() << " \u03BCs" << std::endl;
        std::cout << "\u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;

    }

	// Recover LWE
	void Recover_Extraction(LweSampleQ2 ans[8], int64_t result[8])
    {
        std::cout << "\u2502   \u25BA" << " Recover Time with Blinded Extraction: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t temp;
        for (int i = 0; i < 8; i++)
        {
            intel::hexl::EltwiseMultMod(fhe.temp_prod1, ans[i].lweMask, fhe.secretKey2_ntt, DEGREE2, MODULUS_Q3, 1);
            temp = sum(fhe.temp_prod1, DEGREE2) % MODULUS_Q3 < ans[i].lweVal;
            temp = (temp > ans[i].lweVal) ? (ans[i].lweVal + MODULUS_Q3 - temp) : (ans[i].lweVal - temp);
            modSwitch(MODULUS_Q3, MODULUS_Q4, &temp, &temp, 1);
            result[i] = (temp < HALF_MODULUS_Q4) ? (static_cast<int64_t>(round(static_cast<double>(temp)/static_cast<double>(DELTA_FINAL)))) : (static_cast<int64_t>(round(static_cast<double>(temp - MODULUS_Q4)/static_cast<double>(DELTA_FINAL))));
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << duration.count() << " \u03BCs" << std::endl;
        std::cout << "\u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;
    }


	FHE_LWE fhe;
    intel::hexl::NTT nttQ3;
    uint64_t* temp_poly;
};


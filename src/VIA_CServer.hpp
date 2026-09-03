#pragma once
#include "core.hpp"
#include "functions.hpp"
#include "RgswSample.hpp"
#include "LweSample.hpp"

class VIA_CServer
{
public:
// constructor
	VIA_CServer() : nttQ1_1C(DEGREE1, MODULUS_Q1_1C), nttQ1_2C(DEGREE1, MODULUS_Q1_2C), nttQ2CDeg1(DEGREE1, MODULUS_Q2C), nttQ3C(DEGREE2, MODULUS_Q3C), nttQ3CDeg1(DEGREE1, MODULUS_Q3C)
    {
        FirstDimVec.resize(ROW);
        SecondDimVec.resize(COL);
        temp_message1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        temp_message1_1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        temp_message1_2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        temp_prod1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        gadget_poly1.resize(GADGET_LC_MAX);
        gadget_poly2.resize(GADGET_LC_MAX);
        for (int i = 0; i < GADGET_LC_MAX; i++)
        {
            gadget_poly1[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
            gadget_poly2[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        }

        uint64_t dbSize = (ROW * COL > 1 << 21) ? (1 << 21) : (ROW * COL);
        database.resize(dbSize);
        for (uint64_t i = 0; i < dbSize; i++)
        {
            database[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
            UniformModQ2(database[i], DEGREE1);
        }
        for (uint64_t i = 0; i < LOG_ROW; i++)
        {
            C_ctrl.emplace_back(RgswSampleQ1(GADGET1_LC, GADGET1_LC));
        }
        for (uint64_t i = 0; i < LOG_COL; i++)
        {
            C_selQ1.emplace_back(RgswSampleQ1(GADGET2_LC, GADGET2_LC));
        }
        for (uint64_t i = 0; i < LOG_COL; i++)
        {
            C_sel.emplace_back(RgswSampleQ2Deg1(GADGET2_LC, GADGET2_LC));
        }
        for (uint64_t i = 0; i < 11; i++)
        {
            C_rotQ1.emplace_back(RgswSampleQ1(GADGET2_LC, GADGET2_LC));
        }
        for (uint64_t i = 0; i < 11; i++)
        {
            C_rot.emplace_back(RgswSampleQ2Deg1(GADGET2_LC, GADGET2_LC));
        }
    }

	// destructor
	~VIA_CServer()
    {
        std::free(temp_message1);
        std::free(temp_message1_1);
        std::free(temp_message1_2);
        std::free(temp_prod1);
        for (int i = 0; i < GADGET_LC_MAX; i++)
        {
            std::free(gadget_poly1[i]);
            std::free(gadget_poly2[i]);
        }
        uint64_t dbSize = (ROW * COL > 1 << 21) ? (1 << 21) : (ROW * COL);
        for (uint64_t i = 0; i < dbSize; i++) 
        {
            std::free(database[i]);
        }
    }

	// init_Database
	void SetupDatabase()
    {
        std::cout << "\u2502   \u25BA" << " Setup Database: ";
        uint64_t dbSize = ROW * COL;
        uint64_t loop = 1;
        if (LOG_ROW + LOG_COL > 21)
        {
            dbSize = (1 << 21);
            loop = (ROW * COL) / dbSize;
        }
        auto start = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < loop; j++)
        {
            for (uint64_t i = 0; i < dbSize; i++)
            {
                nttQ2CDeg1.ComputeForward(database[i], database[i], 1, 1);
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << static_cast<double>(duration) / 1000000.0 << " s" << std::endl;
    }


    void Answer(std::vector<RlevSampleQ1> ksk1, std::vector<RlevSampleQ1> ksk2, RingSwitchingKey& rsk, RlevSampleQ1& toRGSW, LweSampleQ1* qu, RlweSampleQ2& result)
    {
        std::cout << "\u2502   \u25BA" << " Answer Time: " << std::endl;
        // decompress the query
        std::cout << "\u2502     \u25C7" << " Decompress the query: ";
        auto start_Decompress = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < LOG_ROW*GADGET1_LC + GADGET2_LC*(LOG_COL + 2); i++)
        {
            CRT_C(temp_message1, qu[i].lweMask1, qu[i].lweMask2, DEGREE1);
            halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[0], temp_message1, temp_rlweSampleQ1_1);
            halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[0], temp_message1, temp_rlweSampleQ1_2);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
            for (int j = 0; j < 9; j++)
            {
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[j], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[j], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
            }
            if (i < LOG_ROW*GADGET1_LC)
            {
                uint64_t m = i / GADGET1_LC;
                uint64_t n = i % GADGET1_LC;
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[10], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[10], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_1.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, toRGSW, temp_message1, C_ctrl[m].keyRlev.rlweSamples[n]);
                intel::hexl::EltwiseAddMod(C_ctrl[m].keyRlev.rlweSamples[n].rlweMask1, C_ctrl[m].keyRlev.rlweSamples[n].rlweMask1, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].keyRlev.rlweSamples[n].rlweMask2, C_ctrl[m].keyRlev.rlweSamples[n].rlweMask2, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask2, DEGREE1, MODULUS_Q1_2C);
            }
            else if (i < LOG_ROW*GADGET1_LC + GADGET2_LC*LOG_COL)
            {
                uint64_t m = (i - LOG_ROW*GADGET1_LC) / GADGET2_LC;
                uint64_t n = (i - LOG_ROW*GADGET1_LC) % GADGET2_LC;
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[10], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[10], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_1.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, toRGSW, temp_message1, C_selQ1[m].keyRlev.rlweSamples[n]);
                intel::hexl::EltwiseAddMod(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, 1, 1);
                CRT_C(temp_message1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_sel[m].keyRlev.rlweSamples[n].rlweVal, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_sel[m].keyRlev.rlweSamples[n].rlweVal, C_sel[m].keyRlev.rlweSamples[n].rlweVal, 1, 1);
                nttQ1_1C.ComputeForward(C_selQ1[m].keyRlev.rlweSamples[n].rlweVal1, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal1, 1, 1);
                nttQ1_2C.ComputeForward(C_selQ1[m].keyRlev.rlweSamples[n].rlweVal2, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal2, 1, 1);
                CRT_C(temp_message1, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal1, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_sel[m].keyRlev.rlweSamples[n].rlweMask, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_sel[m].keyRlev.rlweSamples[n].rlweMask, C_sel[m].keyRlev.rlweSamples[n].rlweMask, 1, 1);
            }
            else
            {
                uint64_t m = (i - LOG_ROW*GADGET1_LC - GADGET2_LC*LOG_COL) / GADGET2_LC;
                uint64_t n = (i - LOG_ROW*GADGET1_LC - GADGET2_LC*LOG_COL) % GADGET2_LC;
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[10], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[10], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_1.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, toRGSW, temp_message1, C_rotQ1[m].keyRlev.rlweSamples[n]);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, 1, 1);
                CRT_C(temp_message1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_rot[m].keyRlev.rlweSamples[n].rlweMask, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_rot[m].keyRlev.rlweSamples[n].rlweMask, C_rot[m].keyRlev.rlweSamples[n].rlweMask, 1, 1);
                nttQ1_1C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, 1, 1);
                CRT_C(temp_message1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_rot[m].keyRlev.rlweSamples[n].rlweVal, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_rot[m].keyRlev.rlweSamples[n].rlweVal, C_rot[m].keyRlev.rlweSamples[n].rlweVal, 1, 1);
            }
        }
        auto end_Decompress = std::chrono::high_resolution_clock::now();
        auto duration_Decompress = std::chrono::duration_cast<std::chrono::microseconds>(end_Decompress - start_Decompress).count();
        std::cout << static_cast<double>(duration_Decompress) / 1000.0 << " ms" << std::endl;

        // DMux
        std::cout << "\u2502     \u25C7" << " DMux: ";
        auto start_DMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < LOG_ROW; i++)
        for (uint64_t i = 0; i < LOG_ROW; i++)
        {
            for (int64_t j = (1 << i) - 1; j >= 0; j--)
            {
                DMuxQ1(CEIL_LOG_Q1C, GADGET1_LC, GADGET1_LC, LOG_GADGET1_BASEC, LOG_GADGET1_BASEC, C_ctrl[i], FirstDimVec[j], FirstDimVec[j << 1], FirstDimVec[(j << 1) + 1]);
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
            nttQ1_1C.ComputeInverse(FirstDimVec[i].rlweMask1, FirstDimVec[i].rlweMask1, 1, 1);
            nttQ1_2C.ComputeInverse(FirstDimVec[i].rlweMask2, FirstDimVec[i].rlweMask2, 1, 1);
            CRT_C(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, FirstDimVec[i].rlweMask1, DEGREE1);
            nttQ2CDeg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_1C.ComputeInverse(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_2C.ComputeInverse(FirstDimVec[i].rlweVal2, FirstDimVec[i].rlweVal2, 1, 1);
            CRT_C(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, FirstDimVec[i].rlweVal1, DEGREE1);
            nttQ2CDeg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
        }
        auto end_ModSwitch = std::chrono::high_resolution_clock::now();
        auto duration_ModSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_ModSwitch - start_ModSwitch).count();
        std::cout << static_cast<double>(duration_ModSwitch) / 1000.0 << " ms" << std::endl;

        // First Dimension
        std::cout << "\u2502     \u25C7" << " First Dimension: ";
        auto start_FirstDim = std::chrono::high_resolution_clock::now();
        for (uint64_t j = 0; j < COL; j++)
        {
            intel::hexl::EltwiseMultMod(SecondDimVec[0].rlweMask, FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2C, 1);
            intel::hexl::EltwiseMultMod(SecondDimVec[0].rlweVal, FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2C, 1);
            for (uint64_t i = 1; i < ROW; i++)
            {
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2C, 1);
                intel::hexl::EltwiseAddMod(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, temp_prod1, DEGREE1, MODULUS_Q2C);
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2C, 1);
                intel::hexl::EltwiseAddMod(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, temp_prod1, DEGREE1, MODULUS_Q2C);
            }
        }
        auto end_FirstDim = std::chrono::high_resolution_clock::now();
        auto duration_FirstDim = std::chrono::duration_cast<std::chrono::microseconds>(end_FirstDim - start_FirstDim).count();
        std::cout << static_cast<double>(duration_FirstDim) / 1000.0 << " ms" << std::endl;

        // CMux
        std::cout << "\u2502     \u25C7" << " CMux: ";
        auto start_CMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = LOG_COL; i > 3; i--)
        {
            for (uint64_t j = 0; j < (1 << (i - 1)); j++)
            {
                CMuxQ2Deg1(CEIL_LOG_Q2C, GADGET2_LC, GADGET2_LC, LOG_GADGET2_BASEC, LOG_GADGET2_BASEC, C_sel[LOG_COL - i], SecondDimVec[j << 1], SecondDimVec[(j << 1) + 1], SecondDimVec[j]);
            }
        }
        auto end_CMux = std::chrono::high_resolution_clock::now();
        auto duration_CMux = std::chrono::duration_cast<std::chrono::microseconds>(end_CMux - start_CMux).count();
        std::cout << static_cast<double>(duration_CMux) / 1000.0 << " ms" << std::endl;

        // CRot
        std::cout << "\u2502     \u25C7" << " CRot: ";
        auto start_CRot = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < DEGREE1 - 1; i++)
        {
            temp_rlweSampleQ2Deg1.rlweMask[0] = SecondDimVec[0].rlweMask[i+1];
            temp_rlweSampleQ2Deg1.rlweVal[0] = SecondDimVec[0].rlweVal[i+1];
        }
        temp_rlweSampleQ2Deg1.rlweMask[DEGREE1 - 1] = -SecondDimVec[0].rlweMask[0];
        temp_rlweSampleQ2Deg1.rlweVal[DEGREE1 - 1] = -SecondDimVec[0].rlweVal[0];
        CMuxQ2Deg1(CEIL_LOG_Q2, GADGET2_LC, GADGET2_LC, LOG_GADGET2_BASEC, LOG_GADGET2_BASEC, C_rot[0], SecondDimVec[0], temp_rlweSampleQ2Deg1, SecondDimVec[0]);
        for (uint64_t i = 0; i < DEGREE1 - 2; i++)
        {
            temp_rlweSampleQ2Deg1.rlweMask[0] = SecondDimVec[0].rlweMask[i+1];
            temp_rlweSampleQ2Deg1.rlweVal[0] = SecondDimVec[0].rlweVal[i+1];
        }
        temp_rlweSampleQ2Deg1.rlweMask[DEGREE1 - 1] = -SecondDimVec[0].rlweMask[1];
        temp_rlweSampleQ2Deg1.rlweVal[DEGREE1 - 2] = -SecondDimVec[0].rlweVal[0];
        temp_rlweSampleQ2Deg1.rlweVal[DEGREE1 - 1] = -SecondDimVec[0].rlweVal[1];
        temp_rlweSampleQ2Deg1.rlweVal[DEGREE1 - 2] = -SecondDimVec[0].rlweVal[0];
        CMuxQ2Deg1(CEIL_LOG_Q2, GADGET2_LC, GADGET2_LC, LOG_GADGET2_BASEC, LOG_GADGET2_BASEC, C_rot[1], SecondDimVec[0], temp_rlweSampleQ2Deg1, SecondDimVec[0]);
        auto end_CRot = std::chrono::high_resolution_clock::now();
        auto duration_CRot = std::chrono::duration_cast<std::chrono::microseconds>(end_CRot - start_CRot).count();
        std::cout << static_cast<double>(duration_CRot) / 1000.0 << " ms" << std::endl;


        // RingSwitch
        std::cout << "\u2502     \u25C7" << " RingSwitch: ";
        auto start_RingSwitch = std::chrono::high_resolution_clock::now();

        nttQ2CDeg1.ComputeInverse(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, 1, 1);
        modSwitch(MODULUS_Q2C, MODULUS_Q3C, SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, DEGREE1);
        nttQ2CDeg1.ComputeInverse(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, 1, 1); 
        modSwitch(MODULUS_Q2C, MODULUS_Q3C, SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, DEGREE1);
        GadgetDecompositionDeg1(SecondDimVec[0].rlweMask, gadget_poly1, CEIL_LOG_Q3C, GADGET_RSK_LC, LOG_GADGET_RSK_BASEC);
        for (int i = 0; i < GADGET_RSK_LC; i++)
        {
            nttQ3CDeg1.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
        }
        intel::hexl::EltwiseMultMod(SecondDimVec[0].rlweMask, gadget_poly1[0], rsk.rlweMasks[0], DEGREE1, MODULUS_Q3C, 1);
        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[0], rsk.rlweVals[0], DEGREE1, MODULUS_Q3C, 1);
        intel::hexl::EltwiseSubMod(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, temp_prod1, DEGREE1, MODULUS_Q3C);
        for (int i = 1; i < GADGET_RSK_LC; i++)
        {
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rsk.rlweMasks[i], DEGREE1, MODULUS_Q3C, 1);
            intel::hexl::EltwiseAddMod(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, temp_prod1, DEGREE1, MODULUS_Q3C);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rsk.rlweVals[i], DEGREE1, MODULUS_Q3C, 1);
            intel::hexl::EltwiseSubMod(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, temp_prod1, DEGREE1, MODULUS_Q3C);
        }
        std::fill(temp_message1, temp_message1 + DEGREE1, 0);
        intel::hexl::EltwiseSubMod(SecondDimVec[0].rlweMask, temp_message1, SecondDimVec[0].rlweMask, DEGREE1, MODULUS_Q3C);
        nttQ3CDeg1.ComputeInverse(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, 1, 1);
        ring_project(0, SecondDimVec[0].rlweMask, result.rlweMask);
        nttQ3CDeg1.ComputeInverse(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, 1, 1);
        ring_project(0, SecondDimVec[0].rlweVal, result.rlweVal);
        modSwitch(MODULUS_Q3C, MODULUS_Q4C, result.rlweVal, result.rlweVal, DEGREE2);
        auto end_RingSwitch = std::chrono::high_resolution_clock::now();
        auto duration_RingSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_RingSwitch - start_RingSwitch).count();
        std::cout << static_cast<double>(duration_RingSwitch) / 1000.0 << " ms" << std::endl;

        std::cout << "\u2502   \u25BA" << " Total Answer Time: "  <<   static_cast<double>(duration_Decompress + duration_DMux + duration_FirstDim + duration_CMux + duration_CRot + duration_RingSwitch) / 1000.0 << " ms" << std::endl;
        std::cout << "\u2502 " << std::endl;
    }

    void Answer_Extraction(std::vector<RlevSampleQ1> ksk1, std::vector<RlevSampleQ1> ksk2, RingSwitchingKey& rsk, RlevSampleQ1& toRGSW, LweSampleQ1* qu, LweSampleQ2& result)
    {
        std::cout << "\u2502   \u25BA" << " Answer Time with Blinded Extraction: " << std::endl;
        // decompress the query
        std::cout << "\u2502     \u25C7" << " Decompress the query: ";
        auto start_Decompress = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < LOG_ROW*GADGET1_LC + GADGET2_LC*(LOG_COL + 11); i++)
        {
            CRT_C(temp_message1, qu[i].lweMask1, qu[i].lweMask2, DEGREE1);
            halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[0], temp_message1, temp_rlweSampleQ1_1);
            halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[0], temp_message1, temp_rlweSampleQ1_2);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
            intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
            for (int j = 0; j < 9; j++)
            {
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[j], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[j], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1.rlweVal1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1.rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
            }
            if (i < LOG_ROW*GADGET1_LC)
            {
                uint64_t m = i / GADGET1_LC;
                uint64_t n = i % GADGET1_LC;
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[10], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[10], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_1.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, C_ctrl[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, toRGSW, temp_message1, C_ctrl[m].keyRlev.rlweSamples[n]);
                intel::hexl::EltwiseAddMod(C_ctrl[m].keyRlev.rlweSamples[n].rlweMask1, C_ctrl[m].keyRlev.rlweSamples[n].rlweMask1, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_ctrl[m].keyRlev.rlweSamples[n].rlweMask2, C_ctrl[m].keyRlev.rlweSamples[n].rlweMask2, C_ctrl[m].msgRlev.rlweSamples[n].rlweMask2, DEGREE1, MODULUS_Q1_2C);
            }
            else if (i < LOG_ROW*GADGET1_LC + GADGET2_LC*LOG_COL)
            {
                uint64_t m = (i - LOG_ROW*GADGET1_LC) / GADGET2_LC;
                uint64_t n = (i - LOG_ROW*GADGET1_LC) % GADGET2_LC;
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[10], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[10], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_1.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_selQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, toRGSW, temp_message1, C_selQ1[m].keyRlev.rlweSamples[n]);
                intel::hexl::EltwiseAddMod(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_selQ1[m].msgRlev.rlweSamples[n].rlweMask2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, 1, 1);
                CRT_C(temp_message1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_selQ1[m].keyRlev.rlweSamples[n].rlweMask2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_sel[m].keyRlev.rlweSamples[n].rlweVal, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_sel[m].keyRlev.rlweSamples[n].rlweVal, C_sel[m].keyRlev.rlweSamples[n].rlweVal, 1, 1);
                nttQ1_1C.ComputeForward(C_selQ1[m].keyRlev.rlweSamples[n].rlweVal1, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal1, 1, 1);
                nttQ1_2C.ComputeForward(C_selQ1[m].keyRlev.rlweSamples[n].rlweVal2, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal2, 1, 1);
                CRT_C(temp_message1, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal1, C_selQ1[m].keyRlev.rlweSamples[n].rlweVal2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_sel[m].keyRlev.rlweSamples[n].rlweMask, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_sel[m].keyRlev.rlweSamples[n].rlweMask, C_sel[m].keyRlev.rlweSamples[n].rlweMask, 1, 1);
            }
            else
            {
                uint64_t m = (i - LOG_ROW*GADGET1_LC - GADGET2_LC*LOG_COL) / GADGET2_LC;
                uint64_t n = (i - LOG_ROW*GADGET1_LC - GADGET2_LC*LOG_COL) % GADGET2_LC;
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk1[10], temp_message1, temp_rlweSampleQ1_1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, ksk2[10], temp_message1, temp_rlweSampleQ1_2);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweMask1, temp_rlweSampleQ1_2.rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, temp_rlweSampleQ1_2.rlweVal1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask2, temp_rlweSampleQ1_1.rlweMask2, temp_rlweSampleQ1_2.rlweMask2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_1.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, C_rotQ1[m].msgRlev.rlweSamples[n].rlweVal2, temp_rlweSampleQ1_2.rlweVal2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(temp_rlweSampleQ1.rlweMask2, temp_rlweSampleQ1.rlweMask2, 1, 1);
                CRT_C(temp_message1, temp_rlweSampleQ1.rlweMask1, temp_rlweSampleQ1.rlweMask2, DEGREE1);
                halfExternalProdQ1(CEIL_LOG_Q1C, GADGET_L_CONV, LOG_GADGET_BASE_CONV, toRGSW, temp_message1, C_rotQ1[m].keyRlev.rlweSamples[n]);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask1, DEGREE1, MODULUS_Q1_1C);
                intel::hexl::EltwiseAddMod(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].msgRlev.rlweSamples[n].rlweMask2, DEGREE1, MODULUS_Q1_2C);
                nttQ1_1C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, 1, 1);
                CRT_C(temp_message1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_rot[m].keyRlev.rlweSamples[n].rlweMask, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_rot[m].keyRlev.rlweSamples[n].rlweMask, C_rot[m].keyRlev.rlweSamples[n].rlweMask, 1, 1);
                nttQ1_1C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, 1, 1);
                nttQ1_2C.ComputeInverse(C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, 1, 1);
                CRT_C(temp_message1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask1, C_rotQ1[m].keyRlev.rlweSamples[n].rlweMask2, DEGREE1);
                modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, C_rot[m].keyRlev.rlweSamples[n].rlweVal, DEGREE1);
                nttQ2CDeg1.ComputeForward(C_rot[m].keyRlev.rlweSamples[n].rlweVal, C_rot[m].keyRlev.rlweSamples[n].rlweVal, 1, 1);
            }
        }
        auto end_Decompress = std::chrono::high_resolution_clock::now();
        auto duration_Decompress = std::chrono::duration_cast<std::chrono::microseconds>(end_Decompress - start_Decompress).count();
        std::cout << static_cast<double>(duration_Decompress) / 1000.0 << " ms" << std::endl;

        // DMux
        std::cout << "\u2502     \u25C7" << " DMux: ";
        auto start_DMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < LOG_ROW; i++)
        for (uint64_t i = 0; i < LOG_ROW; i++)
        {
            for (int64_t j = (1 << i) - 1; j >= 0; j--)
            {
                DMuxQ1(CEIL_LOG_Q1C, GADGET1_LC, GADGET1_LC, LOG_GADGET1_BASEC, LOG_GADGET1_BASEC, C_ctrl[i], FirstDimVec[j], FirstDimVec[j << 1], FirstDimVec[(j << 1) + 1]);
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
            nttQ1_1C.ComputeInverse(FirstDimVec[i].rlweMask1, FirstDimVec[i].rlweMask1, 1, 1);
            nttQ1_2C.ComputeInverse(FirstDimVec[i].rlweMask2, FirstDimVec[i].rlweMask2, 1, 1);
            CRT_C(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, FirstDimVec[i].rlweMask1, DEGREE1);
            nttQ2CDeg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_1C.ComputeInverse(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
            nttQ1_2C.ComputeInverse(FirstDimVec[i].rlweVal2, FirstDimVec[i].rlweVal2, 1, 1);
            CRT_C(temp_message1, FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal2, DEGREE1);
            modSwitch(MODULUS_Q1C, MODULUS_Q2C, temp_message1, FirstDimVec[i].rlweVal1, DEGREE1);
            nttQ2CDeg1.ComputeForward(FirstDimVec[i].rlweVal1, FirstDimVec[i].rlweVal1, 1, 1);
        }
        auto end_ModSwitch = std::chrono::high_resolution_clock::now();
        auto duration_ModSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_ModSwitch - start_ModSwitch).count();
        std::cout << static_cast<double>(duration_ModSwitch) / 1000.0 << " ms" << std::endl;

        // First Dimension
        std::cout << "\u2502     \u25C7" << " First Dimension: ";
        auto start_FirstDim = std::chrono::high_resolution_clock::now();
        for (uint64_t j = 0; j < COL; j++)
        {
            intel::hexl::EltwiseMultMod(SecondDimVec[0].rlweMask, FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2C, 1);
            intel::hexl::EltwiseMultMod(SecondDimVec[0].rlweVal, FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2C, 1);
            for (uint64_t i = 1; i < ROW; i++)
            {
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2C, 1);
                intel::hexl::EltwiseAddMod(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, temp_prod1, DEGREE1, MODULUS_Q2C);
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2C, 1);
                intel::hexl::EltwiseAddMod(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, temp_prod1, DEGREE1, MODULUS_Q2C);
            }
        }
        auto end_FirstDim = std::chrono::high_resolution_clock::now();
        auto duration_FirstDim = std::chrono::duration_cast<std::chrono::microseconds>(end_FirstDim - start_FirstDim).count();
        std::cout << static_cast<double>(duration_FirstDim) / 1000.0 << " ms" << std::endl;

        // CMux
        std::cout << "\u2502     \u25C7" << " CMux: ";
        auto start_CMux = std::chrono::high_resolution_clock::now();
        for (uint64_t i = LOG_COL; i > 3; i--)
        {
            for (uint64_t j = 0; j < (1 << (i - 1)); j++)
            {
                CMuxQ2Deg1(CEIL_LOG_Q2C, GADGET2_LC, GADGET2_LC, LOG_GADGET2_BASEC, LOG_GADGET2_BASEC, C_sel[LOG_COL - i], SecondDimVec[j << 1], SecondDimVec[(j << 1) + 1], SecondDimVec[j]);
            }
        }
        auto end_CMux = std::chrono::high_resolution_clock::now();
        auto duration_CMux = std::chrono::duration_cast<std::chrono::microseconds>(end_CMux - start_CMux).count();
        std::cout << static_cast<double>(duration_CMux) / 1000.0 << " ms" << std::endl;

        // CRot
        std::cout << "\u2502     \u25C7" << " CRot: ";
        auto start_CRot = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < 11; j++)
        {
            for (uint64_t i = 0; i < DEGREE1 - (1 << j); i++)
            {
                temp_rlweSampleQ2Deg1.rlweMask[i] = SecondDimVec[0].rlweMask[i + (1 << j)];
                temp_rlweSampleQ2Deg1.rlweVal[i] = SecondDimVec[0].rlweVal[i + (1 << j)];
            }
            for (uint64_t i = DEGREE1 - (1 << j); i < DEGREE1; i++)
            {
                temp_rlweSampleQ2Deg1.rlweMask[i] = SecondDimVec[1].rlweMask[i - (DEGREE1 - (1 << j))];
                temp_rlweSampleQ2Deg1.rlweVal[i] = SecondDimVec[1].rlweVal[i - (DEGREE1 - (1 << j))];
            }
            CMuxQ2Deg1(CEIL_LOG_Q2, GADGET2_LC, GADGET2_LC, LOG_GADGET2_BASEC, LOG_GADGET2_BASEC, C_rot[j], SecondDimVec[0], temp_rlweSampleQ2Deg1, SecondDimVec[0]);
        }
        auto end_CRot = std::chrono::high_resolution_clock::now();
        auto duration_CRot = std::chrono::duration_cast<std::chrono::microseconds>(end_CRot - start_CRot).count();
        std::cout << static_cast<double>(duration_CRot) / 1000.0 << " ms" << std::endl;


        // RingSwitch
        std::cout << "\u2502     \u25C7" << " RingSwitch: ";
        auto start_RingSwitch = std::chrono::high_resolution_clock::now();

        nttQ2CDeg1.ComputeInverse(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, 1, 1);
        modSwitch(MODULUS_Q2C, MODULUS_Q3C, SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, DEGREE1);
        nttQ2CDeg1.ComputeInverse(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, 1, 1); 
        modSwitch(MODULUS_Q2C, MODULUS_Q3C, SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, DEGREE1);
        GadgetDecompositionDeg1(SecondDimVec[0].rlweMask, gadget_poly1, CEIL_LOG_Q3C, GADGET_RSK_LC, LOG_GADGET_RSK_BASEC);
        for (int i = 0; i < GADGET_RSK_LC; i++)
        {
            nttQ3CDeg1.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
        }
        intel::hexl::EltwiseMultMod(SecondDimVec[0].rlweMask, gadget_poly1[0], rsk.rlweMasks[0], DEGREE1, MODULUS_Q3C, 1);
        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[0], rsk.rlweVals[0], DEGREE1, MODULUS_Q3C, 1);
        intel::hexl::EltwiseSubMod(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, temp_prod1, DEGREE1, MODULUS_Q3C);
        for (int i = 1; i < GADGET_RSK_LC; i++)
        {
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rsk.rlweMasks[i], DEGREE1, MODULUS_Q3C, 1);
            intel::hexl::EltwiseAddMod(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, temp_prod1, DEGREE1, MODULUS_Q3C);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rsk.rlweVals[i], DEGREE1, MODULUS_Q3C, 1);
            intel::hexl::EltwiseSubMod(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, temp_prod1, DEGREE1, MODULUS_Q3C);
        }
        std::fill(temp_message1, temp_message1 + DEGREE1, 0);
        intel::hexl::EltwiseSubMod(SecondDimVec[0].rlweMask, temp_message1, SecondDimVec[0].rlweMask, DEGREE1, MODULUS_Q3C);
        nttQ3CDeg1.ComputeInverse(SecondDimVec[0].rlweMask, SecondDimVec[0].rlweMask, 1, 1);
        ring_project(0, SecondDimVec[0].rlweMask, result.lweMask);
        nttQ3CDeg1.ComputeInverse(SecondDimVec[0].rlweVal, SecondDimVec[0].rlweVal, 1, 1);
        modSwitch(MODULUS_Q3C, MODULUS_Q4C, SecondDimVec[0].rlweVal, &result.lweVal, 1);
        auto end_RingSwitch = std::chrono::high_resolution_clock::now();
        auto duration_RingSwitch = std::chrono::duration_cast<std::chrono::microseconds>(end_RingSwitch - start_RingSwitch).count();
        std::cout << static_cast<double>(duration_RingSwitch) / 1000.0 << " ms" << std::endl;

        std::cout << "\u2502   \u25BA" << " Total Answer Time with Blinded Extraction: "  <<   static_cast<double>(duration_Decompress + duration_DMux + duration_FirstDim + duration_CMux + duration_CRot + duration_RingSwitch) / 1000.0 << " ms" << std::endl;
        std::cout << "\u2502 " << std::endl;
    }

    // half extetnal product for Q1
    void halfExternalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ1& rlevSample, uint64_t* poly, RlweSampleQ1& result)
    {
        GadgetDecompositionDeg1(poly, gadget_poly1, log_modulus, gadget_l, log_gadget_base);
        for (int i = 0; i < gadget_l; i++)
        {
            nttQ1_2C.ComputeForward(gadget_poly2[i], gadget_poly1[i], 1, 1);
            nttQ1_1C.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
        }
        intel::hexl::EltwiseMultMod(result.rlweMask1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweMask1, DEGREE1, MODULUS_Q1_1C, 1);
        intel::hexl::EltwiseMultMod(result.rlweVal1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweVal1, DEGREE1, MODULUS_Q1_1C, 1);

        intel::hexl::EltwiseMultMod(result.rlweMask2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweMask2, DEGREE1, MODULUS_Q1_1C, 2);
        intel::hexl::EltwiseMultMod(result.rlweVal2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweVal2, DEGREE1, MODULUS_Q1_1C, 2);

        for (int i = 1; i < gadget_l; i++)
        {
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweMask1, DEGREE1, MODULUS_Q1_1C, 1);
            intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_prod1, DEGREE1, MODULUS_Q1_1C);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweVal1, DEGREE1, MODULUS_Q1_1C, 1);
            intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_prod1, DEGREE1, MODULUS_Q1_1C);

            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweMask2, DEGREE1, MODULUS_Q1_1C, 2);
            intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_prod1, DEGREE1, MODULUS_Q1_2C);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweVal2, DEGREE1, MODULUS_Q1_1C, 2);
            intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_prod1, DEGREE1, MODULUS_Q1_2C);
        }
    }

    // external product for Q1
    void externalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result)
    {
        nttQ1_1C.ComputeInverse(temp_message1_1, rlweSample.rlweMask1, 1, 1);
        nttQ1_2C.ComputeInverse(temp_message1_2, rlweSample.rlweMask2, 1, 1);
        CRT_C(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
        halfExternalProdQ1(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1, result);
        nttQ1_1C.ComputeInverse(temp_message1_1, rlweSample.rlweVal1, 1, 1);
        nttQ1_2C.ComputeInverse(temp_message1_2, rlweSample.rlweVal2, 1, 1);
        CRT_C(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
        halfExternalProdQ1(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1, temp_rlweSampleQ1);
        intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_rlweSampleQ1.rlweMask1, DEGREE1, MODULUS_Q1_1C);
        intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_rlweSampleQ1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
        intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_rlweSampleQ1.rlweMask2, DEGREE1, MODULUS_Q1_2C);
        intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_rlweSampleQ1.rlweVal2, DEGREE1, MODULUS_Q1_2C);
    }

    // DMux for Q1
    void DMuxQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result0, RlweSampleQ1& result1)
    {
        externalProdQ1(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, rlweSample, temp_rlweSampleQ1);
        intel::hexl::EltwiseSubMod(result0.rlweMask1, rlweSample.rlweMask1, temp_rlweSampleQ1.rlweMask1, DEGREE1, MODULUS_Q1_1C);
        intel::hexl::EltwiseSubMod(result0.rlweVal1, rlweSample.rlweVal1, temp_rlweSampleQ1.rlweVal1, DEGREE1, MODULUS_Q1_1C);
        intel::hexl::EltwiseSubMod(result0.rlweMask2, rlweSample.rlweMask2, temp_rlweSampleQ1.rlweMask2, DEGREE1, MODULUS_Q1_2C);
        intel::hexl::EltwiseSubMod(result0.rlweVal2, rlweSample.rlweVal2, temp_rlweSampleQ1.rlweVal2, DEGREE1, MODULUS_Q1_2C);

        intel::hexl::EltwiseSubMod(result1.rlweMask1, result0.rlweMask1, rlweSample.rlweMask1, DEGREE1, MODULUS_Q1_1C);
        intel::hexl::EltwiseSubMod(result1.rlweVal1, result0.rlweVal1, rlweSample.rlweVal1, DEGREE1, MODULUS_Q1_1C);
        intel::hexl::EltwiseSubMod(result1.rlweMask2, result0.rlweMask2, rlweSample.rlweMask2, DEGREE1, MODULUS_Q1_2C);    
        intel::hexl::EltwiseSubMod(result1.rlweVal2, result0.rlweVal2, rlweSample.rlweVal2, DEGREE1, MODULUS_Q1_2C);
    }

    // half product for Q2 degree1
    void halfExternalProdQ2Deg1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ2Deg1& rlevSample, uint64_t* poly, RlweSampleQ2Deg1& result)
    {
        GadgetDecompositionDeg1(poly, gadget_poly1, log_modulus, gadget_l, log_gadget_base);
        for (int i = 0; i < gadget_l; i++)
        {
            nttQ2CDeg1.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
        }

        intel::hexl::EltwiseMultMod(result.rlweMask, gadget_poly1[0], rlevSample.rlweSamples[0].rlweMask, DEGREE1, MODULUS_Q2C, 1);
        intel::hexl::EltwiseMultMod(result.rlweVal, gadget_poly1[0], rlevSample.rlweSamples[0].rlweVal, DEGREE1, MODULUS_Q2C, 1);

        for (int i = 1; i < gadget_l; i++)
        {
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweMask, DEGREE1, MODULUS_Q2C, 1);
            intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, temp_prod1, DEGREE1, MODULUS_Q2C);
            intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweVal, DEGREE1, MODULUS_Q2C, 1);
            intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, temp_prod1, DEGREE1, MODULUS_Q2C);
        }
    }

    // external product for Q2
    void externalProdQ2Deg1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2Deg1& rgswSample, RlweSampleQ2Deg1& rlweSample, RlweSampleQ2Deg1& result)
    {
        nttQ2CDeg1.ComputeInverse(temp_message1_1, rlweSample.rlweMask, 1, 1);
        halfExternalProdQ2Deg1(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1_1, result);
        nttQ2CDeg1.ComputeInverse(temp_message1_1, rlweSample.rlweVal, 1, 1);
        halfExternalProdQ2Deg1(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1_1, temp_rlweSampleQ2Deg1);
        intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, temp_rlweSampleQ2Deg1.rlweMask, DEGREE1, MODULUS_Q2C);
        intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, temp_rlweSampleQ2Deg1.rlweVal, DEGREE1, MODULUS_Q2C);
    }

    // CMux for Q2
    void CMuxQ2Deg1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2Deg1& rgswSample, RlweSampleQ2Deg1& rlweSample0, RlweSampleQ2Deg1& rlweSample1, RlweSampleQ2Deg1& result)
    {
        intel::hexl::EltwiseSubMod(temp_rlweSampleQ2Deg1.rlweMask, rlweSample1.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2C);
        intel::hexl::EltwiseSubMod(temp_rlweSampleQ2Deg1.rlweVal, rlweSample1.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2C);
        externalProdQ2Deg1(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, temp_rlweSampleQ2Deg1, result);
        intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2C);
        intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2C);
    }

	std::vector<uint64_t*> database;
    std::vector<RlweSampleQ1> FirstDimVec;
    std::vector<RlweSampleQ2Deg1> SecondDimVec;
    intel::hexl::NTT nttQ1_1C;
    intel::hexl::NTT nttQ1_2C;
    intel::hexl::NTT nttQ2CDeg1;
    intel::hexl::NTT nttQ3C;
    intel::hexl::NTT nttQ3CDeg1;
    intel::hexl::NTT nttQ4C;
    uint64_t* temp_prod1;
    uint64_t* temp_message1;
    uint64_t* temp_message1_1;
    uint64_t* temp_message1_2;
    std::vector<uint64_t*> gadget_poly1;
    std::vector<uint64_t*> gadget_poly2;
    RlweSampleQ1 temp_rlweSampleQ1;
    RlweSampleQ1 temp_rlweSampleQ1_1;
    RlweSampleQ1 temp_rlweSampleQ1_2;
    RlweSampleQ2Deg1 temp_rlweSampleQ2Deg1;
    std::vector<RgswSampleQ1> C_ctrl;
    std::vector<RgswSampleQ1> C_selQ1;
    std::vector<RgswSampleQ2Deg1> C_sel;
    std::vector<RgswSampleQ1> C_rotQ1;
    std::vector<RgswSampleQ2Deg1> C_rot;
};
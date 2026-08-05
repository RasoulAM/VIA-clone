#pragma once
#include <fstream>
#include <iostream>
#include <chrono>
#include "core.hpp"
#include "functions.hpp"
#include "RgswSample.hpp"
#include "LweSample.hpp"

class VIAServer
{
public:
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

        uint64_t dbSize = (ROW * COL > 1 << 21) ? (1 << 21) : (ROW * COL);
        database.resize(dbSize);
        for (uint64_t i = 0; i < dbSize; i++)
        {
            database[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
            UniformModQ2(database[i], DEGREE1);
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
        uint64_t dbSize = (ROW * COL > 1 << 21) ? (1 << 21) : (ROW * COL);
        for (uint64_t i = 0; i < dbSize; i++) 
        {
            if (database[i] != nullptr) {
                std::free(database[i]);
                database[i] = nullptr;
            }
        }
    }

	// init_Database
	void SetupDatabase()
    {
        std::cout << "│   ►" << " Setup Database: ";
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
                nttQ2Deg1.ComputeForward(database[i], database[i], 1, 1);
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << static_cast<double>(duration) / 1000000.0 << " s" << std::endl;
        std::cout << "│ " <<std::endl;

        // Write database to disk
        std::ofstream outfile("/tmp/via_db.bin", std::ios::binary);
        if (outfile.is_open()) {
            uint64_t element_size = DEGREE1 * sizeof(uint64_t);
            for (uint64_t i = 0; i < dbSize; i++) {
                outfile.write(reinterpret_cast<const char*>(database[i]), element_size);
            }
            outfile.close();
            std::cout << "Server Long term Storage: " << (dbSize * element_size) << " bytes" << std::endl;
        } else {
            std::cerr << "Failed to open /tmp/via_db.bin for writing!" << std::endl;
        }
        
        // Free the memory to emulate disk-backed storage
        for (uint64_t i = 0; i < dbSize; i++) {
            std::free(database[i]);
            database[i] = nullptr;
        }
    }

    void LoadDatabase()
    {
        uint64_t dbSize = ROW * COL;
        if (LOG_ROW + LOG_COL > 21)
        {
            dbSize = (1 << 21);
        }
        auto start = std::chrono::high_resolution_clock::now();
        
        std::ifstream infile("/tmp/via_db.bin", std::ios::binary);
        if (infile.is_open()) {
            uint64_t element_size = DEGREE1 * sizeof(uint64_t);
            for (uint64_t i = 0; i < dbSize; i++) {
                if (database[i] == nullptr) {
                    database[i] = static_cast<uint64_t*>(std::aligned_alloc(64, element_size));
                }
                infile.read(reinterpret_cast<char*>(database[i]), element_size);
            }
            infile.close();
        } else {
            std::cerr << "Failed to open /tmp/via_db.bin for reading!" << std::endl;
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Server Data Load Time: " << duration << " us" << std::endl;
    }

	// Answer
	void Answer(RlweSampleQ1 c_rot[4], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel, RlweSampleQ2 result[8])
    {
        auto latency_start = std::chrono::high_resolution_clock::now();
        LoadDatabase();
        std::cout << "│   ►" << " Answer Time: " << std::endl;
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
        for (uint64_t j = 0; j < COL; j++)
        {
            intel::hexl::EltwiseMultMod(FirstDimResMasks[0], FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseMultMod(FirstDimResVals[0], FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2, 1);
            for (uint64_t i = 1; i < ROW; i++)
            {
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseAddMod(FirstDimResMasks[0], FirstDimResMasks[0], temp_prod1, DEGREE1, MODULUS_Q2);
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseAddMod(FirstDimResVals[0], FirstDimResVals[0], temp_prod1, DEGREE1, MODULUS_Q2);
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

        std::cout << "│   ►" << " Total Answer Time: " << static_cast<double>(duration_DMux + duration_ModSwitch + duration_FirstDim + duration_RingSwitch + duration_CMux + duration_FinalModSwitch) / 1000.0 << " ms" << std::endl;
        
        auto latency_end = std::chrono::high_resolution_clock::now();
        auto latency_duration = std::chrono::duration_cast<std::chrono::microseconds>(latency_end - latency_start).count();
        std::cout << "Server Response Latency: " << latency_duration << " us" << std::endl;
        std::cout << "│ " << std::endl;
    }

	void Answer_Extraction(RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel, LweSampleQ2 result[8])
    {
        auto latency_start = std::chrono::high_resolution_clock::now();
        LoadDatabase();
        std::cout << "│   ►" << " Answer Time with blind extraction:: " << std::endl;
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
        for (uint64_t j = 0; j < COL; j++)
        {
            intel::hexl::EltwiseMultMod(FirstDimResMasks[0], FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2, 1);
            intel::hexl::EltwiseMultMod(FirstDimResVals[0], FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2, 1);
            for (uint64_t i = 1; i < ROW; i++)
            {
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweMask1, database[0], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseAddMod(FirstDimResMasks[0], FirstDimResMasks[0], temp_prod1, DEGREE1, MODULUS_Q2);
                intel::hexl::EltwiseMultMod(temp_prod1, FirstDimVec[0].rlweVal1, database[0], DEGREE1, MODULUS_Q2, 1);
                intel::hexl::EltwiseAddMod(FirstDimResVals[0], FirstDimResVals[0], temp_prod1, DEGREE1, MODULUS_Q2);
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
        std::cout << "Server Response Latency: " << latency_duration << " us" << std::endl;
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


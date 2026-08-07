#include <iostream>
#include <cstdint>
#include <vector>
#include <array>
#include <chrono>
#include <memory>
#include <string>
#include "core.hpp"
#include "VIAClient.hpp"
#include "VIAServer.hpp"
#include "VIA_CClient.hpp"
#include "VIA_CServer.hpp"

int main()
{
    uint64_t  IsVIA_C, IsBlindedExtraction, ExecutionCount, WriteToDisk;
    std::cout << "Test VIA or VIA-C (0 or 1): " << std::endl;
    std::cin >> IsVIA_C;
    std::cout << "Test with or without Blinded Extraction (0 or 1): " << std::endl;
    std::cin >> IsBlindedExtraction;
    std::cout << "Execution count: " << std::endl;
    std::cin >> ExecutionCount;
    std::cout << "Write to disk (0 or 1): " << std::endl;
    std::cin >> WriteToDisk;

    if (IsVIA_C == 1)
    {
        static VIA_CClient client;
        static VIA_CServer server;
        printInfoVIA_C(IsBlindedExtraction);
        server.SetupDatabase();
        if (IsBlindedExtraction == 1)
        {
            std::vector<RlevSampleQ1> ksk1, ksk2;
            RingSwitchingKey rsk(GADGET_RSK_LC);
            RlevSampleQ1 toRGSW(GADGET_L_CONV);
            uint64_t index[2] = {20471, 2};
            LweSampleQ1 qu[LOG_ROW*GADGET1_LC + GADGET2_LC*(LOG_COL + 11)];
            LweSampleQ2 ans;
            int64_t result;
            for (int i = 0; i < ExecutionCount; i++)
            {
                std::cout << "| \u2756" << " Round " << i + 1 << ": computation time" << std::endl;
                server.SetupDatabase();
                for (int i = 0; i < 11; i++)
                {
                    ksk1.emplace_back(RlevSampleQ1(GADGET_L_CONV));
                    ksk2.emplace_back(RlevSampleQ1(GADGET_L_CONV));
                }

                client.Setup(ksk1, ksk2, rsk, toRGSW);
                client.Query_Extraction(index, qu);

                server.Answer_Extraction(ksk1, ksk2, rsk, toRGSW, qu, ans);

                client.Recover_Extraction(ans, result);
            }
        }
        else
        {
            std::vector<RlevSampleQ1> ksk1, ksk2;
            for (int i = 0; i < 11; i++)
            {
                ksk1.emplace_back(RlevSampleQ1(GADGET_L_CONV));
                ksk2.emplace_back(RlevSampleQ1(GADGET_L_CONV));
            }
            RingSwitchingKey rsk(GADGET_RSK_LC);
            RlevSampleQ1 toRGSW(GADGET_L_CONV);
            LweSampleQ1 qu[LOG_ROW*GADGET1_LC + GADGET2_LC*(LOG_COL + 2)];
            RlweSampleQ2 ans;
            int64_t* result = static_cast<int64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));
            for (int i = 0; i < ExecutionCount; i++)
            {
                client.Setup(ksk1, ksk2, rsk, toRGSW);

                client.Query(20471, qu);

                server.Answer(ksk1, ksk2, rsk, toRGSW, qu, ans);

                client.Recover(ans, result);
            }
        }
    }
    else
    {
        static VIAClient client;
        static VIAServer server;
        server.setWriteToDisk(WriteToDisk == 1);
        printInfoVIA(IsBlindedExtraction);
        server.SetupDatabase();
        if (IsBlindedExtraction == 1)
        {
            uint64_t index[2] = {20471, 2};
            RlweSampleQ1 c_rot[8];
            RingSwitchingKey rsk(GADGET_RSK_L);
            std::vector<RgswSampleQ1> C_ctrl;
            for (uint64_t i = 0; i < LOG_ROW - 3; i++)
            {
                C_ctrl.emplace_back(RgswSampleQ1(GADGET1_L, GADGET1_L));
            }
            std::vector<RgswSampleQ2> C_sel;
            for (uint64_t i = 0; i < LOG_COL - 3; i++)
            {
                C_sel.emplace_back(RgswSampleQ2(GADGET2_L_1, GADGET2_L_2));
            }
            LweSampleQ2 ans[8];
            int64_t result[8];

            for (int k = 0; k<ExecutionCount; k++)
            { 
                std::cout << "| \u2756" << " Round " << k + 1 << ": computation time" << std::endl;

                client.Query_Extraction(index, c_rot, rsk, C_ctrl, C_sel);

                server.Answer_Extraction(c_rot, rsk, C_ctrl, C_sel, ans);

                client.Recover_Extraction(ans, result);

                std::vector<uint64_t> expected(2048);
                server.getExpectedPolynomial(index[0], expected.data());
                if (result[0] != expected[index[1]]) {
                    std::cout << "│   ► Verification: FAILED!" << std::endl;
                } else {
                    std::cout << "│   ► Verification: SUCCESS!" << std::endl;
                }
            }
        }
        else
        {
            uint64_t index = 20471;
            RlweSampleQ1 c_rot[8];
            RingSwitchingKey rsk(GADGET_RSK_L);
            std::vector<RgswSampleQ1> C_ctrl;
            RlweSampleQ2 ans[8];
            int64_t* result[8];
            for (uint64_t i = 0; i < 8; i++)
            {
                result[i] = static_cast<int64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));
            }
            for (uint64_t i = 0; i < LOG_ROW - 3; i++)
            {
                C_ctrl.emplace_back(RgswSampleQ1(GADGET1_L, GADGET1_L));
            }
            std::vector<RgswSampleQ2> C_sel;
            for (uint64_t i = 0; i < LOG_COL - 3; i++)
            {
                C_sel.emplace_back(RgswSampleQ2(GADGET2_L_1, GADGET2_L_2));
            }
            for (int k = 0; k<ExecutionCount; k++)
            { 
                std::cout << "| \u2756" << " Round " << k + 1 << ": computation time" << std::endl;

                client.Query(index, c_rot, rsk, C_ctrl, C_sel);

                server.Answer(c_rot, rsk, C_ctrl, C_sel, ans);

                client.Recover(ans, result);

                std::vector<uint64_t> expected(2048);
                server.getExpectedPolynomial(index, expected.data());
                
                bool success = true;
                for (int i = 0; i < 2048; i++) { // DEGREE1 is 2048
                    if (result[0][i] != expected[i]) {
                        success = false;
                        break;
                    }
                }
                if (success) {
                    std::cout << "│   ► Verification: SUCCESS!" << std::endl;
                } else {
                    std::cout << "│   ► Verification: FAILED!" << std::endl;
                }
            }
        }
    }

    return 0;
}

#include <iostream>
#include <cstdint>
#include <vector>
#include <array>
#include <chrono>
#include <memory>
#include <string>
#include <cstdlib>
#include "core.hpp"
#include "VIAClient.hpp"
#include "VIAServer.hpp"
#include "VIA_CClient.hpp"
#include "VIA_CServer.hpp"

uint64_t LOG_ROW = 6;
uint64_t LOG_COL = 9;
uint64_t ROW = 1 << 6;
uint64_t COL = 1 << 9;

int main()
{
    std::cout << "LOG_ROW: " << std::endl;
    std::cin >> LOG_ROW;
    std::cout << "LOG_COL: " << std::endl;
    std::cin >> LOG_COL;
    ROW = 1ULL << LOG_ROW;
    COL = 1ULL << LOG_COL;
    uint64_t  IsVIA_C, IsBlindedExtraction, ExecutionCount;
    std::cout << "Test VIA or VIA-C (0 or 1): " << std::endl;
    std::cin >> IsVIA_C;
    std::cout << "Test with or without Blinded Extraction (0 or 1): " << std::endl;
    std::cin >> IsBlindedExtraction;
    std::cout << "Execution count: " << std::endl;
    std::cin >> ExecutionCount;
    int write_to_file;
    std::cout << "Write to file (0 or 1): ";
    std::cin >> write_to_file;
    // Optional: SSD format (0 = legacy NTT-form uint64, 1 = compact plaintext
    // bytes + NTT on load [default], 2 = compact, chunked load with prefetch).
    int ssd_format = 1;
    if (!(std::cin >> ssd_format)) ssd_format = 1;

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
            std::vector<LweSampleQ1> qu(LOG_ROW*GADGET1_LC + GADGET2_LC*(LOG_COL + 11));
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
                client.Query_Extraction(index, qu.data());

                server.Answer_Extraction(ksk1, ksk2, rsk, toRGSW, qu.data(), ans);

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
            std::vector<LweSampleQ1> qu(LOG_ROW*GADGET1_LC + GADGET2_LC*(LOG_COL + 2));
            RlweSampleQ2 ans;
            int64_t* result = static_cast<int64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));
            for (int i = 0; i < ExecutionCount; i++)
            {
                client.Setup(ksk1, ksk2, rsk, toRGSW);

                client.Query(20471, qu.data());

                server.Answer(ksk1, ksk2, rsk, toRGSW, qu.data(), ans);

                client.Recover(ans, result);
            }
        }
    }
    else
    {
        static VIAClient client;
        static VIAServer server;
        server.setWriteToFile(write_to_file == 1);
        server.setSsdFormat(ssd_format);
        printInfoVIA(IsBlindedExtraction);

        // Target record and the 8 columns whose polynomials come back in the
        // response (the CMux tree resolves only the low LOG_COL-3 column bits).
        const uint64_t target = 20471;
        const uint64_t coeff_idx = 2;  // blinded extraction: coefficient index
        const uint64_t idx_i = (target % (ROW * COL)) / COL;
        const uint64_t idx_j = target % COL;
        const uint64_t low_mask = (UINT64_C(1) << (LOG_COL - 3)) - 1;
        const uint64_t t_target = idx_j >> (LOG_COL - 3);
        std::vector<uint64_t> resp_polys(8);
        for (uint64_t t = 0; t < 8; t++)
            resp_polys[t] = idx_i * COL + ((t << (LOG_COL - 3)) | (idx_j & low_mask));
        server.setRetainPolys(resp_polys);

        // Negative tests: VIA_CORRUPT=db (perturb the stored DB entry) or
        // VIA_CORRUPT=resp (perturb the server response before decryption).
        const char* corrupt_env = std::getenv("VIA_CORRUPT");
        const std::string corrupt = corrupt_env ? corrupt_env : "";

        server.SetupDatabase();
        if (corrupt == "db")
            server.corruptDatabaseForTest(resp_polys[t_target], IsBlindedExtraction == 1 ? coeff_idx : 0);

        std::vector<std::vector<uint8_t>> truth(8);
        for (uint64_t t = 0; t < 8; t++)
            if (!server.getRetainedPlaintext(resp_polys[t], truth[t])) std::cout << "│   ► missing ground truth for poly " << resp_polys[t] << std::endl;

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
        RlweSampleQ1 c_rot[8];

        if (IsBlindedExtraction == 1)
        {
            uint64_t index[2] = {target, coeff_idx};
            LweSampleQ2 ans[8];
            int64_t result[8];

            for (int k = 0; k<ExecutionCount; k++)
            {
                std::cout << "| \u2756" << " Round " << k + 1 << ": computation time" << std::endl;

                client.Query_Extraction(index, c_rot, rsk, C_ctrl, C_sel);

                server.Answer_Extraction(c_rot, rsk, C_ctrl, C_sel, ans);
                server.VerifyLoadedDatabase();
                if (corrupt == "resp")
                {
                    ans[t_target].lweVal = (ans[t_target].lweVal + DELTA_FINAL) % MODULUS_Q4;
                    std::cout << "│   ► TEST: corrupted response " << t_target << std::endl;
                }

                client.Recover_Extraction(ans, result);

                uint64_t bad = 0;
                for (uint64_t t = 0; t < 8; t++)
                    if (truth[t].empty() || result[t] != static_cast<int64_t>(truth[t][coeff_idx])) bad++;
                std::cout << "│   ► Verification: checked 8 records (coefficient " << coeff_idx << " of 8 polys), mismatches " << bad << std::endl;
                std::cout << (bad == 0 ? "│   ► Verification: SUCCESS!" : "│   ► Verification: FAILED!") << std::endl;
            }
        }
        else
        {
            uint64_t index = target;
            RlweSampleQ2 ans[8];
            int64_t* result[8];
            for (uint64_t i = 0; i < 8; i++)
            {
                result[i] = static_cast<int64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));
            }
            for (int k = 0; k<ExecutionCount; k++)
            {
                std::cout << "| \u2756" << " Round " << k + 1 << ": computation time" << std::endl;

                client.Query(index, c_rot, rsk, C_ctrl, C_sel);

                server.Answer(c_rot, rsk, C_ctrl, C_sel, ans);
                server.VerifyLoadedDatabase();
                if (corrupt == "resp")
                {
                    ans[t_target].rlweVal[0] = (ans[t_target].rlweVal[0] + DELTA_FINAL) % MODULUS_Q4;
                    std::cout << "│   ► TEST: corrupted response " << t_target << std::endl;
                }

                client.Recover(ans, result);

                // Response t holds coefficients 4m (m < 512) of poly resp_polys[t].
                uint64_t bad = 0;
                for (uint64_t t = 0; t < 8; t++)
                    for (uint64_t m = 0; m < DEGREE2; m++)
                        if (truth[t].empty() || result[t][m] != static_cast<int64_t>(truth[t][m * (DEGREE1 / DEGREE2)])) bad++;
                std::cout << "│   ► Verification: checked 8x512 coefficients, mismatches " << bad << std::endl;
                std::cout << (bad == 0 ? "│   ► Verification: SUCCESS!" : "│   ► Verification: FAILED!") << std::endl;
            }
        }
    }

    return 0;
}

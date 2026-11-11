#pragma once
#include "VIACrypto.hpp"

class VIAClient
{
public:
	// constructor
	VIAClient()
    {
        temp_poly = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    }
	// destructor
	~VIAClient()
    {
        std::free(temp_poly);
    }

    // Row index layout used by the server's DMux tree: the 8 initial
    // ciphertexts select the 3 most significant row bits (ControlBits[0] is the
    // MSB), each further DMux level appends the next bit as LSB.
    static uint64_t rotSlot(const std::vector<uint8_t>& ControlBits)
    {
        return 4 * ControlBits[0] + 2 * ControlBits[1] + ControlBits[2];
    }

	// Query: record index -> (row, col, k); the response holds coefficients
	// 4m + k (m < 512) of the selected polynomial in 8 columns.
	void Query(uint64_t& index, RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel)
    {
        std::cout << "\u2502   \u25BA" << " Query Time: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t idx_i = (index % (ROW * COL)) / COL;
        uint64_t idx_j = index % COL;
        uint64_t idx_k = index / (ROW * COL);
        buildQuery(idx_i, idx_j, idx_k == 0 ? 0 : (DEGREE1 - idx_k), idx_k != 0, c_rot, rsk, C_ctrl, C_sel);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << static_cast<double>(duration.count())/1000.0 << " ms" << std::endl;
    }

	// Query for blinded extraction: index[0] selects the polynomial, index[1]
	// the coefficient that is rotated into the constant term (X^{-index[1]}).
	void Query_Extraction(uint64_t index[2], RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1> C_ctrl, std::vector<RgswSampleQ2> C_sel)
    {
        std::cout << "\u2502   \u25BA" << " Query Time with Blinded Extraction: ";
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t idx_i = (index[0] % (ROW * COL)) / COL;
        uint64_t idx_j = index[0] % COL;
        uint64_t shift = index[1] % DEGREE1;
        buildQuery(idx_i, idx_j, shift == 0 ? 0 : (DEGREE1 - shift), shift != 0, c_rot, rsk, C_ctrl, C_sel);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << static_cast<double>(duration.count())/1000.0 << " ms" << std::endl;
    }

	// Recover RLWE
	void Recover(RlweSampleQ2 ans[8], int64_t* result[8])
    {
        std::cout << "\u2502   \u25BA" << " Recover Time: ";
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 8; i++) crypto.decryptResponse(ans[i], result[i]);
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
        for (int i = 0; i < 8; i++) result[i] = crypto.decryptResponseLwe(ans[i]);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << duration.count() << " \u03BCs" << std::endl;
        std::cout << "\u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;
    }

    VIACrypto crypto;
    uint64_t* temp_poly;

private:
    // c_rot[slot] = RLWE(+-Delta1 * X^pos), others RLWE(0); C_ctrl = RGSW(row
    // bits 3..LOG_ROW-1); C_sel = RGSW(col bits 0..LOG_COL-4); rsk.
    void buildQuery(uint64_t idx_i, uint64_t idx_j, uint64_t pos, bool negate, RlweSampleQ1 c_rot[8], RingSwitchingKey& rsk, std::vector<RgswSampleQ1>& C_ctrl, std::vector<RgswSampleQ2>& C_sel)
    {
        std::vector<uint8_t> ControlBits(LOG_ROW);
        for (uint64_t i = 0; i < LOG_ROW; i++)
        {
            ControlBits[LOG_ROW - i - 1] = (idx_i >> i) & 1;
        }
        // (was a heap overflow: the blinded variant wrote LOG_COL entries)
        std::vector<uint8_t> SelectBits(LOG_COL - 3);
        for (uint64_t i = 0; i < LOG_COL - 3; i++)
        {
            SelectBits[i] = (idx_j >> i) & 1;
        }
        const uint64_t slot = rotSlot(ControlBits);
        std::fill(temp_poly, temp_poly + DEGREE1, 0);
        for (uint64_t i = 0; i < 8; i++)
        {
            if (i != slot) crypto.encRlweQ1(temp_poly, c_rot[i]);
        }
        temp_poly[pos] = negate ? (MODULUS_Q1 - DELTA1) : DELTA1;
        crypto.encRlweQ1(temp_poly, c_rot[slot]);

        for (uint64_t i = 3; i < LOG_ROW; i++)
        {
            crypto.encRgswQ1Const(ControlBits[i], GADGET1_L, LOG_GADGET1_BASE, C_ctrl[i - 3]);
        }
        for (uint64_t i = 0; i < LOG_COL - 3; i++)
        {
            crypto.encRgswQ2Const(SelectBits[i], GADGET2_L_1, LOG_GADGET2_BASE_1, GADGET2_L_2, LOG_GADGET2_BASE_2, C_sel[i]);
        }
        crypto.genRsk(GADGET_RSK_L, LOG_GADGET_RSK_BASE, rsk);
    }
};

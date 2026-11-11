#include "functions.hpp"

// Discrete Gaussian Distribution

void DiscreteGuassian(float mean, float stddev, uint64_t modulus, uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = static_cast<int64_t>(std::round(normal_dist(gen) * stddev + mean));
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + modulus) : temp_int);
	}	
}

// KGen_1
void KGen_1(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = static_cast<int64_t>(KGen_1_dist(gen));
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q1) : temp_int);
	}
}


// KGen_2
void KGen_2(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = static_cast<int64_t>(KGen_2_dist(gen));
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + 2) : temp_int);
	}
}

// secret key generator 1
void KGen_1C(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = static_cast<int64_t>(KGen_1C_dist(gen));
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q1) : temp_int);
	}
}

// secret key generator 2
void KGen_2C(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = static_cast<int64_t>(KGen_2C_dist(gen));
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + 2) : temp_int);
	}
}

// UniformModQ1_1
void UniformModQ1_1(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = UniformModQ1_1_dist(gen);
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q1_1) : temp_int);
	}
}

// UniformModQ1_2
void UniformModQ1_2(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = UniformModQ1_2_dist(gen);
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q1_2) : temp_int);
	}
}


// uniformModQ2
void UniformModQ2(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = UniformModQ2_dist(gen);
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q2) : temp_int);
	}
}

// uniformModP
void UniformModP(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = UniformModP_dist(gen);
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + HALF_MODULUS_P) : temp_int);
	}
}

//  uniform mod Q1_1C
void UniformModQ1_1C(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = UniformModQ1_1C_dist(gen);
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q1_1) : temp_int);
	}
}

//  uniform mod Q1_2C
void UniformModQ1_2C(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = UniformModQ1_2C_dist(gen);
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q1_2) : temp_int);
	}
}

//  uniform mod Q2C
void UniformModQ2C(uint64_t* result, uint64_t size)
{
	int64_t temp_int;
	for (uint64_t i = 0; i < size; ++i)
	{
		temp_int = UniformModQ2C_dist(gen);
		result[i] = static_cast<uint64_t>((temp_int < 0) ? (temp_int + MODULUS_Q2) : temp_int);
	}
}

// the sum of elements in Src
uint64_t sum(const uint64_t* Src, uint64_t size)
{
	uint64_t res = 0;
	for (uint64_t i = 0; i < size; ++i)
	{
		res += Src[i];
	}
	return res;
}

// solve CRT equations
void CRT(uint64_t& result, const uint64_t& op1, const uint64_t& op2)
{
	result = (op1 * CRT_COEFFICIENT1 + op2 * CRT_COEFFICIENT2) % MODULUS_Q1; 
}

void CRT(uint64_t* result, uint64_t* op1, uint64_t* op2, uint64_t size)
{
	intel::hexl::EltwiseMultMod(op1, op1, CRT_ARRAY1.data(), DEGREE1, MODULUS_Q1, 1);
	intel::hexl::EltwiseMultMod(op2, op2, CRT_ARRAY2.data(), DEGREE1, MODULUS_Q1, 1);
	intel::hexl::EltwiseAddMod(result, op1, op2, DEGREE1, MODULUS_Q1);
}

// solve CRT equations for VIA-C
void CRT_C(uint64_t& result, const uint64_t& op1, const uint64_t& op2)
{
	result = (op1 * CRT_COEFFICIENT1C + op2 * CRT_COEFFICIENT2C) % MODULUS_Q1C; 
}

void CRT_C(uint64_t* result, uint64_t* op1, uint64_t* op2, uint64_t size)
{
	intel::hexl::EltwiseMultMod(op1, op1, CRT_ARRAY1C.data(), DEGREE1, MODULUS_Q1C, 1);
	intel::hexl::EltwiseMultMod(op2, op2, CRT_ARRAY2C.data(), DEGREE2, MODULUS_Q1C, 1);
	intel::hexl::EltwiseAddMod(result, op1, op2, DEGREE1, MODULUS_Q1C);
}

// Gadget Decomposition for poly. of degree1
void GadgetDecompositionDeg1(uint64_t* poly, std::vector<uint64_t*> result, const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base)
{
	const uint64_t HALF_BASE_Pow_LEV = UINT64_C(1) << (log_gadget_base * gadget_l - 1);
	const uint64_t BASE = UINT64_C(1) << log_gadget_base;
	addC_I_Deg1(HALF_BASE_Pow_LEV, poly);
	rshiftC_I_Deg1(poly, log_modulus - log_gadget_base * gadget_l);
	for (uint64_t i = 0; i < gadget_l; ++i)
	{
		andC_Deg1(BASE - 1, poly, result[i]);
		lshiftC_I_Deg1(result[i], 64 - log_gadget_base);
		rshiftC_Signed_I_Deg1(reinterpret_cast<int64_t*>(result[i]), 64 - log_gadget_base);
		addC_I_Deg1(BASE >> 1, poly);
		rshiftC_I_Deg1(poly, log_gadget_base);
	}
}

// Gadget Decomposition for poly. of degree2
void GadgetDecompositionDeg2(uint64_t* poly, std::vector<uint64_t*> result, const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base)
{
	const uint64_t HALF_BASE_Pow_LEV = UINT64_C(1) << (log_gadget_base * gadget_l - 1);
	const uint64_t BASE = UINT64_C(1) << log_gadget_base;
	addC_I_Deg2(HALF_BASE_Pow_LEV, poly);
	rshiftC_I_Deg2(poly, log_modulus - log_gadget_base * gadget_l);
	for (uint64_t i = 0; i < gadget_l; ++i)
	{
		andC_Deg2(BASE - 1, poly, result[i]);
		lshiftC_I_Deg2(result[i], 64 - log_gadget_base);
		rshiftC_Signed_I_Deg2(reinterpret_cast<int64_t*>(result[i]), 64 - log_gadget_base);
		addC_I_Deg2(BASE >> 1, poly);
		rshiftC_I_Deg2(poly, log_gadget_base);
	}
}

// Gadget Decomposition for poly.
void GadgetDecompositionDegAll(uint64_t* poly, std::vector<uint64_t*> result, const uint64_t log_modulus ,const uint64_t gadget_l, const uint64_t log_gadget_base, const uint64_t deg)
{
	const uint64_t HALF_BASE_Pow_LEV = UINT64_C(1) << (log_gadget_base * gadget_l - 1);
	const uint64_t BASE = UINT64_C(1) << log_gadget_base;
	if (deg >= 8)
	{
		addC_I_DegAll(HALF_BASE_Pow_LEV, poly, deg);
		rshiftC_I_DegAll(poly, log_modulus - log_gadget_base * gadget_l, deg);
		for (uint64_t i = 0; i < gadget_l; ++i)
		{
			andC_DegAll(BASE - 1, poly, result[i], deg);
			lshiftC_I_DegAll(result[i], 64 - log_gadget_base, deg);
			rshiftC_Signed_I_DegAll(reinterpret_cast<int64_t*>(result[i]), 64 - log_gadget_base, deg);
			addC_I_DegAll(BASE >> 1, poly, deg);
			rshiftC_I_DegAll(poly, log_gadget_base, deg);
		}
	}
	else
	{
		for (uint64_t j = 0; j < deg; j++)
		{
			poly[j] += HALF_BASE_Pow_LEV;
			poly[j] >>= log_modulus - log_gadget_base * gadget_l;
			for (uint64_t i = 0; i < gadget_l; ++i)
			{
				result[i][j] = poly[j] & (BASE - 1);
				result[i][j] <<= 64 - log_gadget_base;
				result[i][j] >>= 64 - log_gadget_base;
				result[i][j] += BASE >> 1;
				result[i][j] >>= log_gadget_base;
			}
		}
	}
}

// ring embedding
void ring_embed(const uint64_t deg_src, const uint64_t deg_dst, const uint64_t idx, const uint64_t* src, uint64_t* dst)
{
	uint64_t temp = deg_dst / deg_src;
	std::fill(dst, dst + deg_dst, 0);
	for (uint64_t i = 0; i < deg_src; ++i)
	{
		dst[idx + i * temp] = src[i];
	}
}

// AVX512 ring projection loop
template <uint64_t Iter = 0>
inline void RingProj_loop(const uint64_t idx, const uint64_t* src, uint64_t* dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2)
	{  
		constexpr uint64_t offset = Iter * 8;
		__m512i indices = _mm512_setr_epi64(Iter * 32 + idx, 4 +Iter * 32 + idx, 8 + Iter * 32 + idx, 12 + Iter * 32 + idx, 16 + Iter * 32 + idx, 20 + Iter * 32 + idx, 24 + Iter * 32 + idx, 28 + Iter * 32 + idx);

		// Fully unrolled AVX instruction sequence
		_mm512_store_si512((__m512i*)(dst + offset), _mm512_i64gather_epi64(indices, src, 8));

		// Recursive template expansion followed by iteration
		RingProj_loop<Iter + 1>(idx, src, dst);
	}
}	

// ring projection
void ring_project(const uint64_t idx, const uint64_t* src, uint64_t* dst)
{
	RingProj_loop<0>(idx, src, dst);
}

// modulus switching
void modSwitch(const uint64_t modulus_src, const uint64_t modulus_dst, uint64_t* Src, uint64_t* Dst, const uint64_t size)
{
	// round(x * dst / src) mod dst. (The original computed the ratio with
	// integer division, which is 0 whenever dst < src, i.e. always here.)
	const double ratio = static_cast<double>(modulus_dst) / static_cast<double>(modulus_src);
	for (uint64_t i = 0; i < size; ++i)
	{
		uint64_t r = static_cast<uint64_t>(static_cast<double>(Src[i]) * ratio + 0.5);
		Dst[i] = (r >= modulus_dst) ? r - modulus_dst : r;
	}
}

// printInfoVIA
void printInfoVIA(uint64_t IsBlindedExtraction)
{
    std::cout << "\u250C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510" << std::endl;
	std::cout << "\u2502             Test for VIA Protocol             \u2502" << std::endl;
    std::cout << "\u251C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518" << std::endl;
	std::cout << "\u2502 Database parameters:" << std::endl;
    std::cout << "\u2502   \u25BA" << " Length of Database(N) : " << ROW*COL*DEGREE1/DEGREE2 << std::endl;
    std::cout << "\u2502   \u25BA" << " Row: " << ROW << std::endl;
    std::cout << "\u2502   \u25BA" << " Col: " << COL << std::endl;
	std::cout << "\u2502   \u25BA" << " Database size: " << ((ROW*COL*DEGREE1*LOG_MODULUS_P) >> 23) << " MB" << std::endl;
	std::cout << "\u2502   " << std::endl;
    std::cout << "\u2502 LWE parameters(modulus):" << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q1: " << MODULUS_Q1 << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q1_1: " << MODULUS_Q1_1 << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q1_2: " << MODULUS_Q1_2 << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q2: " << MODULUS_Q2 << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q3: " << MODULUS_Q3 << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q4: " << MODULUS_Q4 << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus P: " << MODULUS_P << std::endl;

    std::cout << "\u2502 LWE parameters(degree):" << std::endl;
    std::cout << "\u2502   \u25BA" << " Degree of RLWE1: " << DEGREE1 << std::endl;
    std::cout << "\u2502   \u25BA" << " Degree of RLWE2: " << DEGREE2 << std::endl;

    std::cout << "\u2502 Gadget parameters:" << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget1_L: " << GADGET1_L << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget1_Base: " << GADGET1_BASE << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget2_L_1: " << GADGET2_L_1 << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget2_Base_1: " << GADGET2_BASE_1 << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget2_L_2: " << GADGET2_L_2 << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget2_Base_2: " << GADGET2_BASE_2 << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget_RSK_L: " << GADGET_RSK_L << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget_RSK_Base: " << GADGET_RSK_BASE << std::endl;
    std::cout << "\u251C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;


	double OnlineUpload = static_cast<double>(2*GADGET1_L*(LOG_ROW-3)+8)*sizePolyDeg1Q1 + static_cast<double>(GADGET_RSK_L)*sizePolyDeg1Q2 + static_cast<double>((GADGET2_L_1+GADGET2_L_2)*(LOG_COL-3)*sizePolyDeg2Q2);
	double OnlineDownload = (IsBlindedExtraction == 1) ? 8*sizePolyDeg2Q3 + LOG_MODULUS_Q4/1024.0 : 8*(sizePolyDeg2Q3 + sizePolyDeg2Q4);
	std::cout << "\u2502 Communication overhead:" << std::endl;
	std::cout << "\u2502   \u25BA" << " Offline upload: " << 0 << std::endl;
	std::cout << "\u2502   \u25BA" << " Offline download: " << 0 << std::endl;
	std::cout << "\u2502   \u25BA" << " Online upload: " << OnlineUpload <<  " KiB" <<std::endl;
	std::cout << "\u2502   \u25BA" << " Online download: " << OnlineDownload << " KiB" << std::endl;
    std::cout << "\u251C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;

}

// printInfoVIA_C
void printInfoVIA_C(uint64_t IsBlindedExtraction)
{
    std::cout << "\u250C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510" << std::endl;
	std::cout << "\u2502            Test for VIA-C Protocol            \u2502" << std::endl;
    std::cout << "\u251C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518" << std::endl;
	std::cout << "\u2502 Database parameters:" << std::endl;
    std::cout << "\u2502   \u25BA" << " Length of Database(N) : " << ROW*COL*DEGREE1/DEGREE2 << std::endl;
    std::cout << "\u2502   \u25BA" << " Row: " << ROW << std::endl;
    std::cout << "\u2502   \u25BA" << " Col: " << COL << std::endl;
	std::cout << "\u2502   \u25BA" << " Database size: " << ((ROW*COL*DEGREE1*LOG_MODULUS_P) >> 24) << " MB" << std::endl;
	std::cout << "\u2502   " << std::endl;
    std::cout << "\u2502 LWE parameters(modulus):" << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q1: " << MODULUS_Q1C << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q1_1: " << MODULUS_Q1_1C << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q1_2: " << MODULUS_Q1_2C << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q2: " << MODULUS_Q2C << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q3: " << MODULUS_Q3C << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus Q4: " << MODULUS_Q4C << std::endl;
    std::cout << "\u2502   \u25BA" << " Modulus P: " << MODULUS_PC << std::endl;

    std::cout << "\u2502 LWE parameters(degree):" << std::endl;
    std::cout << "\u2502   \u25BA" << " Degree of LWE: " << DEGREE1 << std::endl;
    std::cout << "\u2502   \u25BA" << " Degree of RLWE1: " << DEGREE1 << std::endl;
    std::cout << "\u2502   \u25BA" << " Degree of RLWE2: " << DEGREE2 << std::endl;

    std::cout << "\u2502 Gadget parameters:" << std::endl;
	std::cout << "\u2502   \u25BA" << " Gadget_L_Conv: " << GADGET_L_CONV << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget_Base_Conv: " << GADGET_BASE_CONV << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget1_L: " << GADGET1_LC << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget1_Base: " << GADGET1_BASEC << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget2_L: " << GADGET2_LC << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget2_Base: " << GADGET2_BASEC << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget_RSK_L: " << GADGET_RSK_L << std::endl;
    std::cout << "\u2502   \u25BA" << " Gadget_RSK_Base: " << GADGET_RSK_BASE << std::endl;
    std::cout << "\u251C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;



	double OfflineUpload = 22*2*GADGET_L_CONV*sizePolyDeg1Q1C + GADGET_L_CONV*sizePolyDeg1Q1C + GADGET_RSK_LC*sizePolyDeg2Q3C;

	double OnlineUpload = (IsBlindedExtraction == 1) ? GADGET1_LC*(LOG_ROW + LOG_COL + 11)*CEIL_LOG_Q1C / 8192.0 : GADGET1_LC*(LOG_ROW + LOG_COL + 2)*CEIL_LOG_Q1C / 8192.0;
	double OnlineDownload = (IsBlindedExtraction == 1) ?  (DEGREE2 * CEIL_LOG_Q3C + LOG_MODULUS_Q4C) / 8192.0 : DEGREE2 * (CEIL_LOG_Q3C + LOG_MODULUS_Q4C) /8192.0;
	std::cout << "\u2502 Communication overhead:" << std::endl;
	std::cout << "\u2502   \u25BA" << " Offline upload: " << OfflineUpload / 1024.0 << " MiB" << std::endl;
	std::cout << "\u2502   \u25BA" << " Offline download: " << 0 << std::endl;
	std::cout << "\u2502   \u25BA" << " Online upload: " << OnlineUpload <<  " KiB" <<std::endl;
	std::cout << "\u2502   \u25BA" << " Online download: " << OnlineDownload << " KiB" << std::endl;
    std::cout << "\u251C\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500" << std::endl;

}


// AVX512 multiply Constant loop
template <uint64_t Iter = 0>
inline void Deg1_MulC_loop(const uint64_t val, uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1)
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_m = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_mullo_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), vval_m);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_MulC_loop<Iter + 1>(val, Src, Dst);
	}
}


// AVX512 multiply Constant
void mulC_Deg1(const uint64_t val, uint64_t* Src, uint64_t* Dst)
{
	Deg1_MulC_loop<0>(val, Src, Dst);
}

// AVX512 In-place multiply Constant
template <uint64_t Iter = 0>
inline void Deg1_MulC_I_loop(const uint64_t val, uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_mul = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_mullo_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), vval_mul);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_MulC_I_loop<Iter + 1>(val, SrcDst);
	}
}

// AVX512 In-place multiply Constant
void mulC_I_Deg1(const uint64_t val, uint64_t* SrcDst)
{
	Deg1_MulC_I_loop<0>(val, SrcDst);
}

// AVX512 addition loop
template <uint64_t Iter = 0>
inline void Deg1_add_loop(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v1 = _mm512_load_si512(
			(__m512i*)(Src1 + offset));
		__m512i v2 = _mm512_load_si512(
			(__m512i*)(Src2 + offset));
		__m512i vdst = _mm512_add_epi64(v1, v2);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vdst);

		// Recursive template expansion followed by iteration
		Deg1_add_loop<Iter + 1>(Src1, Src2, Dst);
	}
}


// AVX512 addition
void add_Deg1(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	Deg1_add_loop<0>(Src1, Src2, Dst);
}

// AVX512 In-place addition loop
template <uint64_t Iter = 0>
inline void Deg1_add_I_loop(const uint64_t* Src, uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vs = _mm512_load_si512(
			(__m512i*)(Src + offset));
		__m512i vd = _mm512_load_si512(
			(__m512i*)(SrcDst + offset));
		vs = _mm512_add_epi64(vs, vd);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vs);

		// Recursive template expansion followed by iteration
		Deg1_add_I_loop<Iter + 1>(Src, SrcDst);
	}
}


// AVX512 In-place addition
void add_I_Deg1(const uint64_t* Src, uint64_t* SrcDst)
{
	Deg1_add_I_loop<0>(Src, SrcDst);
}

// AVX512 subtraction loop
template <uint64_t Iter = 0>
inline void Deg1_sub_loop(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v1 = _mm512_load_si512(
			(__m512i*)(Src1 + offset));
		__m512i v2 = _mm512_load_si512(
			(__m512i*)(Src2 + offset));
		__m512i vd = _mm512_sub_epi64(v1, v2);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_sub_loop<Iter + 1>(Src1, Src2, Dst);
	}
}

// AVX512 subtraction
void sub_Deg1(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	Deg1_sub_loop<0>(Src1, Src2, Dst);
}

// AVX512 In-place subtraction loop
template <uint64_t Iter = 0>
inline void Deg1_sub_I_loop(uint64_t* SrcDst, const uint64_t* Src)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vs = _mm512_load_si512(
			(__m512i*)(Src + offset));
		__m512i vd = _mm512_load_si512(
			(__m512i*)(SrcDst + offset));
		vs = _mm512_sub_epi64(vd, vs);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vs);

		// Recursive template expansion followed by iteration
		Deg1_sub_I_loop<Iter + 1>(SrcDst, Src);
	}
}

// AVX512 In-place subtraction
void sub_I_Deg1(uint64_t* SrcDst, const uint64_t* Src)
{
	Deg1_sub_I_loop<0>(SrcDst, Src);
}

// AVX512 negative loop
template <uint64_t Iter = 0>
inline void Deg1_negative_loop(uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{ 
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v = _mm512_load_si512(
			(__m512i*)(Src + offset));
		v = _mm512_sub_epi64(_mm512_setzero_si512(), v);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), v);

		// Recursive template expansion followed by iteration
		Deg1_negative_loop<Iter + 1>(Src, Dst);
	}
}

// AVX512 negative
void negative_Deg1(uint64_t* Src, uint64_t*Dst)
{
	Deg1_negative_loop<0>(Src, Dst);
}

// AVX512 In-place negative loop
template <uint64_t Iter = 0>
inline void Deg1_negative_I_loop(uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v = _mm512_load_si512(
			(__m512i*)(SrcDst + offset));
		v = _mm512_sub_epi64(_mm512_setzero_si512(), v);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), v);

		// Recursive template expansion followed by iteration
		Deg1_negative_I_loop<Iter + 1>(SrcDst);
	}
}

// AVX512 In-place negative
void negative_I_Deg1(uint64_t* SrcDst)
{
	Deg1_negative_I_loop<0>(SrcDst);
}

// AVX512 copy loop
template <uint64_t Iter = 0>
inline void Deg1_copy_loop(const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		_mm512_store_si512((__m512i*)(Dst + offset), _mm512_load_si512((__m512i*)(Src + offset)));

		// Recursive template expansion followed by iteration
		Deg1_copy_loop<Iter + 1>(Src, Dst);
	}
}

// AVX512 copy
void copy_Deg1(const uint64_t* Src, uint64_t* Dst)
{
	Deg1_copy_loop<0>(Src, Dst);
}

// AVX512 add Constant loop
template <uint64_t Iter = 0>
inline void Deg1_AddC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_addC = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_add_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), vval_addC);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_AddC_loop<Iter + 1>(val, Src, Dst);
	}
}

// AVX512 add Constant
void addC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg1_AddC_loop<0>(val, Src, Dst);
}

// AVX512 In-place add Constant loop
template <uint64_t Iter = 0>
inline void Deg1_AddC_I_loop(const uint64_t val, uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_add_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), vval);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_AddC_I_loop<Iter + 1>(val, SrcDst);
	}
}

// AVX512 In-place add Constant
void addC_I_Deg1(const uint64_t val, uint64_t* SrcDst)
{
	Deg1_AddC_I_loop<0>(val, SrcDst);
}

// AVX512 shift right loop
template <uint64_t Iter = 0>
inline void Deg1_rshiftC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_srli_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), val);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_rshiftC_loop<Iter + 1>(val, Src, Dst);
	}
}


// AVX512 shift right
void rshiftC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg1_rshiftC_loop<0>(val, Src, Dst);
}


// AVX512 In-place shift right loop
template <uint64_t Iter = 0>
inline void Deg1_rshiftC_I_loop(uint64_t* SrcDst, const uint64_t val)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1)
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_srli_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), val);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_rshiftC_I_loop<Iter + 1>(SrcDst, val);
	}
}

// AVX512 In-place shift right
void rshiftC_I_Deg1(uint64_t* SrcDst, uint64_t shift)
{
	Deg1_rshiftC_I_loop<0>(SrcDst, shift);
}

// AVX512 In-place shift right signed loop
template <uint64_t Iter = 0>
inline void Deg1_rshiftC_Signed_I_loop(int64_t* SrcDst, const uint64_t val)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1)
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_srai_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), val);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_rshiftC_Signed_I_loop<Iter + 1>(SrcDst, val);
	}
}

// AVX512 In-place shift right signed
void rshiftC_Signed_I_Deg1(int64_t* SrcDst, uint64_t shift)
{
	Deg1_rshiftC_Signed_I_loop<0>(SrcDst, shift);
}


// AVX512 In-place shift left loop
template <uint64_t Iter = 0>
inline void Deg1_lshiftC_I_loop(uint64_t* SrcDst, const uint64_t val)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1)
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_slli_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), val);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_lshiftC_I_loop<Iter + 1>(SrcDst, val);
	}
}

// AVX512 shift left loop
template <uint64_t Iter = 0>
inline void Deg1_lshiftC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_slli_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), val);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_lshiftC_loop<Iter + 1>(val, Src, Dst);
	}
}


// AVX512 shift left
void lshiftC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg1_lshiftC_loop<0>(val, Src, Dst);
}

// AVX512 In-place shift left
void lshiftC_I_Deg1(uint64_t* SrcDst, uint64_t shift)
{
	Deg1_lshiftC_I_loop<0>(SrcDst, shift);
}

// AVX512 and Constant loop
template <uint64_t Iter = 0>
inline void Deg1_andC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg1) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_andC = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_and_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), vval_andC);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg1_andC_loop<Iter + 1>(val, Src, Dst);
	}
}

// AVX512 logic and constant
void andC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg1_andC_loop<0>(val, Src, Dst);
}


// AVX512 multiply Constant loop
template <uint64_t Iter = 0>
inline void Deg2_MulC_loop(const uint64_t val, uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2)
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_m = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_mullo_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), vval_m);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_MulC_loop<Iter + 1>(val, Src, Dst);
	}
}


// AVX512 multiply Constant
void mulC_Deg2(const uint64_t val, uint64_t* Src, uint64_t* Dst)
{
	Deg2_MulC_loop<0>(val, Src, Dst);
}

// AVX512 In-place multiply Constant
template <uint64_t Iter = 0>
inline void Deg2_MulC_I_loop(const uint64_t val, uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_mul = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_mullo_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), vval_mul);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_MulC_I_loop<Iter + 1>(val, SrcDst);
	}
}

// AVX512 In-place multiply Constant
void mulC_I_Deg2(const uint64_t val, uint64_t* SrcDst)
{
	Deg2_MulC_I_loop<0>(val, SrcDst);
}

// AVX512 addition loop
template <uint64_t Iter = 0>
inline void Deg2_add_loop(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v1 = _mm512_load_si512(
			(__m512i*)(Src1 + offset));
		__m512i v2 = _mm512_load_si512(
			(__m512i*)(Src2 + offset));
		__m512i vdst = _mm512_add_epi64(v1, v2);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vdst);

		// Recursive template expansion followed by iteration
		Deg2_add_loop<Iter + 1>(Src1, Src2, Dst);
	}
}


// AVX512 addition
void add_Deg2(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	Deg2_add_loop<0>(Src1, Src2, Dst);
}

// AVX512 In-place addition loop
template <uint64_t Iter = 0>
inline void Deg2_add_I_loop(const uint64_t* Src, uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vs = _mm512_load_si512(
			(__m512i*)(Src + offset));
		__m512i vd = _mm512_load_si512(
			(__m512i*)(SrcDst + offset));
		vs = _mm512_add_epi64(vs, vd);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vs);

		// Recursive template expansion followed by iteration
		Deg2_add_I_loop<Iter + 1>(Src, SrcDst);
	}
}


// AVX512 In-place addition
void add_I_Deg2(const uint64_t* Src, uint64_t* SrcDst)
{
	Deg2_add_I_loop<0>(Src, SrcDst);
}

// AVX512 subtraction loop
template <uint64_t Iter = 0>
inline void Deg2_sub_loop(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v1 = _mm512_load_si512(
			(__m512i*)(Src1 + offset));
		__m512i v2 = _mm512_load_si512(
			(__m512i*)(Src2 + offset));
		__m512i vd = _mm512_sub_epi64(v1, v2);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_sub_loop<Iter + 1>(Src1, Src2, Dst);
	}
}

// AVX512 subtraction
void sub_Deg2(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst)
{
	Deg2_sub_loop<0>(Src1, Src2, Dst);
}

// AVX512 In-place subtraction loop
template <uint64_t Iter = 0>
inline void Deg2_sub_I_loop(uint64_t* SrcDst, const uint64_t* Src)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vs = _mm512_load_si512(
			(__m512i*)(Src + offset));
		__m512i vd = _mm512_load_si512(
			(__m512i*)(SrcDst + offset));
		vs = _mm512_sub_epi64(vd, vs);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vs);

		// Recursive template expansion followed by iteration
		Deg2_sub_I_loop<Iter + 1>(SrcDst, Src);
	}
}

// AVX512 In-place subtraction
void sub_I_Deg2(uint64_t* SrcDst, const uint64_t* Src)
{
	Deg2_sub_I_loop<0>(SrcDst, Src);
}

// AVX512 negative loop
template <uint64_t Iter = 0>
inline void Deg2_negative_loop(uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{ 
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v = _mm512_load_si512(
			(__m512i*)(Src + offset));
		v = _mm512_sub_epi64(_mm512_setzero_si512(), v);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), v);

		// Recursive template expansion followed by iteration
		Deg2_negative_loop<Iter + 1>(Src, Dst);
	}
}

// AVX512 negative
void negative_Deg2(uint64_t* Src, uint64_t*Dst)
{
	Deg2_negative_loop<0>(Src, Dst);
}

// AVX512 In-place negative loop
template <uint64_t Iter = 0>
inline void Deg2_negative_I_loop(uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i v = _mm512_load_si512(
			(__m512i*)(SrcDst + offset));
		v = _mm512_sub_epi64(_mm512_setzero_si512(), v);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), v);

		// Recursive template expansion followed by iteration
		Deg2_negative_I_loop<Iter + 1>(SrcDst);
	}
}

// AVX512 In-place negative
void negative_I_Deg2(uint64_t* SrcDst)
{
	Deg2_negative_I_loop<0>(SrcDst);
}

// AVX512 copy loop
template <uint64_t Iter = 0>
inline void Deg2_copy_loop(const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		_mm512_store_si512((__m512i*)(Dst + offset), _mm512_load_si512((__m512i*)(Src + offset)));

		// Recursive template expansion followed by iteration
		Deg1_copy_loop<Iter + 1>(Src, Dst);
	}
}

// AVX512 copy
void copy_Deg2(const uint64_t* Src, uint64_t* Dst)
{
	Deg1_copy_loop<0>(Src, Dst);
}

// AVX512 add Constant loop
template <uint64_t Iter = 0>
inline void Deg2_AddC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_addC = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_add_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), vval_addC);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_AddC_loop<Iter + 1>(val, Src, Dst);
	}
}

// AVX512 add Constant
void addC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg2_AddC_loop<0>(val, Src, Dst);
}

// AVX512 In-place add Constant loop
template <uint64_t Iter = 0>
inline void Deg2_AddC_I_loop(const uint64_t val, uint64_t* SrcDst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_add_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), vval);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_AddC_I_loop<Iter + 1>(val, SrcDst);
	}
}

// AVX512 In-place add Constant
void addC_I_Deg2(const uint64_t val, uint64_t* SrcDst)
{
	Deg2_AddC_I_loop<0>(val, SrcDst);
}

// AVX512 shift right loop
template <uint64_t Iter = 0>
inline void Deg2_rshiftC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_srli_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), val);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_rshiftC_loop<Iter + 1>(val, Src, Dst);
	}
}


// AVX512 shift right
void rshiftC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg2_rshiftC_loop<0>(val, Src, Dst);
}

// AVX512 In-place shift right loop
template <uint64_t Iter = 0>
inline void Deg2_rshiftC_I_loop(uint64_t* SrcDst, const uint64_t val)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2)
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_srli_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), val);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_rshiftC_I_loop<Iter + 1>(SrcDst, val);
	}
}

// AVX512 In-place shift right
void rshiftC_I_Deg2(uint64_t* SrcDst, uint64_t shift)
{
	Deg2_rshiftC_I_loop<0>(SrcDst, shift);
}

// AVX512 In-place shift right signed loop
template <uint64_t Iter = 0>
inline void Deg2_rshiftC_Signed_I_loop(int64_t* SrcDst, const uint64_t val)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2)
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_srai_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), val);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_rshiftC_Signed_I_loop<Iter + 1>(SrcDst, val);
	}
}

// AVX512 In-place shift right signed
void rshiftC_Signed_I_Deg2(int64_t* SrcDst, uint64_t shift)
{
	Deg2_rshiftC_Signed_I_loop<0>(SrcDst, shift);
}

// AVX512 shift left loop
template <uint64_t Iter = 0>
inline void Deg2_lshiftC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_slli_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), val);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_lshiftC_loop<Iter + 1>(val, Src, Dst);
	}
}


// AVX512 shift left
void lshiftC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg2_lshiftC_loop<0>(val, Src, Dst);
}

// AVX512 In-place shift left loop
template <uint64_t Iter = 0>
inline void Deg2_lshiftC_I_loop(uint64_t* SrcDst, const uint64_t val)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2)
	{  
		constexpr uint64_t offset = Iter * 8;

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_slli_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + offset)), val);
		_mm512_store_si512(
			(__m512i*)(SrcDst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_lshiftC_I_loop<Iter + 1>(SrcDst, val);
	}
}

// AVX512 In-place shift left
void lshiftC_I_Deg2(uint64_t* SrcDst, uint64_t shift)
{
	Deg2_lshiftC_I_loop<0>(SrcDst, shift);
}

// AVX512 and Constant loop
template <uint64_t Iter = 0>
inline void Deg2_andC_loop(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	if constexpr (Iter < AVX512_LOOP_SIZE_Deg2) 
	{  
		constexpr uint64_t offset = Iter * 8;
		const __m512i vval_andC = _mm512_set1_epi64(val);

		// Fully unrolled AVX instruction sequence
		__m512i vd = _mm512_and_epi64(_mm512_load_si512(
			(__m512i*)(Src + offset)), vval_andC);
		_mm512_store_si512(
			(__m512i*)(Dst + offset), vd);

		// Recursive template expansion followed by iteration
		Deg2_andC_loop<Iter + 1>(val, Src, Dst);
	}
}

// AVX512 logic and constant
void andC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst)
{
	Deg2_andC_loop<0>(val, Src, Dst);
}

//--------------------------------------
// AVX512 multiply Constant
void mulC_DegAll(const uint64_t val, uint64_t* Src, uint64_t* Dst, const uint64_t size)
{
	const __m512i vval_m = _mm512_set1_epi64(val);
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_mullo_epi64(_mm512_load_si512(
			(__m512i*)(Src + 8*i)), vval_m);
		_mm512_store_si512(
			(__m512i*)(Dst + 8*i), vd);
	}
}

// AVX512 In-place multiply Constant
void mulC_I_DegAll(const uint64_t val, uint64_t* SrcDst, const uint64_t size)
{
	const __m512i vval_mul = _mm512_set1_epi64(val);
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_mullo_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + 8*i)), vval_mul);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), vd);
	}
}

// AVX512 addition
void add_DegAll(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i v1 = _mm512_load_si512(
			(__m512i*)(Src1 + 8*i));
		__m512i v2 = _mm512_load_si512(
			(__m512i*)(Src2 + 8*i));
		__m512i vdst = _mm512_add_epi64(v1, v2);
		_mm512_store_si512(
			(__m512i*)(Dst + 8*i), vdst);
	}
}

// AVX512 In-place addition
void add_I_DegAll(const uint64_t* Src, uint64_t* SrcDst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vs = _mm512_load_si512(
			(__m512i*)(Src + 8*i));
		__m512i vd = _mm512_load_si512(
			(__m512i*)(SrcDst + 8*i));
		vs = _mm512_add_epi64(vs, vd);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), vs);
	}
}

// AVX512 subtraction
void sub_DegAll(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i v1 = _mm512_load_si512(
			(__m512i*)(Src1 + 8*i));
		__m512i v2 = _mm512_load_si512(
			(__m512i*)(Src2 + 8*i));
		__m512i vd = _mm512_sub_epi64(v1, v2);
		_mm512_store_si512(
			(__m512i*)(Dst + 8*i), vd);
	}
}

// AVX512 In-place subtraction
void sub_I_DegAll(uint64_t* SrcDst, const uint64_t* Src, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vs = _mm512_load_si512(
			(__m512i*)(Src + 8*i));
		__m512i vd = _mm512_load_si512(
			(__m512i*)(SrcDst + 8*i));
		vs = _mm512_sub_epi64(vd, vs);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), vs);
	}
}

// AVX512 negative
void negative_DegAll(uint64_t* Src, uint64_t*Dst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i v = _mm512_load_si512(
		(__m512i*)(Src + 8*i));
		v = _mm512_sub_epi64(_mm512_setzero_si512(), v);
		_mm512_store_si512(
		(__m512i*)(Dst + 8*i), v);
	}
}

// AVX512 In-place negative
void negative_I_DegAll(uint64_t* SrcDst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i v = _mm512_load_si512(
			(__m512i*)(SrcDst + 8*i));
		v = _mm512_sub_epi64(_mm512_setzero_si512(), v);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), v);
	}
}

// AVX512 copy
void copy_DegAll(const uint64_t* Src, uint64_t* Dst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		_mm512_store_si512((__m512i*)(Dst + 8*i), _mm512_load_si512((__m512i*)(Src + 8*i)));
	}
}

// AVX512 add Constant
void addC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size)
{
	const __m512i vval_addC = _mm512_set1_epi64(val);
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_add_epi64(_mm512_load_si512(
			(__m512i*)(Src + 8*i)), vval_addC);
		_mm512_store_si512(
			(__m512i*)(Dst + 8*i), vd);
	}
}

// AVX512 In-place add Constant
void addC_I_DegAll(const uint64_t val, uint64_t* SrcDst, const uint64_t size)
{
	const __m512i vval = _mm512_set1_epi64(val);
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_add_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + 8*i)), vval);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), vd);
	}
}

// AVX512 shift right
void rshiftC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_srli_epi64(_mm512_load_si512(
			(__m512i*)(Src + 8*i)), val);
		_mm512_store_si512(
			(__m512i*)(Dst + 8*i), vd);
	}
}

// AVX512 In-place shift right
void rshiftC_I_DegAll(uint64_t* SrcDst, uint64_t shift, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_srli_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + 8*i)), shift);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), vd);
	}
}

// AVX512 In-place shift right signed
void rshiftC_Signed_I_DegAll(int64_t* SrcDst, uint64_t shift, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_srai_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + 8*i)), shift);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), vd);
	}
}

// AVX512 shift left
void lshiftC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_slli_epi64(_mm512_load_si512(
			(__m512i*)(Src + 8*i)), val);
		_mm512_store_si512(
			(__m512i*)(Dst + 8*i), vd);
	}
}

// AVX512 In-place shift left
void lshiftC_I_DegAll(uint64_t* SrcDst, uint64_t shift, const uint64_t size)
{
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_slli_epi64(_mm512_load_si512(
			(__m512i*)(SrcDst + 8*i)), shift);
		_mm512_store_si512(
			(__m512i*)(SrcDst + 8*i), vd);
	}
}

// AVX512 logic and constant
void andC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size)
{
	const __m512i vval_andC = _mm512_set1_epi64(val);
	for (uint64_t i = 0; i < size >> 3; i++)
	{
		__m512i vd = _mm512_and_epi64(_mm512_load_si512(
			(__m512i*)(Src + 8*i)), vval_andC);
		_mm512_store_si512(
			(__m512i*)(Dst + 8*i), vd);
	}
}

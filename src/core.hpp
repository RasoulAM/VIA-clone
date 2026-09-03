#pragma once
#include <hexl/hexl.hpp>
#include <iostream>
#include <cstdint>
#include <stdexcept>
#include <thread>
#include <algorithm>
#include <array>
#include <vector>
#include <chrono>
#include <cmath>

template<typename T, T val, std::size_t... I>
constexpr std::array<T, sizeof...(I)> make_array_impl(std::index_sequence<I...>)
{
    return { ((void)I, val)... };
}

template<typename T, T val, std::size_t N>
constexpr std::array<T, N> make_array()
{
    return make_array_impl<T, val>(std::make_index_sequence<N>{});
}

// parameters for database size
extern uint64_t LOG_ROW;
extern uint64_t LOG_COL;

// parameters for modulus
constexpr static const uint64_t MODULUS_Q1_1 = 268369921;
constexpr static const uint64_t MODULUS_Q1_2 = 536608769;
constexpr static const uint64_t MODULUS_Q1 = MODULUS_Q1_1 * MODULUS_Q1_2;
constexpr static const uint64_t CEIL_LOG_Q1 = 57;
// using in CRT
constexpr static const uint64_t INV1 = 13834852; // inverse of Q1_2 modulo Q1_1
constexpr static const uint64_t INV2 = 508945822; // inverse of Q1_1 modulo Q1_2

constexpr static const uint64_t MODULUS_Q2 = 34359214081;
constexpr static const uint64_t CEIL_LOG_Q2 = 35;

constexpr static const uint64_t LOG_MODULUS_Q3 = 31;
constexpr static const uint64_t MODULUS_Q3 = 2147352577;
constexpr static const uint64_t LOG_MODULUS_Q4 = 15;
constexpr static const uint64_t MODULUS_Q4 = UINT64_C(1) << LOG_MODULUS_Q4;

constexpr static const uint64_t LOG_MODULUS_P = 8;
constexpr static const uint64_t MODULUS_P = 1 << LOG_MODULUS_P;

// parameters for modulus (VIA-C)
constexpr static const uint64_t MODULUS_Q1_1C = 137438822401;
constexpr static const uint64_t MODULUS_Q1_2C = 274810798081;
constexpr static const uint64_t MODULUS_Q1C = MODULUS_Q1_1C * MODULUS_Q1_2C;
constexpr static const uint64_t CEIL_LOG_Q1C = 75;
// using in CRT
constexpr static const uint64_t INV1C = 111955861741; // inverse of Q1_2 modulo Q1_1
constexpr static const uint64_t INV2C = 50953527061; // inverse of Q1_1 modulo Q1_2

constexpr static const uint64_t MODULUS_Q2C = 17175674881;
constexpr static const uint64_t CEIL_LOG_Q2C = 34;

constexpr static const uint64_t MODULUS_Q3C = 8380417;
constexpr static const uint64_t CEIL_LOG_Q3C = 23;
constexpr static const uint64_t LOG_MODULUS_Q4C = 12;
constexpr static const uint64_t MODULUS_Q4C = UINT64_C(1) << LOG_MODULUS_Q4C;

constexpr static const uint64_t LOG_MODULUS_PC = 8;
constexpr static const uint64_t MODULUS_PC = 1 << LOG_MODULUS_PC;

// parameters for LWE degree
constexpr static const uint64_t DEGREE1 = 2048;
constexpr static const uint64_t DEGREE2 = 512;
constexpr static const uint64_t DEGREE_LWE1 = 2048;
constexpr static const uint64_t DEGREE_LWE2 = 512;

// parameters for secret key generation
constexpr static const int64_t MODULUS_KGEN1 = 3;
constexpr static const int64_t MODULUS_KGEN2 = 759;
constexpr static const int64_t MODULUS_KGEN1C = 3;
constexpr static const int64_t MODULUS_KGEN2C = 759;

// parameters for error distribution
constexpr static const double STDDEVLWE1 = 1.0;
constexpr static const double STDDEVLWE2 = 1.0;
constexpr static const double STDDEV1 = 1.0;
constexpr static const double STDDEV2 = 1.0;
constexpr static const double STDDEV1C = 1.0;
constexpr static const double STDDEV2C = 1.0;

// gadget parameters
constexpr static const uint64_t GADGET1_L = 2;  
constexpr static const uint64_t LOG_GADGET1_BASE = 19;
constexpr static const uint64_t GADGET1_BASE = 370758;

constexpr static const uint64_t GADGET2_L_1 = 4;
constexpr static const uint64_t LOG_GADGET2_BASE_1 = 5;
constexpr static const uint64_t GADGET2_BASE_1 = 24;

constexpr static const uint64_t GADGET2_L_2 = 3;
constexpr static const uint64_t LOG_GADGET2_BASE_2 = 5;
constexpr static const uint64_t GADGET2_BASE_2 = 24;

constexpr static const uint64_t GADGET_RSK_L = 4;
constexpr static const uint64_t LOG_GADGET_RSK_BASE = 5;
constexpr static const uint64_t GADGET_RSK_BASE = 24;

constexpr static const uint64_t GADGET_L_CONV = 18;  
constexpr static const uint64_t LOG_GADGET_BASE_CONV = 5;
constexpr static const uint64_t GADGET_BASE_CONV = 18;

constexpr static const uint64_t GADGET1_LC = 2;  
constexpr static const uint64_t LOG_GADGET1_BASEC = 16;
constexpr static const uint64_t GADGET1_BASEC = 55879;

constexpr static const uint64_t GADGET2_LC = 2;
constexpr static const uint64_t LOG_GADGET2_BASEC = 7;
constexpr static const uint64_t GADGET2_BASEC = 81;

constexpr static const uint64_t GADGET_RSK_LC = 8;
constexpr static const uint64_t LOG_GADGET_RSK_BASEC = 3;
constexpr static const uint64_t GADGET_RSK_BASEC = 8;

constexpr static const uint64_t GADGET_L_MAX = 4;
constexpr static const uint64_t GADGET_LC_MAX = 18;

// using in pir
extern uint64_t ROW;
extern uint64_t COL;

constexpr static const uint64_t DELTA1 = MODULUS_Q1 / MODULUS_P;
constexpr static const uint64_t DELTA_FINAL = MODULUS_Q4 / MODULUS_P;
constexpr static const uint64_t DELTAC_FINAL = MODULUS_Q4C / MODULUS_P;


constexpr static const int64_t HALF_MODULUS_KGEN1 = (MODULUS_KGEN1 - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_KGEN2 = (MODULUS_KGEN2 - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_KGEN1C = (MODULUS_KGEN1C - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_KGEN2C = (MODULUS_KGEN2C - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q1 = (MODULUS_Q1 - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q1_1 = (MODULUS_Q1_1 - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q1_2 = (MODULUS_Q1_2 - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q2 = (MODULUS_Q2 - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q4 = MODULUS_Q4 >> 1;
constexpr static const int64_t HALF_MODULUS_P = MODULUS_P >> 1;
constexpr static const int64_t HALF_MODULUS_Q1C = (MODULUS_Q1C - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q1_1C = (MODULUS_Q1_1C - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q1_2C = (MODULUS_Q1_2C - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q2C = (MODULUS_Q2 - 1) >> 1;
constexpr static const int64_t HALF_MODULUS_Q4C = MODULUS_Q4C >> 1;

constexpr static const uint64_t CRT_COEFFICIENT1 = (MODULUS_Q1_2 * INV1) % MODULUS_Q1;
constexpr static const uint64_t CRT_COEFFICIENT2 = (MODULUS_Q1_1 * INV2) % MODULUS_Q1;
constexpr static const uint64_t CRT_COEFFICIENT1C = (MODULUS_Q1_2C * INV1C) % MODULUS_Q1C;
constexpr static const uint64_t CRT_COEFFICIENT2C = (MODULUS_Q1_1C * INV2C) % MODULUS_Q1C;

constexpr static auto CRT_ARRAY1 = make_array<uint64_t, CRT_COEFFICIENT1, DEGREE1>();
constexpr static auto CRT_ARRAY2 = make_array<uint64_t, CRT_COEFFICIENT2, DEGREE1>();
constexpr static auto CRT_ARRAY1C = make_array<uint64_t, CRT_COEFFICIENT1C, DEGREE1>();
constexpr static auto CRT_ARRAY2C = make_array<uint64_t, CRT_COEFFICIENT2C, DEGREE1>();

constexpr static const uint64_t AVX512_LOOP_SIZE_Deg1 = DEGREE1 >> 3;
constexpr static const uint64_t AVX512_LOOP_SIZE_Deg2 = DEGREE2 >> 3;

constexpr static const double sizePolyDeg1Q1 = static_cast<double>(DEGREE1 * CEIL_LOG_Q1) / 8192.0;
constexpr static const double sizePolyDeg1Q2 = static_cast<double>(DEGREE1 * CEIL_LOG_Q2) / 8192.0;
constexpr static const double sizePolyDeg2Q2 = static_cast<double>(DEGREE2 * CEIL_LOG_Q2) / 8192.0;
constexpr static const double sizePolyDeg2Q3 = static_cast<double>(DEGREE2 * LOG_MODULUS_Q3) / 8192.0;
constexpr static const double sizePolyDeg2Q4 = static_cast<double>(DEGREE2 * LOG_MODULUS_Q4) / 8192.0;
constexpr static const double sizePolyDeg1Q1C = static_cast<double>(DEGREE1 * CEIL_LOG_Q1C) / 8192.0;
constexpr static const double sizePolyDeg1Q2C = static_cast<double>(DEGREE1 * CEIL_LOG_Q2C) / 8192.0;
constexpr static const double sizePolyDeg2Q2C = static_cast<double>(DEGREE2 * CEIL_LOG_Q2C) / 8192.0;
constexpr static const double sizePolyDeg2Q3C = static_cast<double>(DEGREE2 * CEIL_LOG_Q3C) / 8192.0;
constexpr static const double sizePolyDeg2Q4C = static_cast<double>(DEGREE2 * LOG_MODULUS_Q4C) / 8192.0;

#pragma once

// GMP values with RAII, and the one rounding step both references share.

#include "types.h"

#include <cmath>
#include <numbers>

#include <gmp.h>

namespace exact {

class Z {
public:
    mpz_t v;

    Z() { mpz_init(v); }
    ~Z() { mpz_clear(v); }
    Z(Z const&)            = delete;
    Z& operator=(Z const&) = delete;
    Z(Z&& other) noexcept
    {
        mpz_init(v);
        mpz_swap(v, other.v);
    }
    Z& operator=(Z&& other) noexcept
    {
        mpz_swap(v, other.v);
        return *this;
    }
};

class Q {
public:
    mpq_t v;

    Q() { mpq_init(v); }
    ~Q() { mpq_clear(v); }
    Q(Q const&)            = delete;
    Q& operator=(Q const&) = delete;
    Q(Q&& other) noexcept
    {
        mpq_init(v);
        mpq_swap(v, other.v);
    }
    Q& operator=(Q&& other) noexcept
    {
        mpq_swap(v, other.v);
        return *this;
    }
};

static_assert(sizeof(unsigned long) == sizeof(u64));

// truncated to 64 bits: under one ulp of an 80-bit long double
[[nodiscard]] inline fext ToExt(mpq_srcptr q)
{
    i32 const sign { mpq_sgn(q) };
    if (sign == 0)
        return 0.L;
    Z num;
    Z den;
    Z quo;
    mpz_abs(num.v, mpq_numref(q));
    mpz_set(den.v, mpq_denref(q));
    i64 shift { 66 - (static_cast<i64>(mpz_sizeinbase(num.v, 2)) - static_cast<i64>(mpz_sizeinbase(den.v, 2))) };
    if (shift > 0)
        mpz_mul_2exp(num.v, num.v, static_cast<mp_bitcnt_t>(shift));
    else
        mpz_mul_2exp(den.v, den.v, static_cast<mp_bitcnt_t>(-shift));
    mpz_tdiv_q(quo.v, num.v, den.v);
    i64 const extra { static_cast<i64>(mpz_sizeinbase(quo.v, 2)) - 64 };
    if (extra > 0) {
        mpz_tdiv_q_2exp(quo.v, quo.v, static_cast<mp_bitcnt_t>(extra));
        shift -= extra;
    }
    fext const magnitude { std::ldexp(static_cast<fext>(mpz_get_ui(quo.v)), static_cast<i32>(-shift)) };
    return sign < 0 ? -magnitude : magnitude;
}

// from twice the axis faces' areas and twice the diagonal faces' areas
// projected onto an axis plane (sqrt(3) too small); equal rationals give equal
// areas
[[nodiscard]] inline fext Area(mpq_srcptr twiceAxis, mpq_srcptr twiceDiag)
{
    return (ToExt(twiceAxis) + std::numbers::sqrt3_v<fext> * ToExt(twiceDiag)) / 2.L;
}

}

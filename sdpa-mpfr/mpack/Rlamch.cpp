/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlamch_mpfr.cpp,v 1.7 2009/09/18 23:01:08 nakatamaho Exp $ 
 *
 * MPACK - multiple precision arithmetic library
 *
 * This file is part of MPACK.
 *
 * MPACK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License version 3
 * only, as published by the Free Software Foundation.
 *
 * MPACK is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License version 3 for more details
 * (a copy is included in the LICENSE file that accompanied this code).
 *
 * You should have received a copy of the GNU Lesser General Public License
 * version 3 along with MPACK.  If not, see
 * <http://www.gnu.org/licenses/lgpl.html>
 * for a copy of the LGPLv3 License.
 *
 ************************************************************************/

#include <mblas_mpfr.h>
#include <mlapack_mpfr.h>

//exponent bound used for "S", "M", "U", "L" and "O".
//2^(mp_bits_per_limb-8)-1 as in the GMP version, clamped so that
//2^-exp2, 2^exp2 and their reciprocals stay inside the MPFR exponent range.
static unsigned long Rlamch_exp_limit(void)
{
    unsigned long exp2 = (1UL << (mp_bits_per_limb - 8)) - 1;
    unsigned long emax = (unsigned long)mpfr_get_emax() - 2;
    unsigned long emin = (unsigned long)(-(mpfr_get_emin() + 2));
    if (exp2 > emax)
	exp2 = emax;
    if (exp2 > emin)
	exp2 = emin;
    return exp2;
}

//"E" denots we always calculate relative machine precision (e).
//where 1+e > 1, minimum of e.
mpfr_class RlamchE_mpfr(void)
{
    static mpfr_class eps;
    static int called = 0;
    if (called)
	return eps;
    mpfr_class one;
    unsigned long exp2;
    one = 1.0;
    exp2 = mpfr_get_prec(one.get_mpfr_t());
    mpfr_div_2ui(eps.get_mpfr_t(), one.get_mpfr_t(), exp2, MPFR_RNDN);
    called = 1;
    return eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchS_mpfr(void)
{
    mpfr_class sfmin;
    mpfr_class one = 1.0;
    unsigned long exp2;

    exp2 = Rlamch_exp_limit();
    mpfr_div_2ui(sfmin.get_mpfr_t(), one.get_mpfr_t(), exp2, MPFR_RNDN);
    return sfmin;


}
//"B" base  = base of the machine
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchB_mpfr(void)
{
    mpfr_class two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchP_mpfr(void)
{
    mpfr_class base, eps, prec;

    base = RlamchB_mpfr();
    eps = RlamchE_mpfr();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchN_mpfr(void)
{
    unsigned long int tmp;
    mpfr_class mtmp;
    mpfr_class mtmp2;
    tmp = mpfr_get_prec(mtmp.get_mpfr_t());
    mtmp2 = tmp;
    return mtmp2;


}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchR_mpfr(void)
{
//always rounding in addition on MPFR.
    mpfr_class mtmp;

    mtmp = 1.0;
    return mtmp;
}

//"M"
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchM_mpfr(void)
{
    unsigned long exp2;
    mpfr_class tmp;
    mpfr_class uflowmin, one=1.0; 
    exp2 = Rlamch_exp_limit();
    tmp = exp2; 
    return -tmp;

}

//"U"
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchU_mpfr(void)
{
    mpfr_class underflowmin;
    mpfr_class one = 1.0;
    unsigned long exp2;
    exp2 = Rlamch_exp_limit();
    mpfr_div_2ui(underflowmin.get_mpfr_t(), one.get_mpfr_t(), exp2, MPFR_RNDN);
    return underflowmin;
}

//"L"
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchL_mpfr(void)
{
    mpfr_class maxexp;
    unsigned long exp2;
    exp2 = Rlamch_exp_limit();
    maxexp = exp2;
    return maxexp;
}

//"O"
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchO_mpfr(void)
{
    mpfr_class overflowmax;
    mpfr_class one = 1.0;
    unsigned long exp2;

    exp2 = Rlamch_exp_limit();
    mpfr_mul_2ui(overflowmax.get_mpfr_t(), one.get_mpfr_t(), exp2, MPFR_RNDN);

    return overflowmax;
}

//"Z" :dummy
//cf.http://www.netlib.org/blas/dlamch.f
mpfr_class RlamchZ_mpfr(void)
{
    mpfr_class mtemp = 0.0;
    return mtemp;
}

mpfr_class Rlamch_mpfr(const char *cmach)
{
    if (Mlsame_mpfr(cmach, "E"))
	return RlamchE_mpfr();
    if (Mlsame_mpfr(cmach, "S"))
	return RlamchS_mpfr();
    if (Mlsame_mpfr(cmach, "B"))
	return RlamchB_mpfr();
    if (Mlsame_mpfr(cmach, "P"))
	return RlamchP_mpfr();
    if (Mlsame_mpfr(cmach, "N"))
	return RlamchN_mpfr();
    if (Mlsame_mpfr(cmach, "R"))
	return RlamchR_mpfr();
    if (Mlsame_mpfr(cmach, "M"))
	return RlamchM_mpfr();
    if (Mlsame_mpfr(cmach, "U"))
	return RlamchU_mpfr();
    if (Mlsame_mpfr(cmach, "L"))
	return RlamchL_mpfr();
    if (Mlsame_mpfr(cmach, "O"))
	return RlamchO_mpfr();

    Mxerbla_mpfr("Rlamch", 1);
    return RlamchZ_mpfr();
}

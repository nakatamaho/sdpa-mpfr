/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2009 by Nakata, Maho
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

#ifndef _MUTILS_MPFR_H_
#define _MUTILS_MPFR_H_

using std::max;
using std::min;

mpfr_class Msign(mpfr_class a, mpfr_class b);
double cast2double(mpfr_class a);
int M2int(mpfr_class a);

//implementation of sign transfer function.
//returns |a| with the sign of b; |a| when b is zero.
inline mpfr_class
Msign(mpfr_class a, mpfr_class b)
{
    mpfr_class mtmp;
    mpfr_abs(mtmp.get_mpfr_t(), a.get_mpfr_t(), MPFR_RNDN);
    if (mpfr_sgn(b.get_mpfr_t()) < 0) {
	mpfr_neg(mtmp.get_mpfr_t(), mtmp.get_mpfr_t(), MPFR_RNDN);
    }
    return mtmp;
}

inline double
cast2double(mpfr_class a)
{
    return a.get_d();
}

inline int
M2int(mpfr_class a)
{
    a = a + 0.5;
    return (int)mpfr_get_si(a.get_mpfr_t(), MPFR_RNDD);
}

#endif

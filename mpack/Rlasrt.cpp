/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlasrt.cpp,v 1.7 2009/09/25 04:00:39 nakatamaho Exp $ 
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
#include <stdlib.h>

int compare_mpf_gt(const mpfr_class* a, const mpfr_class* b)
{
    if (*a > *b)
        return 1;
    if (*a == *b)
        return 0;
    if (*a < *b)
        return -1;
    return 0; //never occurs
}

int compare_mpf_lt(const mpfr_class* a, const mpfr_class* b)
{
    if (*a > *b)
        return -1;
    if (*a == *b)
        return 0;
    if (*a < *b)
        return 1;
    return 0; //never occurs
}

void Rlasrt(const char* id, mpackint n, mpfr_class* d, mpackint* info)
{
    //Error check
    if (!Mlsame_mpfr(id, "I") && !Mlsame_mpfr(id, "D")) {
        *info = -1;
        Mxerbla_mpfr("Rlasrt", -(*info));
        return;
    }
    if (n < 0) {
        *info = -2;
        Mxerbla_mpfr("Rlasrt", -(*info));
        return;
    }
    if (Mlsame_mpfr(id, "I")) {
        qsort(d, n, sizeof(mpfr_class), (int (*)(const void*, const void*))compare_mpf_gt);
    }
    if (Mlsame_mpfr(id, "d")) {
        qsort(d, n, sizeof(mpfr_class), (int (*)(const void*, const void*))compare_mpf_lt);
    }
    *info = 0;
}

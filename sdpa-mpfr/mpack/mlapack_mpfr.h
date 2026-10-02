/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: mlapack_mpfr.h,v 1.6 2009/09/22 20:27:18 nakatamaho Exp $ 
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

#ifndef _MLAPACK_MPFR_H_
#define _MLAPACK_MPFR_H_

/* this is a subset of mpack for SDPA-MPFR only */
/* http://mplapack.sourceforge.net/ */

/* mlapack prototypes */
void Rsteqr(const char* compz, mpackint n, mpfr_class* d, mpfr_class* e, mpfr_class* Z,
            mpackint ldz, mpfr_class* work, mpackint* info);
void Rsyev(const char* jobz, const char* uplo, mpackint n, mpfr_class* A, mpackint lda,
           mpfr_class* w, mpfr_class* work, mpackint* lwork, mpackint* info);
void Rpotrf(const char* uplo, mpackint n, mpfr_class* A, mpackint lda, mpackint* info);
mpackint iMlaenv_mpfr(mpackint ispec, const char* name, const char* opts, mpackint n1, mpackint n2,
                      mpackint n3, mpackint n4);
mpfr_class Rlamch_mpfr(const char* cmach);
mpfr_class Rlansy(const char* norm, const char* uplo, mpackint n, mpfr_class* A, mpackint lda,
                  mpfr_class* work);
void Rlascl(const char* type, mpackint kl, mpackint ku, mpfr_class cfrom, mpfr_class cto,
            mpackint m, mpackint n, mpfr_class* A, mpackint lda, mpackint* info);
void Rsytrd(const char* uplo, mpackint n, mpfr_class* A, mpackint lda, mpfr_class* d, mpfr_class* e,
            mpfr_class* tau, mpfr_class* work, mpackint lwork, mpackint* info);
void Rsytd2(const char* uplo, mpackint n, mpfr_class* A, mpackint lda, mpfr_class* d, mpfr_class* e,
            mpfr_class* tau, mpackint* info);
mpfr_class Rlanst(const char* norm, mpackint n, mpfr_class* d, mpfr_class* e);
void Rlae2(mpfr_class a, mpfr_class b, mpfr_class c, mpfr_class* rt1, mpfr_class* rt2);
mpfr_class Rlapy2(mpfr_class x, mpfr_class y);
void Rlasrt(const char* id, mpackint n, mpfr_class* d, mpackint* info);
void Rorgql(mpackint m, mpackint n, mpackint k, mpfr_class* A, mpackint lda, mpfr_class* tau,
            mpfr_class* work, mpackint lwork, mpackint* info);
void Rorgqr(mpackint m, mpackint n, mpackint k, mpfr_class* A, mpackint lda, mpfr_class* tau,
            mpfr_class* work, mpackint lwork, mpackint* info);
void Rlarfg(mpackint N, mpfr_class* alpha, mpfr_class* x, mpackint incx, mpfr_class* tau);
void Rlassq(mpackint n, mpfr_class* x, mpackint incx, mpfr_class* scale, mpfr_class* sumsq);
void Rorg2l(mpackint m, mpackint n, mpackint k, mpfr_class* A, mpackint lda, mpfr_class* tau,
            mpfr_class* work, mpackint* info);
void Rlarft(const char* direct, const char* storev, mpackint n, mpackint k, mpfr_class* v,
            mpackint ldv, mpfr_class* tau, mpfr_class* t, mpackint ldt);
void Rlarfb(const char* side, const char* trans, const char* direct, const char* storev, mpackint m,
            mpackint n, mpackint k, mpfr_class* V, mpackint ldv, mpfr_class* T, mpackint ldt,
            mpfr_class* C, mpackint ldc, mpfr_class* work, mpackint ldwork);
void Rorg2r(mpackint m, mpackint n, mpackint k, mpfr_class* A, mpackint lda, mpfr_class* tau,
            mpfr_class* work, mpackint* info);
void Rlarf(const char* side, mpackint m, mpackint n, mpfr_class* v, mpackint incv, mpfr_class tau,
           mpfr_class* C, mpackint ldc, mpfr_class* work);
void Rpotf2(const char* uplo, mpackint n, mpfr_class* A, mpackint lda, mpackint* info);
void Rlaset(const char* uplo, mpackint m, mpackint n, mpfr_class alpha, mpfr_class beta,
            mpfr_class* A, mpackint lda);
void Rlaev2(mpfr_class a, mpfr_class b, mpfr_class c, mpfr_class* rt1, mpfr_class* rt2,
            mpfr_class* cs1, mpfr_class* sn1);
void Rlasr(const char* side, const char* pivot, const char* direct, mpackint m, mpackint n,
           mpfr_class* c, mpfr_class* s, mpfr_class* A, mpackint lda);
void Rlartg(mpfr_class f, mpfr_class g, mpfr_class* cs, mpfr_class* sn, mpfr_class* r);
void Rlatrd(const char* uplo, mpackint n, mpackint nb, mpfr_class* A, mpackint lda, mpfr_class* e,
            mpfr_class* tau, mpfr_class* w, mpackint ldw);
void Rsterf(mpackint n, mpfr_class* d, mpfr_class* e, mpackint* info);
void Rorgtr(const char* uplo, mpackint n, mpfr_class* a, mpackint lda, mpfr_class* tau,
            mpfr_class* work, mpackint lwork, mpackint* info);
#endif

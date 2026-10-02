/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */

#include <sdpa_linear.h>
#include <sdpa_dataset.h>

namespace sdpa {

mpfr_class Lal::getMinEigen(DenseMatrix& lMat, DenseMatrix& xMat, DenseMatrix& Q, Vector& out,
                            Vector& b, Vector& r, Vector& q, Vector& qold, Vector& w, Vector& tmp,
                            Vector& diagVec, Vector& diagVec2, Vector& workVec)
{
    mpfr_class alpha, beta, value;
    mpfr_class min = 1.0e+51, min_old = 1.0e+52, min_min = 1.0e+50;
    mpfr_class error = 1.0e+10;

    int nDim = xMat.nRow;
    int k = 0, kk = 0;

    diagVec.initialize(min_min);
    diagVec2.setZero();
    q.setZero();
    r.initialize(MONE);
    beta = sqrt((mpfr_class)nDim); // norm of "r"

    // nakata 2004/12/12
    while (k < nDim && k < sqrt((mpfr_class)nDim) + 10 && beta > 1.0e-16 &&
           (abs(min - min_old) > (1.0e-5) * abs(min) + (1.0e-8) ||
            abs(error * beta) > (1.0e-2) * abs(min) + (1.0e-4))) {
        qold.copyFrom(q);
        value = MONE / beta;
        Lal::let(q, '=', r, '*', &value);

        // w = (lMat^T)*q
        w.copyFrom(q);
        Rtrmv("Lower", "Transpose", "NotUnit", nDim, lMat.de_ele, nDim, w.ele, 1);

        Lal::let(tmp, '=', xMat, '*', w);
        w.copyFrom(tmp);
        Rtrmv("Lower", "NoTranspose", "NotUnit", nDim, lMat.de_ele, nDim, w.ele, 1);
        // w = lMat*xMat*(lMat^T)*q
        Lal::let(alpha, '=', q, '.', w);
        diagVec.ele[k] = alpha;
        Lal::let(r, '=', w, '-', q, &alpha);
        Lal::let(r, '=', r, '-', qold, &beta);

        if (kk >= sqrt((mpfr_class)k) || k == nDim - 1 || k > sqrt((mpfr_class)nDim + 9)) {
            kk = 0;
            out.copyFrom(diagVec);
            b.copyFrom(diagVec2);
            out.ele[nDim - 1] = diagVec.ele[k];
            b.ele[nDim - 1] = 0.0;

            mpackint info;
            int kp1 = k + 1;
            Rsteqr("I_withEigenvalues", kp1, out.ele, b.ele, Q.de_ele, Q.nRow, workVec.ele, &info);

            if (info < 0) {
                rError(" rLanczos :: bad argument " << -info << " Q.nRow = " << Q.nRow
                                                    << ": nDim = " << nDim << ": kp1 = " << kp1);
            } else if (info > 0) {
                rMessage(" rLanczos :: cannot converge " << info);
                break;
            }

            min_old = min;
            // out have eigen values with ascending order.
            min = out.ele[0];
            error = Q.de_ele[k];

        } // end of 'if ( kk>=sqrt(k) ...)'

        Lal::let(value, '=', r, '.', r);
        beta = sqrt(value);
        diagVec2.ele[k] = beta;
        ++k;
        ++kk;
    } // end of while
    return min - abs(error * beta);
}

mpfr_class Lal::getMinEigenValue(DenseMatrix& aMat, Vector& eigenVec, Vector& workVec)
{
    // aMat is rewritten.
    // aMat must be symmetric.
    // eigenVec is the space of eigen values
    // and needs memory of length aMat.nRow
    // workVec is temporary space and needs
    // 3*aMat.nRow-1 length memory.
    mpackint N = aMat.nRow;
    mpackint LWORK, info;
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        LWORK = 3 * N - 1;
        // "N" means that we need not eigen vectors
        // "L" means that we refer only lower triangular.
        Rsyev("NonVectors", "Lower", N, aMat.de_ele, N, eigenVec.ele, workVec.ele, &LWORK, &info);
        if (info != 0) {
            if (info < 0) {
                rMessage("getMinEigenValue:: info is mistaken " << info);
            } else {
                rMessage("getMinEigenValue:: cannot decomposition");
            }
            exit(0);
            return 0.0;
        }
        return eigenVec.ele[0];
        // Eigen values are sorted by ascending order.
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return 0.0;
}

bool Lal::getInnerProduct(mpfr_class& ret, Vector& aVec, Vector& bVec)
{
    int N = aVec.nDim;
    if (N != bVec.nDim) {
        rError("getInnerProduct:: different memory size");
    }
    ret = Rdot(N, aVec.ele, 1, bVec.ele, 1);

    return _SUCCESS;
}

bool Lal::getInnerProduct(mpfr_class& ret, DenseMatrix& aMat, DenseMatrix& bMat)
{
    if (aMat.nRow != bMat.nRow || aMat.nCol != bMat.nCol) {
        rError("getInnerProduct:: different memory size");
    }
    int length;
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        length = aMat.nRow * aMat.nCol;
        ret = Rdot(length, aMat.de_ele, 1, bMat.de_ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::getInnerProduct(mpfr_class& ret, SparseMatrix& aMat, DenseMatrix& bMat)
{
    if (aMat.nRow != bMat.nRow || aMat.nCol != bMat.nCol) {
        rError("getInnerProduct:: different memory size");
    }
    int length;
    int amari, shou;

    switch (aMat.type) {
    case SparseMatrix::SPARSE:
        // Attension: in SPARSE case, only half elements
        // are stored. And bMat must be DENSE case.
        ret = 0.0;
        amari = aMat.NonZeroCount % 4;
        shou = aMat.NonZeroCount / 4;
        for (int index = 0; index < amari; ++index) {
            int i = aMat.row_index[index];
            int j = aMat.column_index[index];
            mpfr_class value = aMat.sp_ele[index];
            if (i == j) {
                ret += value * bMat.de_ele[i + bMat.nRow * j];
            } else {
                ret += value * (bMat.de_ele[i + bMat.nRow * j] + bMat.de_ele[j + bMat.nRow * i]);
            }
        }
        for (int index = amari, counter = 0; counter < shou; ++counter, index += 4) {
            int i1 = aMat.row_index[index];
            int j1 = aMat.column_index[index];
            mpfr_class value1 = aMat.sp_ele[index];
            mpfr_class ret1;
            if (i1 == j1) {
                ret1 = value1 * bMat.de_ele[i1 + bMat.nRow * j1];
            } else {
                ret1 =
                    value1 * (bMat.de_ele[i1 + bMat.nRow * j1] + bMat.de_ele[j1 + bMat.nRow * i1]);
            }
            int i2 = aMat.row_index[index + 1];
            int j2 = aMat.column_index[index + 1];
            mpfr_class value2 = aMat.sp_ele[index + 1];
            mpfr_class ret2;
            if (i2 == j2) {
                ret2 = value2 * bMat.de_ele[i2 + bMat.nRow * j2];
            } else {
                ret2 =
                    value2 * (bMat.de_ele[i2 + bMat.nRow * j2] + bMat.de_ele[j2 + bMat.nRow * i2]);
            }
            int i3 = aMat.row_index[index + 2];
            int j3 = aMat.column_index[index + 2];
            mpfr_class value3 = aMat.sp_ele[index + 2];
            mpfr_class ret3;
            if (i3 == j3) {
                ret3 = value3 * bMat.de_ele[i3 + bMat.nRow * j3];
            } else {
                ret3 =
                    value3 * (bMat.de_ele[i3 + bMat.nRow * j3] + bMat.de_ele[j3 + bMat.nRow * i3]);
            }
            int i4 = aMat.row_index[index + 3];
            int j4 = aMat.column_index[index + 3];
            mpfr_class value4 = aMat.sp_ele[index + 3];
            mpfr_class ret4;
            if (i4 == j4) {
                ret4 = value4 * bMat.de_ele[i4 + bMat.nRow * j4];
            } else {
                ret4 =
                    value4 * (bMat.de_ele[i4 + bMat.nRow * j4] + bMat.de_ele[j4 + bMat.nRow * i4]);
            }
            ret += (ret1 + ret2 + ret3 + ret4);
        }
        break;
    case SparseMatrix::DENSE:
        length = aMat.nRow * aMat.nCol;
        ret = Rdot(length, aMat.de_ele, 1, bMat.de_ele, 1);
        break;
    }
    return _SUCCESS;
}

bool Lal::getCholesky(DenseMatrix& retMat, DenseMatrix& aMat)
{
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.type != aMat.type) {
        rError("getCholesky:: different memory size");
    }
    int length, shou, amari;
    mpackint info;
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        length = retMat.nRow * retMat.nCol;
        Rcopy(length, aMat.de_ele, 1, retMat.de_ele, 1);
        Rpotrf("Lower", retMat.nRow, retMat.de_ele, retMat.nRow, &info);
        if (info != 0) {
            rMessage("cannot cholesky decomposition");
            rMessage("Could you try with smaller gammaStar?");
            return FAILURE;
        }
        // Make matrix as lower triangular matrix
        for (int j = 0; j < retMat.nCol; ++j) {
            shou = j / 4;
            amari = j % 4;
            for (int i = 0; i < amari; ++i) {
                retMat.de_ele[i + retMat.nCol * j] = 0.0;
            }
            for (int i = amari, count = 0; count < shou; ++count, i += 4) {
                retMat.de_ele[i + retMat.nCol * j] = 0.0;
                retMat.de_ele[i + 1 + retMat.nCol * j] = 0.0;
                retMat.de_ele[i + 2 + retMat.nCol * j] = 0.0;
                retMat.de_ele[i + 3 + retMat.nCol * j] = 0.0;
            }
        }
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

// nakata 2004/12/01
// modified 2008/05/20    "aMat.sp_ele[indexA1] = 0.0;"
// aMat = L L^T
bool Lal::getCholesky(SparseMatrix& aMat, int* diagonalIndex)
{
    int nDIM = aMat.nRow;
    int indexA1, indexA2, indexB2;
    int i, k1, k2, k3;
    mpfr_class tmp, tmp2;
    int tmp3;

    if (aMat.type != SparseMatrix::SPARSE) {
        rError("Lal::getCholesky aMat is not sparse format");
    }

    for (i = 0; i < nDIM; ++i) {
        indexA1 = diagonalIndex[i];
        indexA2 = diagonalIndex[i + 1];
        if (aMat.sp_ele[indexA1] < 0.0) {
            aMat.sp_ele[indexA1] = 0.0;
        } else {
            // inverse diagonal
            aMat.sp_ele[indexA1] = 1.0 / sqrt(aMat.sp_ele[indexA1]);
        }
        for (k1 = indexA1 + 1; k1 < indexA2; ++k1) {
            aMat.sp_ele[k1] *= aMat.sp_ele[indexA1];
        }
        for (k1 = indexA1 + 1; k1 < indexA2; ++k1) {
            tmp = aMat.sp_ele[k1];
            k3 = diagonalIndex[aMat.column_index[k1]];
            indexB2 = diagonalIndex[aMat.column_index[k1] + 1];
            for (k2 = k1; k2 < indexA2; ++k2) {
                tmp2 = aMat.sp_ele[k2];
                tmp3 = aMat.column_index[k2];
                for (; k3 < indexB2; ++k3) {
                    if (aMat.column_index[k3] == tmp3) {
                        aMat.sp_ele[k3] -= tmp * tmp2;
                        k3++;
                        break;
                    }
                }
            }
        }
    }
    return true;
}

bool Lal::getInvLowTriangularMatrix(DenseMatrix& retMat, DenseMatrix& aMat)
{
    // Make inverse with refference only to lower triangular.
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.type != aMat.type) {
        rError("getCholesky:: different memory size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        retMat.setIdentity();
        Rtrsm("Left", "Lower", "NoTraspose", "NonUnitDiagonal", aMat.nRow, aMat.nCol, MONE,
              aMat.de_ele, aMat.nRow, retMat.de_ele, retMat.nRow);
        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::getSymmetrize(DenseMatrix& aMat)
{
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        if (aMat.nRow != aMat.nCol) {
            rError("getSymmetrize:: different memory size");
        }
        for (int index = 0; index < aMat.nRow - 1; ++index) {
            int index1 = index + index * aMat.nRow + 1;
            int index2 = index + (index + 1) * aMat.nRow;
            int length = aMat.nRow - 1 - index;
            // aMat.de_ele[index1] += aMat.de_ele[index2]
            Raxpy(length, MONE, &aMat.de_ele[index2], aMat.nRow, &aMat.de_ele[index1], 1);
            // aMat.de_ele[index1] /= 2.0
            mpfr_class half = 0.5;
            Rscal(length, half, &aMat.de_ele[index1], 1);
            // aMat.de_ele[index2] = aMat.de_ele[index1]
            Rcopy(length, &aMat.de_ele[index1], 1, &aMat.de_ele[index2], aMat.nRow);
        }
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::choleskyFactorWithAdjust(DenseMatrix& aMat)
{
    mpackint info = 0;
    Rpotrf("Lower", aMat.nRow, aMat.de_ele, aMat.nRow, &info);
    if (info < 0) {
        rMessage("cholesky argument is wrong " << -info);
    } else if (info > 0) {
        rMessage("cholesky miss condition :: not positive definite" << " :: info = " << info);

        return FAILURE;
    }
    return _SUCCESS;
}

bool Lal::solveSystems(Vector& xVec, DenseMatrix& aMat, Vector& bVec)
{
    // aMat must have done Cholesky factorized.
    if (aMat.nCol != xVec.nDim || aMat.nRow != bVec.nDim || aMat.nRow != aMat.nCol) {
        rError("solveSystems:: different memory size");
    }
    if (aMat.type != DenseMatrix::DENSE) {
        rError("solveSystems:: matrix type must be DENSE");
    }
    xVec.copyFrom(bVec);
    Rtrsv("Lower", "NoTranspose", "NonUnit", aMat.nRow, aMat.de_ele, aMat.nCol, xVec.ele, 1);
    Rtrsv("Lower", "Transpose", "NonUnit", aMat.nRow, aMat.de_ele, aMat.nCol, xVec.ele, 1);
    return _SUCCESS;
}

// nakata 2004/12/01
bool Lal::solveSystems(Vector& xVec, SparseMatrix& aMat, Vector& bVec)
{
    // Attension: in SPARSE case, only half elements
    // are stored. And bMat must be DENSE case.
    xVec.copyFrom(bVec);
    for (int index = 0; index < aMat.NonZeroCount; ++index) {
        int i = aMat.row_index[index];
        int j = aMat.column_index[index];
        mpfr_class value = aMat.sp_ele[index];
        if (i == j) {
            xVec.ele[i] *= value;
        } else {
            xVec.ele[j] -= value * xVec.ele[i];
        }
    }
    for (int index = aMat.NonZeroCount - 1; index >= 0; --index) {
        int i = aMat.row_index[index];
        int j = aMat.column_index[index];
        mpfr_class value = aMat.sp_ele[index];
        value = aMat.sp_ele[index];
        if (i == j) {
            xVec.ele[i] *= value;
        } else {
            xVec.ele[i] -= value * xVec.ele[j];
        }
    }
    return _SUCCESS;
}

bool Lal::multiply(DenseMatrix& retMat, DenseMatrix& aMat, DenseMatrix& bMat, mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nRow || aMat.nCol != bMat.nRow || bMat.nCol != retMat.nCol ||
        retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("multiply :: different matrix size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
            // attension::scalar is loval variable.
        }
        Rgemm("NoTranspose", "NoTranspose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar,
              aMat.de_ele, aMat.nRow, bMat.de_ele, bMat.nRow, 0.0, retMat.de_ele, retMat.nRow);

        break;
    case DenseMatrix::COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::multiply(DenseMatrix& retMat, SparseMatrix& aMat, DenseMatrix& bMat, mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nRow || aMat.nCol != bMat.nRow || bMat.nCol != retMat.nCol) {
        rError("multiply :: different matrix size");
    }
    retMat.setZero();
    switch (aMat.type) {
    case SparseMatrix::SPARSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("multiply :: different matrix type");
        }
        if (scalar == NULL) {
            for (int index = 0; index < aMat.NonZeroCount; ++index) {
                int i = aMat.row_index[index];
                int j = aMat.column_index[index];
                mpfr_class value = aMat.sp_ele[index];
                if (i != j) {
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + i,
                          retMat.nRow);
                    Raxpy(bMat.nCol, value, bMat.de_ele + i, bMat.nRow, retMat.de_ele + j,
                          retMat.nRow);
                } else {
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + j,
                          retMat.nRow);
                }
            } // end of 'for index'
        } else { // scalar!=NULL
            for (int index = 0; index < aMat.NonZeroCount; ++index) {
                int i = aMat.row_index[index];
                int j = aMat.column_index[index];
                mpfr_class value = aMat.sp_ele[index] * (*scalar);
                if (i != j) {
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + i,
                          retMat.nRow);
                    Raxpy(bMat.nCol, value, bMat.de_ele + i, bMat.nRow, retMat.de_ele + j,
                          retMat.nRow);
                } else {
                    Raxpy(bMat.nCol, value, bMat.de_ele + j, bMat.nRow, retMat.de_ele + j,
                          retMat.nRow);
                }
            } // end of 'for index'
        } // end of 'if (scalar==NULL)
        break;
    case SparseMatrix::DENSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("multiply :: different matrix type");
        }
        if (scalar == NULL) {
            scalar = &MONE;
            // attension:: scalar is local variable.
        }
        Rgemm("NoTranspose", "NoTranspose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar,
              aMat.de_ele, aMat.nRow, bMat.de_ele, bMat.nRow, 0.0, retMat.de_ele, retMat.nRow);
        break;

    } // end of switch

    return _SUCCESS;
}

bool Lal::multiply(DenseMatrix& retMat, DenseMatrix& aMat, mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.type != aMat.type) {
        rError("multiply :: different matrix size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    int length;
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        length = retMat.nRow * retMat.nCol;
        Rcopy(length, aMat.de_ele, 1, retMat.de_ele, 1);
        Rscal(length, *scalar, retMat.de_ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::multiply(Vector& retVec, Vector& aVec, mpfr_class* scalar)
{
    if (retVec.nDim != aVec.nDim) {
        rError("multiply :: different vector size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    Rcopy(retVec.nDim, aVec.ele, 1, retVec.ele, 1);
    Rscal(retVec.nDim, *scalar, retVec.ele, 1);
    return _SUCCESS;
}

bool Lal::multiply(Vector& retVec, DenseMatrix& aMat, Vector& bVec, mpfr_class* scalar)
{
    if (retVec.nDim != aMat.nRow || aMat.nCol != bVec.nDim || bVec.nDim != retVec.nDim) {
        rError("multiply :: different matrix size");
    }
    switch (aMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
        }
        Rgemv("NoTranspose", aMat.nRow, aMat.nCol, *scalar, aMat.de_ele, aMat.nRow, bVec.ele, 1,
              0.0, retVec.ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::tran_multiply(DenseMatrix& retMat, DenseMatrix& aMat, DenseMatrix& bMat,
                        mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nCol || aMat.nRow != bMat.nRow || bMat.nCol != retMat.nCol ||
        retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("multiply :: different matrix size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
            // scalar is local variable
        }
        // The Point is the first argument is "Transpose".
        Rgemm("Transpose", "NoTranspose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar, aMat.de_ele,
              aMat.nCol, bMat.de_ele, bMat.nRow, 0.0, retMat.de_ele, retMat.nRow);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }

    return _SUCCESS;
}

bool Lal::multiply_tran(DenseMatrix& retMat, DenseMatrix& aMat, DenseMatrix& bMat,
                        mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nRow || aMat.nCol != bMat.nCol || bMat.nRow != retMat.nRow ||
        retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("multiply :: different matrix size");
    }
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        if (scalar == NULL) {
            scalar = &MONE;
        }
        // The Point is the first argument is "NoTranspose".
        Rgemm("NoTranspose", "Transpose", retMat.nRow, retMat.nCol, aMat.nCol, *scalar, aMat.de_ele,
              aMat.nRow, bMat.de_ele, bMat.nCol, 0.0, retMat.de_ele, retMat.nRow);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::plus(Vector& retVec, Vector& aVec, Vector& bVec, mpfr_class* scalar)
{
    if (retVec.nDim != aVec.nDim || aVec.nDim != bVec.nDim) {
        rError("plus :: different matrix size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    if (retVec.ele != aVec.ele) {
        Rcopy(retVec.nDim, aVec.ele, 1, retVec.ele, 1);
    }
    Raxpy(retVec.nDim, *scalar, bVec.ele, 1, retVec.ele, 1);
    return _SUCCESS;
}

bool Lal::plus(DenseMatrix& retMat, DenseMatrix& aMat, DenseMatrix& bMat, mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.nRow != bMat.nRow ||
        retMat.nCol != bMat.nCol || retMat.type != aMat.type || retMat.type != bMat.type) {
        rError("plus :: different matrix size");
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    int length;
    switch (retMat.type) {
    case DenseMatrix::DENSE:
        length = retMat.nRow * retMat.nCol;
        if (retMat.de_ele != aMat.de_ele) {
            Rcopy(length, aMat.de_ele, 1, retMat.de_ele, 1);
        }
        Raxpy(length, *scalar, bMat.de_ele, 1, retMat.de_ele, 1);
        break;
    case DenseMatrix::COMPLETION:
        rError("no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

bool Lal::plus(DenseMatrix& retMat, SparseMatrix& aMat, DenseMatrix& bMat, mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.nRow != bMat.nRow ||
        retMat.nCol != bMat.nCol) {
        rError("plus :: different matrix size");
    }
    // ret = (*scalar) * b
    if (multiply(retMat, bMat, scalar) == FAILURE) {
        return FAILURE;
    }
    int length;
    int shou, amari;
    switch (aMat.type) {
    case SparseMatrix::SPARSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
        shou = aMat.NonZeroCount / 4;
        amari = aMat.NonZeroCount % 4;
        for (int index = 0; index < amari; ++index) {
            int i = aMat.row_index[index];
            int j = aMat.column_index[index];
            mpfr_class value = aMat.sp_ele[index];
            if (i != j) {
                retMat.de_ele[i + retMat.nCol * j] += value;
                retMat.de_ele[j + retMat.nCol * i] += value;
            } else {
                retMat.de_ele[i + retMat.nCol * i] += value;
            }
        } // end of 'for index'
        for (int index = amari, counter = 0; counter < shou; ++counter, index += 4) {
            int i1 = aMat.row_index[index];
            int j1 = aMat.column_index[index];
            mpfr_class value1 = aMat.sp_ele[index];
            if (i1 != j1) {
                retMat.de_ele[i1 + retMat.nCol * j1] += value1;
                retMat.de_ele[j1 + retMat.nCol * i1] += value1;
            } else {
                retMat.de_ele[i1 + retMat.nCol * i1] += value1;
            }
            int i2 = aMat.row_index[index + 1];
            int j2 = aMat.column_index[index + 1];
            mpfr_class value2 = aMat.sp_ele[index + 1];
            if (i2 != j2) {
                retMat.de_ele[i2 + retMat.nCol * j2] += value2;
                retMat.de_ele[j2 + retMat.nCol * i2] += value2;
            } else {
                retMat.de_ele[i2 + retMat.nCol * i2] += value2;
            }
            int i3 = aMat.row_index[index + 2];
            int j3 = aMat.column_index[index + 2];
            mpfr_class value3 = aMat.sp_ele[index + 2];
            if (i3 != j3) {
                retMat.de_ele[i3 + retMat.nCol * j3] += value3;
                retMat.de_ele[j3 + retMat.nCol * i3] += value3;
            } else {
                retMat.de_ele[i3 + retMat.nCol * i3] += value3;
            }
            int i4 = aMat.row_index[index + 3];
            int j4 = aMat.column_index[index + 3];
            mpfr_class value4 = aMat.sp_ele[index + 3];
            if (i4 != j4) {
                retMat.de_ele[i4 + retMat.nCol * j4] += value4;
                retMat.de_ele[j4 + retMat.nCol * i4] += value4;
            } else {
                retMat.de_ele[i4 + retMat.nCol * i4] += value4;
            }
        } // end of 'for index'
        break;
    case SparseMatrix::DENSE:
        if (retMat.type != DenseMatrix::DENSE || bMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
        length = retMat.nRow * retMat.nCol;
        Raxpy(length, 1.0, aMat.de_ele, 1, retMat.de_ele, 1);
        break;
    } // end of switch
    return _SUCCESS;
}

bool Lal::plus(DenseMatrix& retMat, DenseMatrix& aMat, SparseMatrix& bMat, mpfr_class* scalar)
{
    if (retMat.nRow != aMat.nRow || retMat.nCol != aMat.nCol || retMat.nRow != bMat.nRow ||
        retMat.nCol != bMat.nCol) {
        rError("plus :: different matrix size");
    }
    // ret = a
    if (retMat.copyFrom(aMat) == FAILURE) {
        return FAILURE;
    }
    if (scalar == NULL) {
        scalar = &MONE;
    }
    int length, shou, amari;
    switch (bMat.type) {
    case SparseMatrix::SPARSE:
        if (retMat.type != DenseMatrix::DENSE || aMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
        shou = bMat.NonZeroCount / 4;
        amari = bMat.NonZeroCount % 4;
        for (int index = 0; index < amari; ++index) {
            int i = bMat.row_index[index];
            int j = bMat.column_index[index];
            mpfr_class value = bMat.sp_ele[index] * (*scalar);
            if (i != j) {
                retMat.de_ele[i + retMat.nCol * j] += value;
                retMat.de_ele[j + retMat.nCol * i] += value;
            } else {
                retMat.de_ele[i + retMat.nCol * i] += value;
            }
        } // end of 'for index'
        for (int index = amari, counter = 0; counter < shou; ++counter, index += 4) {
            int i1 = bMat.row_index[index];
            int j1 = bMat.column_index[index];
            mpfr_class value1 = bMat.sp_ele[index] * (*scalar);
            if (i1 != j1) {
                retMat.de_ele[i1 + retMat.nCol * j1] += value1;
                retMat.de_ele[j1 + retMat.nCol * i1] += value1;
            } else {
                retMat.de_ele[i1 + retMat.nCol * i1] += value1;
            }
            int i2 = bMat.row_index[index + 1];
            int j2 = bMat.column_index[index + 1];
            mpfr_class value2 = bMat.sp_ele[index + 1] * (*scalar);
            if (i2 != j2) {
                retMat.de_ele[i2 + retMat.nCol * j2] += value2;
                retMat.de_ele[j2 + retMat.nCol * i2] += value2;
            } else {
                retMat.de_ele[i2 + retMat.nCol * i2] += value2;
            }
            int i3 = bMat.row_index[index + 2];
            int j3 = bMat.column_index[index + 2];
            mpfr_class value3 = bMat.sp_ele[index + 2] * (*scalar);
            if (i3 != j3) {
                retMat.de_ele[i3 + retMat.nCol * j3] += value3;
                retMat.de_ele[j3 + retMat.nCol * i3] += value3;
            } else {
                retMat.de_ele[i3 + retMat.nCol * i3] += value3;
            }
            int i4 = bMat.row_index[index + 3];
            int j4 = bMat.column_index[index + 3];
            mpfr_class value4 = bMat.sp_ele[index + 3] * (*scalar);
            if (i4 != j4) {
                retMat.de_ele[i4 + retMat.nCol * j4] += value4;
                retMat.de_ele[j4 + retMat.nCol * i4] += value4;
            } else {
                retMat.de_ele[i4 + retMat.nCol * i4] += value4;
            }
        } // end of 'for index'
        break;
    case SparseMatrix::DENSE:
        if (retMat.type != DenseMatrix::DENSE || aMat.type != DenseMatrix::DENSE) {
            rError("plus :: different matrix type");
        }
        length = retMat.nRow * retMat.nCol;
        Raxpy(length, *scalar, bMat.de_ele, 1, retMat.de_ele, 1);
        break;
    } // end of switch
    return _SUCCESS;
}

// ret = a '*' (*scalar)
bool Lal::let(Vector& retVec, const char eq, Vector& aVec, const char op, mpfr_class* scalar)
{
    switch (op) {
    case '*':
        return multiply(retVec, aVec, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool Lal::let(Vector& retVec, const char eq, Vector& aVec, const char op, Vector& bVec,
              mpfr_class* scalar)
{
    mpfr_class minus_scalar;
    switch (op) {
    case '+':
        return plus(retVec, aVec, bVec, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retVec, aVec, bVec, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' '*' 't' 'T' b*(*scalar)
bool Lal::let(DenseMatrix& retMat, const char eq, DenseMatrix& aMat, const char op,
              DenseMatrix& bMat, mpfr_class* scalar)
{
    mpfr_class minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '*':
        return multiply(retMat, aMat, bMat, scalar);
        break;
    case 't':
        // ret = aMat**T * bMat
        return tran_multiply(retMat, aMat, bMat, scalar);
        break;
    case 'T':
        // ret = aMat * bMat**T
        return multiply_tran(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' '*' b*(*scalar)
bool Lal::let(DenseMatrix& retMat, const char eq, SparseMatrix& aMat, const char op,
              DenseMatrix& bMat, mpfr_class* scalar)
{
    mpfr_class minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '*':
        return multiply(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = aMat '*' '/' bVec
bool Lal::let(Vector& rVec, const char eq, DenseMatrix& aMat, const char op, Vector& bVec)
{
    switch (op) {
    case '*':
        return multiply(rVec, aMat, bVec, NULL);
        break;
    case '/':
        // ret = aMat^{-1} * bVec;
        // aMat is positive definite
        // and already colesky factorized.
        return solveSystems(rVec, aMat, bVec);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// nakata 2004/12/01
// ret = aMat '*' '/' bVec
bool Lal::let(Vector& rVec, const char eq, SparseMatrix& aMat, const char op, Vector& bVec)
{
    switch (op) {
    case '/':
        // ret = aMat^{-1} * bVec;
        // aMat is positive definite
        // and already colesky factorized.
        return solveSystems(rVec, aMat, bVec);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(mpfr_class& ret, const char eq, Vector& aVec, const char op, Vector& bVec)
{
    switch (op) {
    case '.':
        return getInnerProduct(ret, aVec, bVec);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(mpfr_class& ret, const char eq, SparseMatrix& aMat, const char op, DenseMatrix& bMat)
{
    switch (op) {
    case '.':
        return getInnerProduct(ret, aMat, bMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

/////////////////////////////////////////////////////////////////////////

bool Lal::getInnerProduct(mpfr_class& ret, DenseLinearSpace& aMat, DenseLinearSpace& bMat)
{
    bool total_judge = _SUCCESS;
    ret = 0.0;
    mpfr_class tmp_ret;

    // for SDP
    if (aMat.SDP_nBlock != bMat.SDP_nBlock) {
        rError("getInnerProduct:: different memory size");
    }
    for (int l = 0; l < aMat.SDP_nBlock; ++l) {
        bool judge = Lal::getInnerProduct(tmp_ret, aMat.SDP_block[l], bMat.SDP_block[l]);
        ret += tmp_ret;
        if (judge == FAILURE) {
            rMessage(" something failed");
            total_judge = FAILURE;
        }
    }

    // for LP
    if (aMat.LP_nBlock != bMat.LP_nBlock) {
        rError("getInnerProduct:: different memory size");
    }
    for (int l = 0; l < aMat.LP_nBlock; ++l) {
        tmp_ret = aMat.LP_block[l] * bMat.LP_block[l];
        ret += tmp_ret;
    }

    return total_judge;
}

bool Lal::getInnerProduct(mpfr_class& ret, SparseLinearSpace& aMat, DenseLinearSpace& bMat)
{
    bool total_judge = _SUCCESS;
    ret = 0.0;
    mpfr_class tmp_ret;

    // for SDP
    for (int l = 0; l < aMat.SDP_sp_nBlock; ++l) {
        int index = aMat.SDP_sp_index[l];
        bool judge = Lal::getInnerProduct(tmp_ret, aMat.SDP_sp_block[l], bMat.SDP_block[index]);
        ret += tmp_ret;
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    for (int l = 0; l < aMat.LP_sp_nBlock; ++l) {
        int index = aMat.LP_sp_index[l];
        tmp_ret = aMat.LP_sp_block[l] * bMat.LP_block[index];
        ret += tmp_ret;
    }

    return total_judge;
}

bool Lal::multiply(DenseLinearSpace& retMat, DenseLinearSpace& aMat, mpfr_class* scalar)
{
    bool total_judge = _SUCCESS;

    // for SDP
    if (retMat.SDP_nBlock != aMat.SDP_nBlock) {
        rError("multiply:: different memory size");
    }
    for (int l = 0; l < aMat.SDP_nBlock; ++l) {
        bool judge = Lal::multiply(retMat.SDP_block[l], aMat.SDP_block[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for LP
    if (retMat.LP_nBlock != aMat.LP_nBlock) {
        rError("multiply:: different memory size");
    }
    for (int l = 0; l < aMat.LP_nBlock; ++l) {
        if (scalar == NULL) {
            retMat.LP_block[l] = aMat.LP_block[l];
        } else {
            retMat.LP_block[l] = aMat.LP_block[l] * (*scalar);
        }
    }

    return total_judge;
}

bool Lal::plus(DenseLinearSpace& retMat, DenseLinearSpace& aMat, DenseLinearSpace& bMat,
               mpfr_class* scalar)
{
    bool total_judge = _SUCCESS;

    // for SDP
    if (retMat.SDP_nBlock != aMat.SDP_nBlock || retMat.SDP_nBlock != bMat.SDP_nBlock) {
        rError("plus:: different nBlock size");
    }
    for (int l = 0; l < retMat.SDP_nBlock; ++l) {
        bool judge = Lal::plus(retMat.SDP_block[l], aMat.SDP_block[l], bMat.SDP_block[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for LP
    if (retMat.LP_nBlock != aMat.LP_nBlock || retMat.LP_nBlock != bMat.LP_nBlock) {
        rError("plus:: different nBlock size");
    }
    for (int l = 0; l < retMat.LP_nBlock; ++l) {
        if (scalar == NULL) {
            retMat.LP_block[l] = aMat.LP_block[l] + bMat.LP_block[l];
        } else {
            retMat.LP_block[l] = aMat.LP_block[l] + bMat.LP_block[l] * (*scalar);
        }
    }

    return total_judge;
}

// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
bool Lal::plus(DenseLinearSpace& retMat, DenseLinearSpace& aMat, SparseLinearSpace& bMat,
               mpfr_class* scalar)
{
    bool total_judge = _SUCCESS;

    // for SDP
    for (int l = 0; l < bMat.SDP_sp_nBlock; ++l) {
        int index = bMat.SDP_sp_index[l];
        bool judge =
            Lal::plus(retMat.SDP_block[index], aMat.SDP_block[index], bMat.SDP_sp_block[l], scalar);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }

    // for LP
    for (int l = 0; l < bMat.LP_sp_nBlock; ++l) {
        int index = bMat.LP_sp_index[l];
        if (scalar == NULL) {
            retMat.LP_block[index] = aMat.LP_block[index] + bMat.LP_sp_block[l];
        } else {
            retMat.LP_block[index] = aMat.LP_block[index] + bMat.LP_sp_block[l] * (*scalar);
        }
    }

    return total_judge;
}

bool Lal::getSymmetrize(DenseLinearSpace& aMat)
{
    bool total_judge = _SUCCESS;
    // for SDP
    for (int l = 0; l < aMat.SDP_nBlock; ++l) {
        bool judge = Lal::getSymmetrize(aMat.SDP_block[l]);
        if (judge == FAILURE) {
            total_judge = FAILURE;
        }
    }
    return total_judge;
}

// ret = a '*' (*scalar)
bool Lal::let(DenseLinearSpace& retMat, const char eq, DenseLinearSpace& aMat, const char op,
              mpfr_class* scalar)
{
    switch (op) {
    case '*':
        return multiply(retMat, aMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool Lal::let(DenseLinearSpace& retMat, const char eq, DenseLinearSpace& aMat, const char op,
              DenseLinearSpace& bMat, mpfr_class* scalar)
{
    mpfr_class minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool Lal::let(DenseLinearSpace& retMat, const char eq, DenseLinearSpace& aMat, const char op,
              SparseLinearSpace& bMat, mpfr_class* scalar)
{
    mpfr_class minus_scalar;
    switch (op) {
    case '+':
        return plus(retMat, aMat, bMat, scalar);
        break;
    case '-':
        if (scalar) {
            minus_scalar = -(*scalar);
            scalar = &minus_scalar;
        } else {
            scalar = &MMONE;
        }
        return plus(retMat, aMat, bMat, scalar);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(mpfr_class& ret, const char eq, DenseLinearSpace& aMat, const char op,
              DenseLinearSpace& bMat)
{
    switch (op) {
    case '.':
        return getInnerProduct(ret, aMat, bMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool Lal::let(mpfr_class& ret, const char eq, SparseLinearSpace& aMat, const char op,
              DenseLinearSpace& bMat)
{
    switch (op) {
    case '.':
        return getInnerProduct(ret, aMat, bMat);
        break;
    default:
        rError("let:: operator error");
        break;
    }
    return FAILURE;
}

} // namespace sdpa

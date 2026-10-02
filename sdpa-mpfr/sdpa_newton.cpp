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

#include <sdpa_newton.h>
#include <sdpa_parts.h>

namespace sdpa {

Newton::Newton(int m, int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock, int* SOCP_blockStruct,
               int LP_nBlock)
{
    initialize(m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
}

Newton::~Newton()
{
    terminate();
}

void Newton::initialize(int m, int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock,
                        int* SOCP_blockStruct, int LP_nBlock)
{
    gVec.initialize(m);

    DxMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    DyVec.initialize(m);
    DzMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    r_zinvMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    x_rd_zinvMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);

    rNewCheck();
    useFormula = new FormulaType[m * SDP_nBlock];
    if (useFormula == NULL) {
        rError("Newton:: memory exhausted ");
    }

    bMat_type = DENSE;

    // Caution: if SDPA doesn't use sparse bMat,
    //          following variables are indefinite.
    this->SDP_nBlock = -1;
    SDP_number = NULL;
    SDP_location_sparse_bMat = NULL;
    SDP_constraint1 = NULL;
    SDP_constraint2 = NULL;
    SDP_blockIndex1 = NULL;
    SDP_blockIndex2 = NULL;
    this->SOCP_nBlock = -1;
    SOCP_number = NULL;
    SOCP_location_sparse_bMat = NULL;
    SOCP_constraint1 = NULL;
    SOCP_constraint2 = NULL;
    SOCP_blockIndex1 = NULL;
    SOCP_blockIndex2 = NULL;
    this->LP_nBlock = -1;
    LP_number = NULL;
    LP_location_sparse_bMat = NULL;
    LP_constraint1 = NULL;
    LP_constraint2 = NULL;
    LP_blockIndex1 = NULL;
    LP_blockIndex2 = NULL;

    ordering = NULL;
    reverse_ordering = NULL;
    diagonalIndex = NULL;
}

void Newton::terminate()
{

    if (bMat_type == SPARSE) {

        if (SDP_location_sparse_bMat && SDP_constraint1 && SDP_constraint2 && SDP_blockIndex1 &&
            SDP_blockIndex2) {
            for (int k = 0; k < SDP_nBlock; ++k) {
                delete[] SDP_location_sparse_bMat[k];
                delete[] SDP_constraint1[k];
                delete[] SDP_constraint2[k];
                delete[] SDP_blockIndex1[k];
                delete[] SDP_blockIndex2[k];
                SDP_location_sparse_bMat[k] = NULL;
                SDP_constraint1[k] = NULL;
                SDP_constraint2[k] = NULL;
                SDP_blockIndex1[k] = NULL;
                SDP_blockIndex2[k] = NULL;
            }
            delete[] SDP_number;
            delete[] SDP_location_sparse_bMat;
            delete[] SDP_constraint1;
            delete[] SDP_constraint2;
            delete[] SDP_blockIndex1;
            delete[] SDP_blockIndex2;
            SDP_number = NULL;
            SDP_location_sparse_bMat = NULL;
            SDP_constraint1 = NULL;
            SDP_constraint2 = NULL;
            SDP_blockIndex1 = NULL;
            SDP_blockIndex2 = NULL;
        }
        if (LP_location_sparse_bMat && LP_constraint1 && LP_constraint2 && LP_blockIndex1 &&
            LP_blockIndex2) {
            for (int k = 0; k < LP_nBlock; ++k) {
                delete[] LP_location_sparse_bMat[k];
                delete[] LP_constraint1[k];
                delete[] LP_constraint2[k];
                delete[] LP_blockIndex1[k];
                delete[] LP_blockIndex2[k];
                LP_location_sparse_bMat[k] = NULL;
                LP_constraint1[k] = NULL;
                LP_constraint2[k] = NULL;
                LP_blockIndex1[k] = NULL;
                LP_blockIndex2[k] = NULL;
            }
            delete[] LP_number;
            delete[] LP_location_sparse_bMat;
            delete[] LP_constraint1;
            delete[] LP_constraint2;
            delete[] LP_blockIndex1;
            delete[] LP_blockIndex2;
            LP_number = NULL;
            LP_location_sparse_bMat = NULL;
            LP_constraint1 = NULL;
            LP_constraint2 = NULL;
            LP_blockIndex1 = NULL;
            LP_blockIndex2 = NULL;
        }

        if (ordering) {
            delete[] ordering;
            ordering = NULL;
        }
        if (reverse_ordering) {
            delete[] reverse_ordering;
            reverse_ordering = NULL;
        }
        if (diagonalIndex) {
            delete[] diagonalIndex;
            diagonalIndex = NULL;
        }
        sparse_bMat.terminate();

    } else { // bMat_type == DENSE
        bMat.terminate();
    }

    gVec.terminate();
    DxMat.terminate();
    DyVec.terminate();
    DzMat.terminate();
    r_zinvMat.terminate();
    x_rd_zinvMat.terminate();

    if (useFormula != NULL) {
        delete[] useFormula;
    }
    useFormula = NULL;
}

void Newton::initialize_dense_bMat(int m)
{
    bMat.initialize(m, m, DenseMatrix::DENSE);
}

// 2008/03/12 kazuhide nakata
void Newton::initialize_sparse_bMat(int m, IV* newToOldIV, IVL* symbfacIVL)
{

    int i, j;
    int* newToOld;

    newToOld = IV_entries(newToOldIV);

    rNewCheck();
    ordering = new int[m];
    if (ordering == NULL) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }
    for (i = 0; i < m; i++) {
        ordering[i] = newToOld[i];
    }

    rNewCheck();
    reverse_ordering = new int[m];
    if (reverse_ordering == NULL) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }
    for (i = 0; i < m; i++) {
        reverse_ordering[ordering[i]] = i;
    }

    // separate front or back node
    int* counter;
    int nClique = IVL_nlist(symbfacIVL);
    int psize;
    int* pivec;
    bool* bnode;
    int* nFront;

    rNewCheck();
    counter = new int[m];
    bnode = new bool[m];
    nFront = new int[nClique];

    if ((counter == NULL) || (bnode == NULL) || (nFront == NULL)) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }

    for (i = 0; i < m; i++) {
        bnode[i] = false;
        counter[i] = -1;
    }

    // search number of front
    for (int l = nClique - 1; l >= 0; l--) {
        IVL_listAndSize(symbfacIVL, l, &psize, &pivec);
        for (i = 0; i < psize; i++) {
            int ii = reverse_ordering[pivec[i]];
            if (bnode[ii] == false) {
                counter[ii] = psize - i;
                bnode[ii] = true;
            } else {
                nFront[l] = i;
                break;
            }
        }
        if (i == psize) {
            nFront[l] = psize;
        }
    }

    // error check
    for (i = 0; i < m; i++) {
        if (counter[i] == -1) {
            rError("Newton::initialize_sparse_bMat: program bug");
        }
    }

    // make index of diagonal
    rNewCheck();
    diagonalIndex = new int[m + 1];
    if (diagonalIndex == NULL) {
        rError("Newton::initialize_sparse_bMat memory exhausted ");
    }

    diagonalIndex[0] = 0;
    for (i = 1; i < m + 1; i++) {
        diagonalIndex[i] = diagonalIndex[i - 1] + counter[i - 1];
    }

    // initialize sparse_bMat
    sparse_bMat.initialize(m, m, SparseMatrix::SPARSE, diagonalIndex[m]);

    // initialize index of sparse_bmat
    int nonzeros = 0;
    for (int l = 0; l < nClique; l++) {
        IVL_listAndSize(symbfacIVL, l, &psize, &pivec);
        for (i = 0; i < nFront[l]; i++) {
            int ii = reverse_ordering[pivec[i]];
            for (j = i; j < psize; j++) {
                int jj = reverse_ordering[pivec[j]];
                int index = diagonalIndex[ii] + j - i;
                sparse_bMat.row_index[index] = ii;
                sparse_bMat.column_index[index] = jj;
                nonzeros++;
            }
        }
    }
    // error check
    if (nonzeros != sparse_bMat.NonZeroNumber) {
        rError("Newton::initialize_sparse_bMat  probram bug");
    }
    sparse_bMat.NonZeroCount = nonzeros;

    delete[] counter;
    delete[] bnode;
    delete[] nFront;
}

// 2008/03/12 kazuhide nakata
void Newton::initialize_bMat(int m, Chordal& chordal, InputData& inputData, FILE* fpOut)
{
    /* Create clique tree */

    switch (chordal.best) {
    case -1: {
        bMat_type = DENSE;
        printf("DENSE computations\n");
        fprintf(fpOut, "DENSE computation\n");
        initialize_dense_bMat(m);
        break;
    }
    case 0: {
        rError("no support for METIS");
        break;
    }
    case 1: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_MMD, chordal.symbfacIVL_MMD);
        make_aggrigateIndex(inputData);
        break;
    }
    case 2: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_ND, chordal.symbfacIVL_ND);
        make_aggrigateIndex(inputData);
        break;
    }
    case 3: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_MS, chordal.symbfacIVL_MS);
        make_aggrigateIndex(inputData);
        break;
    }
    case 4: {
        bMat_type = SPARSE;
        printf("SPARSE computation\n");
        fprintf(fpOut, "SPARSE computation\n");
        initialize_sparse_bMat(m, chordal.newToOldIV_NDMS, chordal.symbfacIVL_NDMS);
        make_aggrigateIndex(inputData);
        break;
    }
    }
}

void Newton::make_aggrigateIndex_SDP(InputData& inputData)
{
    int t, ii, jj;

    SDP_nBlock = inputData.SDP_nBlock;
    rNewCheck();
    SDP_number = new int[SDP_nBlock];
    if (SDP_number == NULL) {
        rError("Newton::make_aggrigateIndex_SDP memory exhausted ");
    }

    // memory allocate for aggrigateIndex
    rNewCheck();
    SDP_constraint1 = new int*[SDP_nBlock];
    SDP_constraint2 = new int*[SDP_nBlock];
    SDP_blockIndex1 = new int*[SDP_nBlock];
    SDP_blockIndex2 = new int*[SDP_nBlock];
    SDP_location_sparse_bMat = new int*[SDP_nBlock];
    if ((SDP_constraint1 == NULL) || (SDP_constraint2 == NULL) || (SDP_blockIndex1 == NULL) ||
        (SDP_blockIndex2 == NULL) || (SDP_location_sparse_bMat == NULL)) {
        rError("Newton::make_aggrigateIndex_SDP memory exhausted ");
    }

    for (int l = 0; l < SDP_nBlock; l++) {
        int tmp = (inputData.SDP_nConstraint[l] + 1) * inputData.SDP_nConstraint[l] / 2;
        rNewCheck();
        SDP_number[l] = tmp;
        SDP_constraint1[l] = new int[tmp];
        SDP_constraint2[l] = new int[tmp];
        SDP_blockIndex1[l] = new int[tmp];
        SDP_blockIndex2[l] = new int[tmp];
        SDP_location_sparse_bMat[l] = new int[tmp];
        if ((SDP_constraint1[l] == NULL) || (SDP_constraint2[l] == NULL) ||
            (SDP_blockIndex1[l] == NULL) || (SDP_blockIndex2[l] == NULL) ||
            (SDP_location_sparse_bMat[l] == NULL)) {
            rError("Newton::make_aggrigateIndex_SDP memory exhausted ");
        }
    }

    for (int l = 0; l < SDP_nBlock; l++) {
        int NonZeroCount = 0;

        for (int k1 = 0; k1 < inputData.SDP_nConstraint[l]; k1++) {
            int i = inputData.SDP_constraint[l][k1];
            int ib = inputData.SDP_blockIndex[l][k1];
            int inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;

            for (int k2 = 0; k2 < inputData.SDP_nConstraint[l]; k2++) {
                int j = inputData.SDP_constraint[l][k2];
                int jb = inputData.SDP_blockIndex[l][k2];
                int jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;

                if ((inz < jnz) || ((inz == jnz) && (i < j))) {
                    continue;
                }

                // set index which A_i and A_j are not zero matrix
                SDP_constraint1[l][NonZeroCount] = i;
                SDP_constraint2[l][NonZeroCount] = j;
                SDP_blockIndex1[l][NonZeroCount] = ib;
                SDP_blockIndex2[l][NonZeroCount] = jb;
                if (reverse_ordering[i] < reverse_ordering[j]) {
                    ii = reverse_ordering[i];
                    jj = reverse_ordering[j];
                } else {
                    jj = reverse_ordering[i];
                    ii = reverse_ordering[j];
                }

                // binary search for index of sparse_bMat
                t = -1;
                int begin = diagonalIndex[ii];
                int end = diagonalIndex[ii + 1] - 1;
                int target = (begin + end) / 2;
                while (end - begin > 1) {
                    if (sparse_bMat.column_index[target] < jj) {
                        begin = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] > jj) {
                        end = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] == jj) {
                        t = target;
                        break;
                    }
                }
                if (t == -1) {
                    if (sparse_bMat.column_index[begin] == jj) {
                        t = begin;
                    } else if (sparse_bMat.column_index[end] == jj) {
                        t = end;
                    } else {
                        rError("Newton::make_aggrigateIndex_SDP  program bug");
                    }
                }

                SDP_location_sparse_bMat[l][NonZeroCount] = t;
                NonZeroCount++;
            }
        } // for k1
    } //for k  kth block
}

void Newton::make_aggrigateIndex_LP(InputData& inputData)
{
    int t, ii, jj;

    LP_nBlock = inputData.LP_nBlock;
    rNewCheck();
    LP_number = new int[LP_nBlock];
    if (LP_number == NULL) {
        rError("Newton::make_aggrigateIndex_LP memory exhausted ");
    }

    // memory allocate for aggrigateIndex
    rNewCheck();
    LP_constraint1 = new int*[LP_nBlock];
    LP_constraint2 = new int*[LP_nBlock];
    LP_blockIndex1 = new int*[LP_nBlock];
    LP_blockIndex2 = new int*[LP_nBlock];
    LP_location_sparse_bMat = new int*[LP_nBlock];
    if ((LP_constraint1 == NULL) || (LP_constraint2 == NULL) || (LP_blockIndex1 == NULL) ||
        (LP_blockIndex2 == NULL) || (LP_location_sparse_bMat == NULL)) {
        rError("Newton::make_aggrigateIndex_LP memory exhausted ");
    }

    for (int l = 0; l < LP_nBlock; l++) {
        int tmp = (inputData.LP_nConstraint[l] + 1) * inputData.LP_nConstraint[l] / 2;
        rNewCheck();
        LP_number[l] = tmp;
        LP_constraint1[l] = new int[tmp];
        LP_constraint2[l] = new int[tmp];
        LP_blockIndex1[l] = new int[tmp];
        LP_blockIndex2[l] = new int[tmp];
        LP_location_sparse_bMat[l] = new int[tmp];
        if ((LP_constraint1[l] == NULL) || (LP_constraint2[l] == NULL) ||
            (LP_blockIndex1[l] == NULL) || (LP_blockIndex2[l] == NULL) ||
            (LP_location_sparse_bMat[l] == NULL)) {
            rError("Newton::make_aggrigateIndex_LP memory exhausted ");
        }
    }

    for (int l = 0; l < LP_nBlock; l++) {
        int NonZeroCount = 0;

        for (int k1 = 0; k1 < inputData.LP_nConstraint[l]; k1++) {
            int i = inputData.LP_constraint[l][k1];
            int ib = inputData.LP_blockIndex[l][k1];

            for (int k2 = 0; k2 < inputData.LP_nConstraint[l]; k2++) {
                int j = inputData.LP_constraint[l][k2];
                int jb = inputData.LP_blockIndex[l][k2];

                if (i < j) {
                    continue;
                }

                // set index which A_i and A_j are not zero matrix
                LP_constraint1[l][NonZeroCount] = i;
                LP_constraint2[l][NonZeroCount] = j;
                LP_blockIndex1[l][NonZeroCount] = ib;
                LP_blockIndex2[l][NonZeroCount] = jb;
                if (reverse_ordering[i] < reverse_ordering[j]) {
                    ii = reverse_ordering[i];
                    jj = reverse_ordering[j];
                } else {
                    jj = reverse_ordering[i];
                    ii = reverse_ordering[j];
                }

                // binary search for index of sparse_bMat
                t = -1;
                int begin = diagonalIndex[ii];
                int end = diagonalIndex[ii + 1] - 1;
                int target = (begin + end) / 2;
                while (end - begin > 1) {
                    if (sparse_bMat.column_index[target] < jj) {
                        begin = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] > jj) {
                        end = target;
                        target = (begin + end) / 2;
                    } else if (sparse_bMat.column_index[target] == jj) {
                        t = target;
                        break;
                    }
                }
                if (t == -1) {
                    if (sparse_bMat.column_index[begin] == jj) {
                        t = begin;
                    } else if (sparse_bMat.column_index[end] == jj) {
                        t = end;
                    } else {
                        rError("Newton::make_aggrigateIndex_SDP  program bug");
                    }
                }

                LP_location_sparse_bMat[l][NonZeroCount] = t;
                NonZeroCount++;
            }
        } // for k1
    } //for k  kth block
}

void Newton::make_aggrigateIndex(InputData& inputData)
{
    make_aggrigateIndex_SDP(inputData);
    make_aggrigateIndex_LP(inputData);
}

void Newton::computeFormula_SDP(InputData& inputData, mpfr_class DenseRatio, mpfr_class Kappa)
{
    int m = inputData.b.nDim;
    int SDP_nBlock = inputData.SDP_nBlock;

    int* upNonZeroCount;
    rNewCheck();
    upNonZeroCount = new int[m * SDP_nBlock];
    if (upNonZeroCount == NULL) {
        rError("Newton:: memory exhausted ");
    }

    // We have no chance to use DenseRatio
    if (upNonZeroCount == NULL || useFormula == NULL) {
        rError("Newton:: failed initialization");
    }

    SparseLinearSpace* A = inputData.A;

    // Count sum of number of elements
    // that each number of elements are less than own.

    for (int iter = 0; iter < m * SDP_nBlock; iter++) {
        upNonZeroCount[iter] = 0;
    }

    for (int l = 0; l < SDP_nBlock; ++l) {
        for (int k1 = 0; k1 < inputData.SDP_nConstraint[l]; k1++) {
            int i = inputData.SDP_constraint[l][k1];
            int ib = inputData.SDP_blockIndex[l][k1];
            int inz = A[i].SDP_sp_block[ib].NonZeroEffect;
            int up = inz;

            for (int k2 = 0; k2 < inputData.SDP_nConstraint[l]; k2++) {
                int j = inputData.SDP_constraint[l][k2];
                int jb = inputData.SDP_blockIndex[l][k2];
                int jnz = A[j].SDP_sp_block[jb].NonZeroEffect;
                if (jnz < inz) {
                    up += jnz;
                } else if ((jnz == inz) && (j < i)) {
                    up += jnz;
                }
            }
            upNonZeroCount[i * SDP_nBlock + l] = up;
        }
    }

    // Determine which formula
    for (int l = 0; l < SDP_nBlock; ++l) {
        int countf1, countf2, countf3;
        countf1 = countf2 = countf3 = 0;
        for (int k = 0; k < inputData.SDP_nConstraint[l]; k++) {
            int i = inputData.SDP_constraint[l][k];
            int ib = inputData.SDP_blockIndex[l][k];
            mpfr_class inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;

            mpfr_class f1, f2, f3;
            mpfr_class n = inputData.A[i].SDP_sp_block[ib].nRow;
            mpfr_class up = upNonZeroCount[i * SDP_nBlock + l];

            f1 = Kappa * n * inz + n * n * n + Kappa * up;
            f2 = Kappa * n * inz + Kappa * (n + 1) * up;
            f3 = Kappa * (2 * Kappa * inz + 1) * up / Kappa;
            if (inputData.A[i].SDP_sp_block[ib].type == SparseMatrix::DENSE) {
                // if DENSE, we use only F1 or F2,
                // that is we don't use F3
                if (f1 < f2) {
                    useFormula[i * SDP_nBlock + l] = F1;
                    countf1++;
                } else {
                    useFormula[i * SDP_nBlock + l] = F2;
                    countf2++;
                }
            } else {
                // this case is SPARSE
                if (f1 < f2 && f1 < f3) {
                    useFormula[i * SDP_nBlock + l] = F1;
                    countf1++;
                } else if (f2 < f3) {
                    useFormula[i * SDP_nBlock + l] = F2;
                    countf2++;
                } else {
                    useFormula[i * SDP_nBlock + l] = F3;
                    countf3++;
                }
            }
        }
    } // end of 'for (int l)'

    if (upNonZeroCount != NULL) {
        delete[] upNonZeroCount;
    }
    upNonZeroCount = NULL;

    return;
}

void Newton::compute_rMat(Newton::WHICH_DIRECTION direction, AverageComplementarity& mu,
                          DirectionParameter& beta, Solutions& currentPt, WorkVariables& work)
{

    //     CORRECTOR ::  r_zinv = (-XZ -dXdZ + mu I)Z^{-1}
    // not CORRECTOR ::  r_zinv = (-XZ + mu I)Z^{-1}
    mpfr_class target = beta.value * mu.current;
    Lal::let(r_zinvMat, '=', currentPt.invzMat, '*', &target);
    Lal::let(r_zinvMat, '=', r_zinvMat, '+', currentPt.xMat, &MMONE);

    if (direction == CORRECTOR) {
        // work.DLS1 = Dx Dz Z^{-1}
        Jal::ns_jordan_triple_product(work.DLS1, DxMat, DzMat, currentPt.invzMat, work.DLS2);
        Lal::let(r_zinvMat, '=', r_zinvMat, '+', work.DLS1, &MMONE);
    }
}

void Newton::Make_gVec(Newton::WHICH_DIRECTION direction, InputData& inputData,
                       Solutions& currentPt, Residuals& currentRes, AverageComplementarity& mu,
                       DirectionParameter& beta, Phase& phase, WorkVariables& work,
                       ComputeTime& com)
{
    TimeStart(START1);
    compute_rMat(direction, mu, beta, currentPt, work);

    TimeEnd(END1);

    com.makerMat += TimeCal(START1, END1);

    TimeStart(START2);
    TimeStart(START_GVEC_MUL);

    // work.DLS1 = R Z^{-1} - X D Z^{-1} = r_zinv - X D Z^{-1}
    if (phase.value == SolveInfo::pFEAS || phase.value == SolveInfo::noINFO) {

        if (direction == CORRECTOR) {
            // x_rd_zinvMat is computed in PREDICTOR step
            Lal::let(work.DLS1, '=', r_zinvMat, '+', x_rd_zinvMat, &MMONE);
        } else {
            // currentPt is infeasilbe, that is the residual
            // dualMat is not 0.
            //      x_rd_zinvMat = X D Z^{-1}
            Jal::ns_jordan_triple_product(x_rd_zinvMat, currentPt.xMat, currentRes.dualMat,
                                          currentPt.invzMat, work.DLS2);
            Lal::let(work.DLS1, '=', r_zinvMat, '+', x_rd_zinvMat, &MMONE);
        } // if (direction == CORRECTOR)

    } else {
        // dualMat == 0
        work.DLS1.copyFrom(r_zinvMat);
    }

    TimeEnd(END_GVEC_MUL);
    com.makegVecMul += TimeCal(START_GVEC_MUL, END_GVEC_MUL);

    inputData.multi_InnerProductToA(work.DLS1, gVec);
    Lal::let(gVec, '=', gVec, '*', &MMONE);

    Lal::let(gVec, '=', gVec, '+', currentRes.primalVec);

    TimeEnd(END2);
    com.makegVec += TimeCal(START2, END2);
}

void Newton::calF1(mpfr_class& ret, DenseMatrix& G, SparseMatrix& Aj)
{
    Lal::let(ret, '=', Aj, '.', G);
}

void Newton::calF2(mpfr_class& ret, DenseMatrix& F, DenseMatrix& G, DenseMatrix& X,
                   SparseMatrix& Aj, bool& hasF2Gcal)
{
    int alpha, beta;
    mpfr_class value1, value2;

    int n = Aj.nRow;
    switch (Aj.type) {
    case SparseMatrix::SPARSE:
        ret = 0.0;
        for (int index = 0; index < Aj.NonZeroCount; ++index) {
            alpha = Aj.row_index[index];
            beta = Aj.column_index[index];
            value1 = Aj.sp_ele[index];

            value2 = Rdot(n, X.de_ele + alpha, n, F.de_ele + (n * beta), 1);
            ret += value1 * value2;
            if (alpha != beta) {
                value2 = Rdot(n, X.de_ele + beta, n, F.de_ele + (n * alpha), 1);
                ret += value1 * value2;
            }
        }
        break;
    case SparseMatrix::DENSE:
        // G is temporary matrix
        if (hasF2Gcal == false) {
            Lal::let(G, '=', X, '*', F);
            hasF2Gcal = true;
        }
        Lal::let(ret, '=', Aj, '.', G);
        break;
    } // end of switch
}

void Newton::calF3(mpfr_class& ret, DenseMatrix& F, DenseMatrix& G, DenseMatrix& X,
                   DenseMatrix& invZ, SparseMatrix& Ai, SparseMatrix& Aj)
{
    // Ai and Aj are SPARSE
    ret = 0.0;
    mpfr_class sum;
    for (int index1 = 0; index1 < Aj.NonZeroCount; ++index1) {
        int alpha = Aj.row_index[index1];
        int beta = Aj.column_index[index1];
        mpfr_class value1 = Aj.sp_ele[index1];
        sum = 0.0;
        for (int index2 = 0; index2 < Ai.NonZeroCount; ++index2) {
            int gamma = Ai.row_index[index2];
            int delta = Ai.column_index[index2];
            mpfr_class value2 = Ai.sp_ele[index2];
            mpfr_class plu =
                value2 * invZ.de_ele[delta + invZ.nCol * beta] * X.de_ele[alpha + X.nCol * gamma];
            sum += plu;
            if (gamma != delta) {
                mpfr_class plu2 = value2 * invZ.de_ele[gamma + invZ.nCol * beta] *
                                  X.de_ele[alpha + X.nCol * delta];
                sum += plu2;
            }
        }
        ret += value1 * sum;
        if (alpha == beta) {
            continue;
        }
        sum = 0.0;
        for (int index2 = 0; index2 < Ai.NonZeroCount; ++index2) {
            int gamma = Ai.row_index[index2];
            int delta = Ai.column_index[index2];
            mpfr_class value2 = Ai.sp_ele[index2];
            mpfr_class plu =
                value2 * invZ.de_ele[delta + invZ.nCol * alpha] * X.de_ele[beta + X.nCol * gamma];
            sum += plu;
            if (gamma != delta) {
                mpfr_class plu2 = value2 * invZ.de_ele[gamma + invZ.nCol * alpha] *
                                  X.de_ele[beta + X.nCol * delta];
                sum += plu2;
            }
        }
        ret += value1 * sum;
    } // end of 'for (index1)'
    return;
}

void Newton::compute_bMat_dense_SDP(InputData& inputData, Solutions& currentPt, WorkVariables& work,
                                    ComputeTime& com)
{
    int m = currentPt.mDim;
    int SDP_nBlock = inputData.SDP_nBlock;

    for (int l = 0; l < SDP_nBlock; ++l) {
        DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
        DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
        DenseMatrix& work1 = work.DLS1.SDP_block[l];
        DenseMatrix& work2 = work.DLS2.SDP_block[l];

        for (int k1 = 0; k1 < inputData.SDP_nConstraint[l]; k1++) {
            int i = inputData.SDP_constraint[l][k1];
            int ib = inputData.SDP_blockIndex[l][k1];
            int inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;
            SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];

            FormulaType formula = useFormula[i * SDP_nBlock + l];
            TimeStart(B_NDIAG_START1);
            TimeStart(B_NDIAG_START2);

            bool hasF2Gcal = false;
            if (formula == F1) {
                Lal::let(work1, '=', Ai, '*', invzMat);
                Lal::let(work2, '=', xMat, '*', work1);
            } else if (formula == F2) {
                Lal::let(work1, '=', Ai, '*', invzMat);
                hasF2Gcal = false;
            }
            TimeEnd(B_NDIAG_END2);
            com.B_PRE += TimeCal(B_NDIAG_START2, B_NDIAG_END2);

            for (int k2 = 0; k2 < inputData.SDP_nConstraint[l]; k2++) {
                int j = inputData.SDP_constraint[l][k2];
                int jb = inputData.SDP_blockIndex[l][k2];
                int jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;
                SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];

                // Select the formula A[i] or the formula A[j].
                // Use formula that has more NonZeroEffects than others.
                // We must calculate i==j.

                if ((inz < jnz) || ((inz == jnz) && (i < j))) {
                    continue;
                }

                mpfr_class value;
                switch (formula) {
                case F1:
                    calF1(value, work2, Aj);
                    break;
                case F2:
                    calF2(value, work1, work2, xMat, Aj, hasF2Gcal);
                    break;
                case F3:
                    calF3(value, work1, work2, xMat, invzMat, Ai, Aj);
                    break;
                } // end of switch
                if (i != j) {
                    bMat.de_ele[i + m * j] += value;
                    bMat.de_ele[j + m * i] += value;
                } else {
                    bMat.de_ele[i + m * i] += value;
                }
            } // end of 'for (int j)'

            TimeEnd(B_NDIAG_END1);
            double t = TimeCal(B_NDIAG_START1, B_NDIAG_END1);
            switch (formula) {
            case F1:
                com.B_F1 += t;
                break;
            case F2:
                com.B_F2 += t;
                break;
            case F3:
                com.B_F3 += t;
                break;
            }
        } // end of 'for (int i)'
    } // end of 'for (int l)'
}

void Newton::compute_bMat_sparse_SDP(InputData& inputData, Solutions& currentPt,
                                     WorkVariables& work, ComputeTime& com)
{
    for (int l = 0; l < SDP_nBlock; ++l) {
        DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
        DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
        DenseMatrix& work1 = work.DLS1.SDP_block[l];
        DenseMatrix& work2 = work.DLS2.SDP_block[l];
        int previous_i = -1;

        for (int iter = 0; iter < SDP_number[l]; iter++) {
            int i = SDP_constraint1[l][iter];
            int ib = SDP_blockIndex1[l][iter];
            SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];
            FormulaType formula = useFormula[i * SDP_nBlock + l];
            bool hasF2Gcal;

            if (i != previous_i) {
                TimeStart(B_NDIAG_START2);

                hasF2Gcal = false;
                if (formula == F1) {
                    Lal::let(work1, '=', Ai, '*', invzMat);
                    Lal::let(work2, '=', xMat, '*', work1);
                } else if (formula == F2) {
                    Lal::let(work1, '=', Ai, '*', invzMat);
                    hasF2Gcal = false;
                }
                TimeEnd(B_NDIAG_END2);
                com.B_PRE += TimeCal(B_NDIAG_START2, B_NDIAG_END2);
            }

            int j = SDP_constraint2[l][iter];
            int jb = SDP_blockIndex2[l][iter];
            SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];

            mpfr_class value;
            switch (formula) {
            case F1:
                calF1(value, work2, Aj);
                break;
            case F2:
                calF2(value, work1, work2, xMat, Aj, hasF2Gcal);
                break;
            case F3:
                calF3(value, work1, work2, xMat, invzMat, Ai, Aj);
                break;
            } // end of switch
            sparse_bMat.sp_ele[SDP_location_sparse_bMat[l][iter]] += value;
            previous_i = i;
        } // end of 'for (int index)'
    } // end of 'for (int l)'
}

void Newton::compute_bMat_dense_LP(InputData& inputData, Solutions& currentPt, WorkVariables& work,
                                   ComputeTime& com)
{
    int m = currentPt.mDim;
    int LP_nBlock = inputData.LP_nBlock;

    TimeEnd(B_DIAG_START1);
    for (int l = 0; l < LP_nBlock; ++l) {
        mpfr_class xMat = currentPt.xMat.LP_block[l];
        mpfr_class invzMat = currentPt.invzMat.LP_block[l];

        for (int k1 = 0; k1 < inputData.LP_nConstraint[l]; k1++) {
            int i = inputData.LP_constraint[l][k1];
            int ib = inputData.LP_blockIndex[l][k1];
            mpfr_class Ai = inputData.A[i].LP_sp_block[ib];

            for (int k2 = k1; k2 < inputData.LP_nConstraint[l]; k2++) {
                int j = inputData.LP_constraint[l][k2];
                int jb = inputData.LP_blockIndex[l][k2];
                mpfr_class Aj = inputData.A[j].LP_sp_block[jb];

                mpfr_class value;
                value = xMat * invzMat * Ai * Aj;

                if (i != j) {
                    bMat.de_ele[i + m * j] += value;
                    bMat.de_ele[j + m * i] += value;
                } else {
                    bMat.de_ele[i + m * i] += value;
                }
            } // end of 'for (int j)'
        } // end of 'for (int i)'
    } // end of 'for (int l)'
    TimeEnd(B_DIAG_END1);
    com.B_DIAG += TimeCal(B_DIAG_START1, B_DIAG_END1);
}

void Newton::compute_bMat_sparse_LP(InputData& inputData, Solutions& currentPt, WorkVariables& work,
                                    ComputeTime& com)
{
    TimeEnd(B_DIAG_START1);
    for (int l = 0; l < LP_nBlock; ++l) {
        mpfr_class xMat = currentPt.xMat.LP_block[l];
        mpfr_class invzMat = currentPt.invzMat.LP_block[l];

        for (int iter = 0; iter < LP_number[l]; iter++) {
            int i = LP_constraint1[l][iter];
            int ib = LP_blockIndex1[l][iter];
            mpfr_class Ai = inputData.A[i].LP_sp_block[ib];

            int j = LP_constraint2[l][iter];
            int jb = LP_blockIndex2[l][iter];
            mpfr_class Aj = inputData.A[j].LP_sp_block[jb];

            mpfr_class value;
            value = xMat * invzMat * Ai * Aj;
            sparse_bMat.sp_ele[LP_location_sparse_bMat[l][iter]] += value;
        } // end of 'for (int iter)
    } // end of 'for (int l)'
    TimeEnd(B_DIAG_END1);
    com.B_DIAG += TimeCal(B_DIAG_START1, B_DIAG_END1);
}

void Newton::Make_bMat(InputData& inputData, Solutions& currentPt, WorkVariables& work,
                       ComputeTime& com)
{
    TimeStart(START3);
    if (bMat_type == SPARSE) {
        // set sparse_bMat zero
        for (int iter = 0; iter < sparse_bMat.NonZeroCount; ++iter) {
            sparse_bMat.sp_ele[iter] = 0.0;
        }
        compute_bMat_sparse_SDP(inputData, currentPt, work, com);
        compute_bMat_sparse_LP(inputData, currentPt, work, com);
    } else {
        bMat.setZero();
        compute_bMat_dense_SDP(inputData, currentPt, work, com);
        compute_bMat_dense_LP(inputData, currentPt, work, com);
    }
    TimeEnd(END3);
    com.makebMat += TimeCal(START3, END3);
}

// nakata 2004/12/01
void Newton::permuteVec(Vector& gVec, Vector& gVec2)
{
    int i, k;
    int mDIM = gVec2.nDim;

    for (k = 0; k < mDIM; k++) {
        i = ordering[k];
        gVec2.ele[k] = gVec.ele[i];
    }
}

// nakata 2004/12/01
void Newton::reverse_permuteVec(Vector& DyVec2, Vector& DyVec)
{
    int i, k;
    int mDIM = DyVec.nDim;

    for (k = 0; k < mDIM; k++) {
        i = ordering[k];
        DyVec.ele[i] = DyVec2.ele[k];
    }
}

bool Newton::compute_DyVec(Newton::WHICH_DIRECTION direction, InputData& inputData,
                           Solutions& currentPt, WorkVariables& work, ComputeTime& com)
{
    if (direction == PREDICTOR) {
        TimeStart(START3_2);

        if (bMat_type == SPARSE) {
            bool ret = Lal::getCholesky(sparse_bMat, diagonalIndex);
            if (ret == FAILURE) {
                return FAILURE;
            }
        } else {
            bool ret = Lal::choleskyFactorWithAdjust(bMat);
            if (ret == FAILURE) {
                return FAILURE;
            }
        }
        TimeEnd(END3_2);
        com.choleskybMat += TimeCal(START3_2, END3_2);
    }
    // bMat is already cholesky factorized.

    TimeStart(START4);
    if (bMat_type == SPARSE) {
        permuteVec(gVec, work.DV1);
        Lal::let(work.DV2, '=', sparse_bMat, '/', work.DV1);
        reverse_permuteVec(work.DV2, DyVec);
    } else {
        Lal::let(DyVec, '=', bMat, '/', gVec);
    }
    TimeEnd(END4);
    com.solve += TimeCal(START4, END4);
    return _SUCCESS;
}

void Newton::compute_DzMat(InputData& inputData, Residuals& currentRes, Phase& phase,
                           ComputeTime& com)
{
    TimeStart(START_SUMDZ);
    inputData.multi_plusToA(DyVec, DzMat);
    Lal::let(DzMat, '=', DzMat, '*', &MMONE);
    if (phase.value == SolveInfo::pFEAS || phase.value == SolveInfo::noINFO) {
        Lal::let(DzMat, '=', DzMat, '+', currentRes.dualMat);
    }
    TimeEnd(END_SUMDZ);
    com.sumDz += TimeCal(START_SUMDZ, END_SUMDZ);
}

void Newton::compute_DxMat(Solutions& currentPt, WorkVariables& work, ComputeTime& com)
{
    TimeStart(START_DX);
    // work.DLS1 = dX dZ Z^{-1}
    Jal::ns_jordan_triple_product(work.DLS1, currentPt.xMat, DzMat, currentPt.invzMat, work.DLS2);
    // dX = R Z^{-1} - dX dZ Z^{-1}
    Lal::let(DxMat, '=', r_zinvMat, '+', work.DLS1, &MMONE);
    TimeEnd(END_DX);
    TimeStart(START_SYMM);
    Lal::getSymmetrize(DxMat);
    TimeEnd(END_SYMM);
    com.makedX += TimeCal(START_DX, END_DX);
    com.symmetriseDx += TimeCal(START_SYMM, END_SYMM);
}

bool Newton::Mehrotra(Newton::WHICH_DIRECTION direction, InputData& inputData, Solutions& currentPt,
                      Residuals& currentRes, AverageComplementarity& mu, DirectionParameter& beta,
                      Switch& reduction, Phase& phase, WorkVariables& work, ComputeTime& com)
{

    Make_gVec(direction, inputData, currentPt, currentRes, mu, beta, phase, work, com);

    if (direction == PREDICTOR) {
        Make_bMat(inputData, currentPt, work, com);
    }

    bool ret = compute_DyVec(direction, inputData, currentPt, work, com);
    if (ret == FAILURE) {
        return FAILURE;
    }

    TimeStart(START5);

    compute_DzMat(inputData, currentRes, phase, com);
    compute_DxMat(currentPt, work, com);

    TimeEnd(END5);
    com.makedXdZ += TimeCal(START5, END5);

    return true;
}

} // namespace sdpa

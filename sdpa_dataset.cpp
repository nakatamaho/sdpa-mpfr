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

#include <sdpa_dataset.h>
#include <sdpa_parts.h>

namespace sdpa {

Solutions::Solutions()
{
    // Nothings needs.
}

Solutions::~Solutions()
{
    terminate();
}

void Solutions::initialize(int m, int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock,
                           int* SOCP_blockStruct, int LP_nBlock, mpfr_class lambda,
                           ComputeTime& com)
{
    mDim = m;
    nDim = 0;
    for (int l = 0; l < SDP_nBlock; ++l) {
        nDim += SDP_blockStruct[l];
    }
    for (int l = 0; l < SOCP_nBlock; ++l) {
        nDim += SOCP_blockStruct[l];
    }
    nDim += LP_nBlock;

    xMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    xMat.setIdentity(lambda);
    zMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    zMat.setIdentity(lambda);
    yVec.initialize(m);
    yVec.setZero();

    invCholeskyX.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    invCholeskyX.setIdentity(1.0 / sqrt(lambda));
    invCholeskyZ.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    invCholeskyZ.setIdentity(1.0 / sqrt(lambda));
    invzMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    invzMat.setIdentity(1.0 / lambda);
}

void Solutions::terminate()
{
    xMat.terminate();
    zMat.terminate();
    yVec.terminate();
    invCholeskyX.terminate();
    invCholeskyZ.terminate();
    invzMat.terminate();
}

bool Solutions::computeInverse(WorkVariables& work, ComputeTime& com)
{

    bool total_judge = _SUCCESS;

    TimeStart(START1_3);
    if (Jal::getInvChol(invCholeskyX, xMat, work.DLS1) == false) {
        total_judge = FAILURE;
    }
    TimeEnd(END1_3);
    com.xMatTime += TimeCal(START1_3, END1_3);

    TimeStart(START1_4);
    if (Jal::getInvCholAndInv(invCholeskyZ, invzMat, zMat, work.DLS2) == false) {
        total_judge = FAILURE;
    }
    TimeEnd(END1_4);
    com.zMatTime += TimeCal(START1_4, END1_4);

    xzMinEigenValue = 1.0;
    return total_judge;
}

bool Solutions::update(StepLength& alpha, Newton& newton, WorkVariables& work, ComputeTime& com)
{

    bool total_judge = _SUCCESS;

    TimeStart(START1_1);
    Lal::let(xMat, '=', xMat, '+', newton.DxMat, &alpha.primal);
    TimeEnd(END1_1);
    com.xMatTime += TimeCal(START1_1, END1_1);
    Lal::let(yVec, '=', yVec, '+', newton.DyVec, &alpha.dual);
    TimeStart(START1_2);
    Lal::let(zMat, '=', zMat, '+', newton.DzMat, &alpha.dual);
    TimeEnd(END1_2);
    com.zMatTime += TimeCal(START1_2, END1_2);

    const mpfr_class cannot_move = 1.0e-4;
    if (alpha.primal < cannot_move && alpha.dual < cannot_move) {
        rMessage("Step length is too small. ");
        return FAILURE;
    }

    total_judge = computeInverse(work, com);

    return total_judge;
}

InputData::InputData()
{
    A = NULL;
    SDP_nBlock = 0;
    SDP_nConstraint = NULL;
    SDP_constraint = NULL;
    SDP_blockIndex = NULL;
    SOCP_nBlock = 0;
    SOCP_nConstraint = NULL;
    SOCP_constraint = NULL;
    SOCP_blockIndex = NULL;
    SDP_nBlock = 0;
    LP_nConstraint = NULL;
    LP_constraint = NULL;
    LP_blockIndex = NULL;
}

InputData::~InputData()
{
    terminate();
}

void InputData::initialize_bVec(int m)
{
    b.initialize(m);
}

void InputData::terminate()
{
    C.terminate();
    if (A) {
        for (int k = 0; k < b.nDim; ++k) {
            A[k].terminate();
        }
        delete[] A;
        A = NULL;
    }
    b.terminate();

    if (SDP_nConstraint && SDP_constraint && SDP_blockIndex) {
        for (int k = 0; k < SDP_nBlock; ++k) {
            delete[] SDP_constraint[k];
            SDP_constraint[k] = NULL;
            delete[] SDP_blockIndex[k];
            SDP_blockIndex[k] = NULL;
        }
        delete[] SDP_nConstraint;
        SDP_nConstraint = NULL;
        delete[] SDP_constraint;
        SDP_constraint = NULL;
        delete[] SDP_blockIndex;
        SDP_blockIndex = NULL;
    }
    if (LP_nConstraint && LP_constraint && LP_blockIndex) {
        for (int k = 0; k < LP_nBlock; ++k) {
            delete[] LP_constraint[k];
            LP_constraint[k] = NULL;
            delete[] LP_blockIndex[k];
            LP_blockIndex[k] = NULL;
        }
        delete[] LP_nConstraint;
        LP_nConstraint = NULL;
        delete[] LP_constraint;
        LP_constraint = NULL;
        delete[] LP_blockIndex;
        LP_blockIndex = NULL;
    }
}

void InputData::initialize_index_SDP(int SDP_nBlock, ComputeTime& com)
{
    int i, k;
    int mDim = b.nDim;
    int index;
    int* SDP_count;

    this->SDP_nBlock = SDP_nBlock;
    rNewCheck();
    SDP_nConstraint = new int[SDP_nBlock];
    if (SDP_nConstraint == NULL) {
        rError("InputData::initialize_index memory exhauseted ");
    }

    // count non-zero block matrix of A
    for (i = 0; i < SDP_nBlock; i++) {
        SDP_nConstraint[i] = 0;
    }
    for (i = 0; i < mDim; i++) {
        for (k = 0; k < A[i].SDP_sp_nBlock; k++) {
            index = A[i].SDP_sp_index[k];
            SDP_nConstraint[index]++;
        }
    }

    // malloc SDP_constraint, SDP_blockIndex
    rNewCheck();
    SDP_constraint = new int*[SDP_nBlock];
    if (SDP_constraint == NULL) {
        rError("InputData::initialize_index memory exhauseted ");
    }
    for (i = 0; i < SDP_nBlock; i++) {
        rNewCheck();
        SDP_constraint[i] = new int[SDP_nConstraint[i]];
        if (SDP_constraint[i] == NULL) {
            rError("InputData::initialize_index memory exhauseted ");
        }
    }
    rNewCheck();
    SDP_blockIndex = new int*[SDP_nBlock];
    if (SDP_blockIndex == NULL) {
        rError("InputData::initialize_index memory exhauseted ");
    }
    for (i = 0; i < SDP_nBlock; i++) {
        rNewCheck();
        SDP_blockIndex[i] = new int[SDP_nConstraint[i]];
        if (SDP_blockIndex[i] == NULL) {
            rError("InputData::initialize_index memory exhauseted ");
        }
    }

    // input index of non-zero block matrix of A
    rNewCheck() SDP_count = new int[SDP_nBlock];
    if (SDP_count == NULL) {
        rError("InputData::initialize_index memory exhauseted ");
    }
    for (i = 0; i < SDP_nBlock; i++) {
        SDP_count[i] = 0;
    }
    for (i = 0; i < mDim; i++) {
        for (k = 0; k < A[i].SDP_sp_nBlock; k++) {
            index = A[i].SDP_sp_index[k];
            SDP_constraint[index][SDP_count[index]] = i;
            SDP_blockIndex[index][SDP_count[index]] = k;
            SDP_count[index]++;
        }
    }

    if (SDP_count) {
        delete[] SDP_count;
        SDP_count = NULL;
    }
}

void InputData::initialize_index_LP(int LP_nBlock, ComputeTime& com)
{
    int i, k;
    int mDim = b.nDim;
    int index;
    int* LP_count;

    // for LP
    this->LP_nBlock = LP_nBlock;
    rNewCheck();
    LP_nConstraint = new int[LP_nBlock];
    if (LP_nConstraint == NULL) {
        rError("rInputData::initialize_index memory exhauseted ");
    }

    // count non-zero block matrix of A
    for (i = 0; i < LP_nBlock; i++) {
        LP_nConstraint[i] = 0;
    }
    for (i = 0; i < mDim; i++) {
        for (k = 0; k < A[i].LP_sp_nBlock; k++) {
            index = A[i].LP_sp_index[k];
            LP_nConstraint[index]++;
        }
    }

    // malloc LP_constraint, LP_blockIndex
    rNewCheck();
    LP_constraint = new int*[LP_nBlock];
    if (LP_constraint == NULL) {
        rError("InputData::initialize_index memory exhauseted ");
    }
    for (i = 0; i < LP_nBlock; i++) {
        rNewCheck();
        LP_constraint[i] = new int[LP_nConstraint[i]];
        if (LP_constraint[i] == NULL) {
            rError("InputData::initialize_index memory exhauseted ");
        }
    }
    rNewCheck();
    LP_blockIndex = new int*[LP_nBlock];
    if (LP_blockIndex == NULL) {
        rError("InputData::initialize_index memory exhauseted ");
    }
    for (i = 0; i < LP_nBlock; i++) {
        rNewCheck();
        LP_blockIndex[i] = new int[LP_nConstraint[i]];
        if (LP_blockIndex[i] == NULL) {
            rError("InputData::initialize_index memory exhauseted ");
        }
    }

    // input index of non-zero block matrix of A
    rNewCheck();
    LP_count = new int[LP_nBlock];
    if (LP_count == NULL) {
        rError("InputData::initialize_index memory exhauseted ");
    }
    for (i = 0; i < LP_nBlock; i++) {
        LP_count[i] = 0;
    }
    for (i = 0; i < mDim; i++) {
        for (k = 0; k < A[i].LP_sp_nBlock; k++) {
            index = A[i].LP_sp_index[k];
            LP_constraint[index][LP_count[index]] = i;
            LP_blockIndex[index][LP_count[index]] = k;
            LP_count[index]++;
        }
    }

    if (LP_count) {
        delete[] LP_count;
        LP_count = NULL;
    }
}

void InputData::initialize_index(int SDP_nBlock, int SOCP_nBlock, int LP_nBlock, ComputeTime& com)
{
    initialize_index_SDP(SDP_nBlock, com);
    initialize_index_LP(LP_nBlock, com);
}

//   retVec_i := A_i bullet xMat (for i)
void InputData::multi_InnerProductToA(DenseLinearSpace& xMat, Vector& retVec)
{
    mpfr_class ip;

    retVec.setZero();
    for (int i = 0; i < retVec.nDim; i++) {
        Lal::let(ip, '=', A[i], '.', xMat);
        retVec.ele[i] = ip;
    }
}

//   retMat := \sum_{i} A_i xVec_i
void InputData::multi_plusToA(Vector& xVec, DenseLinearSpace& retMat)
{
    retMat.setZero();
    for (int i = 0; i < xVec.nDim; i++) {
        Lal::let(retMat, '=', retMat, '+', A[i], &xVec.ele[i]);
    }
}

Residuals::Residuals()
{
    normPrimalVec = 0.0;
    normDualMat = 0.0;
    centerNorm = 0.0;
}

Residuals::Residuals(int m, int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock,
                     int* SOCP_blockStruct, int LP_nBlock, InputData& inputData,
                     Solutions& currentPt)
{
    initialize(m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock, inputData,
               currentPt);
}

Residuals::~Residuals()
{
    terminate();
}

void Residuals::initialize(int m, int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock,
                           int* SOCP_blockStruct, int LP_nBlock, InputData& inputData,
                           Solutions& currentPt)
{
    primalVec.initialize(m);
    dualMat.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    compute(m, inputData, currentPt);
}

void Residuals::terminate()
{
    primalVec.terminate();
    dualMat.terminate();
}

void Residuals::copyFrom(Residuals& other)
{
    if (this == &other) {
        return;
    }
    primalVec.copyFrom(other.primalVec);
    dualMat.copyFrom(other.dualMat);
    normPrimalVec = other.normPrimalVec;
    normDualMat = other.normDualMat;
    centerNorm = other.centerNorm;
}

mpfr_class Residuals::computeMaxNorm(Vector& primalVec)
{
    mpfr_class ret = 0.0;
    for (int k = 0; k < primalVec.nDim; ++k) {
        mpfr_class tmp = abs(primalVec.ele[k]);
        if (tmp > ret) {
            ret = tmp;
        }
    }
    return ret;
}

mpfr_class Residuals::computeMaxNorm(DenseLinearSpace& dualMat)
{
    int SDP_nBlock = dualMat.SDP_nBlock;
    int SOCP_nBlock = dualMat.SOCP_nBlock;
    int LP_nBlock = dualMat.LP_nBlock;
    mpfr_class ret = 0.0;
    mpfr_class tmp;

    for (int l = 0; l < SDP_nBlock; ++l) {
        mpfr_class* target = dualMat.SDP_block[l].de_ele;
        int size = dualMat.SDP_block[l].nRow;
        for (int j = 0; j < size * size; ++j) {
            tmp = abs(target[j]);
            if (tmp > ret) {
                ret = tmp;
            }
        }
    }

    for (int l = 0; l < SOCP_nBlock; ++l) {
        rError("dataset:: current version do not support SOCP");
    }

    for (int l = 0; l < LP_nBlock; ++l) {
        tmp = abs(dualMat.LP_block[l]);
        if (tmp > ret) {
            ret = tmp;
        }
    }

    return ret;
}

void Residuals::update(int m, InputData& inputData, Solutions& currentPt, ComputeTime& com)
{
    TimeStart(UPDATE_START);
    compute(m, inputData, currentPt);
    TimeEnd(UPDATE_END);
    com.updateRes += TimeCal(UPDATE_START, UPDATE_END);
}

void Residuals::compute(int m, InputData& inputData, Solutions& currentPt)
{
    // p[k] = b[k] - A[k].X;
    inputData.multi_InnerProductToA(currentPt.xMat, primalVec);
    Lal::let(primalVec, '=', primalVec, '*', &MMONE);
    Lal::let(primalVec, '=', primalVec, '+', inputData.b);

    // D = C - Z - \sum A[k]y[k]
    inputData.multi_plusToA(currentPt.yVec, dualMat);
    Lal::let(dualMat, '=', dualMat, '*', &MMONE);
    Lal::let(dualMat, '=', dualMat, '+', inputData.C);
    Lal::let(dualMat, '=', dualMat, '-', currentPt.zMat);

    normPrimalVec = computeMaxNorm(primalVec);
    normDualMat = computeMaxNorm(dualMat);
    centerNorm = 0.0;
}

WorkVariables::WorkVariables()
{
    // Nothings needs.
}

WorkVariables::~WorkVariables()
{
    terminate();
}

void WorkVariables::initialize(int m, int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock,
                               int* SOCP_blockStruct, int LP_nBlock)
{
    DLS1.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    DLS2.initialize(SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct, LP_nBlock);
    DV1.initialize(m);
    DV1.initialize(m);

    if (SDP_nBlock > 0) {
        SDP_BV1.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV2.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV3.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV4.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV5.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV6.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV7.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV8.initialize(SDP_nBlock, SDP_blockStruct);
        SDP_BV9.initialize(SDP_nBlock, SDP_blockStruct);

        int* workStruct = NULL;
        rNewCheck();
        workStruct = new int[SDP_nBlock];
        if (workStruct == NULL) {
            rMessage("WorkVariables :: memory exhausted");
        }

        for (int l = 0; l < SDP_nBlock; ++l) {
            workStruct[l] = max(1, 3 * SDP_blockStruct[l] - 1);
        }
        SDP2_BV1.initialize(SDP_nBlock, workStruct);
        if (workStruct) {
            delete[] workStruct;
            workStruct = NULL;
        }
    }
}

void WorkVariables::terminate()
{
    DLS1.terminate();
    DLS2.terminate();
    DV1.terminate();
    DV1.terminate();

    SDP_BV1.terminate();
    SDP_BV2.terminate();
    SDP_BV3.terminate();
    SDP_BV4.terminate();
    SDP_BV5.terminate();
    SDP_BV6.terminate();
    SDP_BV7.terminate();
    SDP_BV8.terminate();
    SDP_BV9.terminate();
    SDP2_BV1.terminate();
};

} // namespace sdpa

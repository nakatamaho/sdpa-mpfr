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

#include <sdpa_struct.h>

#define sdpa_dset(dset_length, dset_value, dset_pointer, dset_step)                                \
    for (int dset_i = 0, dset_index = 0; dset_i < dset_length; ++dset_i) {                         \
        dset_pointer[dset_index] = dset_value;                                                     \
        dset_index += dset_step;                                                                   \
    }

namespace sdpa {

Vector::Vector()
{
    nDim = 0;
    ele = NULL;
}

Vector::~Vector()
{
    terminate();
}

void Vector::initialize(int nDim, mpfr_class value)
{
    if (ele && this->nDim != nDim) {
        if (ele) {
            delete[] ele;
            ele = NULL;
        }
        if (nDim <= 0) {
            rError("Vector:: nDim is nonpositive");
        }
    }
    this->nDim = nDim;
    if (ele == NULL) {
        ele = NULL;
        rNewCheck();
        ele = new mpfr_class[nDim];
        if (ele == NULL) {
            rError("Vector:: memory exhausted");
        }
    }
    sdpa_dset(nDim, value, ele, IONE);
}

void Vector::initialize(mpfr_class value)
{
    if (nDim <= 0) {
        rError("Vector:: nDim is nonpositive");
    }
    if (ele == NULL) {
        rNewCheck();
        ele = new mpfr_class[nDim];
        if (ele == NULL) {
            rError("Vector:: memory exhausted");
        }
    }
    sdpa_dset(nDim, value, ele, IONE);
}

void Vector::terminate()
{
    if (ele) {
        delete[] ele;
    }
    ele = NULL;
}

void Vector::setZero()
{
    mpfr_class zero = 0.0;
    initialize(zero);
}

void Vector::display(FILE* fpout, mpfr_class scalar)
{
    if (fpout == NULL) {
        return;
    }
    fprintf(fpout, "{");
    for (int j = 0; j < nDim - 1; ++j) {
        mpfr_class mtmp = ele[j] * scalar;
        mpfr_fprintf(fpout, P_FORMAT ",", mtmp.get_mpfr_t());
    }
    if (nDim > 0) {
        mpfr_class mtmp = ele[nDim - 1] * scalar;
        mpfr_fprintf(fpout, P_FORMAT "}\n", mtmp.get_mpfr_t());
    } else {
        fprintf(fpout, "  }\n");
    }
}

bool Vector::copyFrom(Vector& other)
{
    if (this == &other) {
        return _SUCCESS;
    }
    if (nDim != other.nDim && ele != NULL) {
        delete[] ele;
        ele = NULL;
    }
    nDim = other.nDim;
    if (nDim <= 0) {
        rError("Vector:: nDim is nonpositive");
    }
    if (ele == NULL) {
        rNewCheck();
        ele = new mpfr_class[nDim];
        if (ele == NULL) {
            rError("Vector:: memory exhausted");
        }
    }
    Rcopy(nDim, other.ele, 1, ele, 1);
    return _SUCCESS;
}

BlockVector::BlockVector()
{
    nBlock = 0;
    blockStruct = NULL;
    ele = NULL;
}

BlockVector::~BlockVector()
{
    terminate();
}

void BlockVector::initialize(int nBlock, int* blockStruct, mpfr_class value)
{
    this->nBlock = nBlock;
    if (nBlock <= 0) {
        rError("BlockVector:: nBlock is nonpositive");
    }
    this->blockStruct = NULL;
    rNewCheck();
    this->blockStruct = new int[nBlock];
    if (this->blockStruct == NULL) {
        rError("BlockVector:: memory exhausted");
    }
    for (int l = 0; l < nBlock; ++l) {
        this->blockStruct[l] = blockStruct[l];
    }

    ele = NULL;
    rNewCheck();
    ele = new Vector[nBlock];
    if (ele == NULL) {
        rError("BlockVector:: memory exhausted");
    }
    for (int l = 0; l < nBlock; ++l) {
        int size = blockStruct[l];
        if (size < 0) {
            size = -size;
        }
        ele[l].initialize(size, value);
    }
}

void BlockVector::terminate()
{
    if (ele && blockStruct && nBlock >= 0) {
        for (int l = 0; l < nBlock; ++l) {
            ele[l].terminate();
        }
        delete[] ele;
        ele = NULL;

        delete[] blockStruct;
        blockStruct = NULL;
    }
}

SparseMatrix::SparseMatrix()
{
    nRow = 0;
    nCol = 0;
    type = SPARSE;

    NonZeroNumber = 0;

    de_ele = NULL;

    row_index = NULL;
    column_index = NULL;
    sp_ele = NULL;
    NonZeroCount = 0;
    NonZeroEffect = 0;
}

SparseMatrix::~SparseMatrix()
{
    terminate();
}

void SparseMatrix::initialize(int nRow, int nCol, SparseMatrix::Type type, int NonZeroNumber)
{

    SparseMatrix();
    if (nRow <= 0 || nCol <= 0) {
        rError("SparseMatrix:: Dimensions are nonpositive");
    }
    this->nRow = nRow;
    this->nCol = nCol;
    this->type = type;

    int length;
    switch (type) {
    case SPARSE:
        this->NonZeroNumber = NonZeroNumber;
        this->NonZeroCount = 0;
        this->NonZeroEffect = 0;
        if (NonZeroNumber > 0) {
            rNewCheck();
            row_index = new int[NonZeroNumber];
            rNewCheck();
            column_index = new int[NonZeroNumber];
            rNewCheck();
            sp_ele = new mpfr_class[NonZeroNumber];
            if (row_index == NULL || column_index == NULL || sp_ele == NULL) {
                rError("SparseMatrix:: memory exhausted");
            }
        }
        break;
    case DENSE:
        this->NonZeroNumber = nRow * nCol;
        this->NonZeroCount = nRow * nCol;
        this->NonZeroEffect = nRow * nCol;
        rNewCheck();
        de_ele = new mpfr_class[NonZeroNumber];
        if (de_ele == NULL) {
            rError("SparseMatrix:: memory exhausted");
        }
        length = nRow * nCol;
        sdpa_dset(length, MZERO, de_ele, IONE);
        // all elements are 0.
        break;
    }
}

void SparseMatrix::terminate()
{
    if (de_ele) {
        delete[] de_ele;
        de_ele = NULL;
    }
    if (row_index) {
        delete[] row_index;
        row_index = NULL;
    }
    if (column_index) {
        delete[] column_index;
        column_index = NULL;
    }
    if (sp_ele) {
        delete[] sp_ele;
        sp_ele = NULL;
    }
}

void SparseMatrix::changeToDense(bool forceChange)
{
    if (type != SPARSE) {
        return;
    }
    if (forceChange == false && NonZeroCount < (nRow * nCol) * 0.20) {
        // if the number of elements are less than 20 percent,
        // we don't change to Dense.
        return;
    }
    type = DENSE;
    de_ele = NULL;
    int length = nRow * nCol;
    rNewCheck();
    de_ele = new mpfr_class[length];
    if (de_ele == NULL) {
        rError("SparseMatrix:: memory exhausted");
    }
    sdpa_dset(length, MZERO, de_ele, 1);
    // all elements are set 0.
    for (int index = 0; index < NonZeroCount; ++index) {
        int i = row_index[index];
        int j = column_index[index];
        mpfr_class value = sp_ele[index];
        if (i == j) {
            de_ele[i + nCol * j] = value;
        } else {
            de_ele[i + nCol * j] = de_ele[j + nCol * i] = value;
        }
    }
    NonZeroCount = NonZeroNumber = NonZeroEffect = length;
    delete[] row_index;
    delete[] column_index;
    delete[] sp_ele;
    row_index = NULL;
    column_index = NULL;
    sp_ele = NULL;
}

DenseMatrix::DenseMatrix()
{
    nRow = 0;
    nCol = 0;
    type = DENSE;

    de_ele = NULL;
}

DenseMatrix::~DenseMatrix()
{
    terminate();
}

void DenseMatrix::initialize(int nRow, int nCol, DenseMatrix::Type type)
{

    DenseMatrix();
    if (nRow <= 0 || nCol <= 0) {
        rError("DenseMatrix:: Dimensions are nonpositive");
    }
    int old_length = this->nRow * this->nCol;
    this->nRow = nRow;
    this->nCol = nCol;

    int length;
    switch (type) {
    case DENSE:
        length = nRow * nCol;
        if (de_ele && old_length != length) {
            delete[] de_ele;
            de_ele = NULL;
        }
        if (de_ele == NULL) {
            rNewCheck();
            de_ele = new mpfr_class[length];
            if (de_ele == NULL) {
                rError("DenseMatrix:: memory exhausted");
            }
        }
        sdpa_dset(length, MZERO, de_ele, IONE);
        break;
    case COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
}

void DenseMatrix::terminate()
{
    if (de_ele) {
        delete[] de_ele;
        de_ele = NULL;
    }
}

void DenseMatrix::display(FILE* fpout)
{
    if (fpout == NULL) {
        return;
    }
    switch (type) {
    case DENSE:
        fprintf(fpout, "{");
        for (int i = 0; i < nRow - 1; ++i) {
            if (i == 0) {
                fprintf(fpout, " ");
            } else {
                fprintf(fpout, "  ");
            }
            fprintf(fpout, "{");
            for (int j = 0; j < nCol - 1; ++j) {
                mpfr_fprintf(fpout, P_FORMAT ",", de_ele[i + nCol * j].get_mpfr_t());
            }
            mpfr_fprintf(fpout, P_FORMAT " },\n", de_ele[i + nCol * (nCol - 1)].get_mpfr_t());
        }
        if (nRow > 1) {
            fprintf(fpout, "  {");
        }
        for (int j = 0; j < nCol - 1; ++j) {
            mpfr_fprintf(fpout, P_FORMAT ",", de_ele[(nRow - 1) + nCol * j].get_mpfr_t());
        }
        mpfr_fprintf(fpout, P_FORMAT " }", de_ele[(nRow - 1) + nCol * (nCol - 1)].get_mpfr_t());
        if (nRow > 1) {
            fprintf(fpout, "   }\n");
        } else {
            fprintf(fpout, "\n");
        }
        break;
    case COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
}

bool DenseMatrix::copyFrom(DenseMatrix& other)
{
    if (this == &other) {
        return _SUCCESS;
    }
    int length;
    switch (other.type) {
    case DENSE:
        type = DENSE;
        if (de_ele && (other.nRow != nRow || other.nCol != nCol)) {
            delete[] de_ele;
            de_ele = NULL;
        }
        nRow = other.nRow;
        nCol = other.nCol;
        if (de_ele == NULL) {
            rNewCheck();
            de_ele = new mpfr_class[nRow * nCol];
            if (de_ele == NULL) {
                rError("DenseMatrix:: memory exhausted");
            }
        }
        length = nRow * nCol;
        Rcopy(length, other.de_ele, 1, de_ele, 1);
        break;
    case COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
    return _SUCCESS;
}

void DenseMatrix::setZero()
{
    int length;
    switch (type) {
    case DENSE:
        length = nRow * nCol;
        sdpa_dset(length, MZERO, de_ele, IONE);
        break;
    case COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
}

void DenseMatrix::setIdentity(mpfr_class scalar)
{
    if (nRow != nCol) {
        rError("SparseMatrix:: Identity matrix must be square matrix");
    }
    int length, step;
    switch (type) {
    case DENSE:
        length = nRow * nCol;
        sdpa_dset(length, MZERO, de_ele, IONE);
        step = nCol + 1;
        sdpa_dset(nCol, scalar, de_ele, step);
        break;
    case COMPLETION:
        rError("DenseMatrix:: no support for COMPLETION");
        break;
    }
}

SparseLinearSpace::SparseLinearSpace()
{
    SDP_sp_nBlock = 0;
    SDP_sp_index = NULL;
    SDP_sp_block = NULL;
    SOCP_sp_nBlock = 0;
    SOCP_sp_index = NULL;
    SOCP_sp_block = NULL;
    LP_sp_nBlock = 0;
    LP_sp_index = NULL;
    LP_sp_block = NULL;
}

SparseLinearSpace::~SparseLinearSpace()
{
    terminate();
}

// sparse form of block index      2008/02/27 kazuhide nakata
void SparseLinearSpace::initialize(int SDP_sp_nBlock, int* SDP_sp_index, int* SDP_sp_blockStruct,
                                   int* SDP_sp_NonZeroNumber, int SOCP_sp_nBlock,
                                   int* SOCP_sp_index, int* SOCP_sp_blockStruct,
                                   int* SOCP_sp_NonZeroNumber, int LP_sp_nBlock, int* LP_sp_index)
{

    // for SDP
    this->SDP_sp_nBlock = SDP_sp_nBlock;
    if (SDP_sp_nBlock > 0) {
        this->SDP_sp_index = NULL;
        rNewCheck();
        this->SDP_sp_index = new int[SDP_sp_nBlock];
        if (this->SDP_sp_index == NULL) {
            rError("SDP_index:: memory exhausted");
        }
        this->SDP_sp_block = NULL;
        rNewCheck();
        this->SDP_sp_block = new SparseMatrix[SDP_sp_nBlock];
        if (this->SDP_sp_block == NULL) {
            rError("SparseLinearSpace:: memory exhausted");
        }
    }
    for (int l = 0; l < SDP_sp_nBlock; ++l) {
        this->SDP_sp_index[l] = SDP_sp_index[l];
        int size = SDP_sp_blockStruct[l];
        SDP_sp_block[l].initialize(size, size, SparseMatrix::SPARSE, SDP_sp_NonZeroNumber[l]);
    }

    // for LP
    this->LP_sp_nBlock = LP_sp_nBlock;
    if (LP_sp_nBlock > 0) {
        this->LP_sp_index = NULL;
        rNewCheck();
        this->LP_sp_index = new int[LP_sp_nBlock];
        if (this->LP_sp_index == NULL) {
            rError("LP_index:: memory exhausted");
        }
        this->LP_sp_block = NULL;
        rNewCheck();
        this->LP_sp_block = new mpfr_class[LP_sp_nBlock];
        if (this->LP_sp_block == NULL) {
            rError("SparseLinearSpace:: memory exhausted");
        }
    }
    for (int l = 0; l < LP_sp_nBlock; ++l) {
        this->LP_sp_index[l] = LP_sp_index[l];
    }
}

void SparseLinearSpace::terminate()
{
    // for SDP
    if (SDP_sp_block && SDP_sp_index && SDP_sp_nBlock >= 0) {
        for (int l = 0; l < SDP_sp_nBlock; ++l) {
            SDP_sp_block[l].terminate();
        }
        delete[] SDP_sp_block;
        SDP_sp_block = NULL;

        delete[] SDP_sp_index;
        SDP_sp_index = NULL;
    }
    // for LP
    if (LP_sp_block && LP_sp_index && LP_sp_nBlock >= 0) {
        delete[] LP_sp_block;
        LP_sp_block = NULL;
        delete[] LP_sp_index;
        LP_sp_index = NULL;
    }
}

void SparseLinearSpace::changeToDense(bool forceChange)
{
    if (SDP_sp_nBlock > 0 && SDP_sp_index && SDP_sp_block) {
        for (int l = 0; l < SDP_sp_nBlock; ++l) {
            SDP_sp_block[l].changeToDense(forceChange);
        }
    }
}

void SparseLinearSpace::setElement_SDP(int block, int i, int j, mpfr_class ele)
{
    int k;

    // seek block
    for (k = 0; k < SDP_sp_nBlock; k++) {
        if (SDP_sp_index[k] == block) {
            break;
        }
    }
    if (k == SDP_sp_nBlock) {
        rError("SparseLinearSpace::setElement no block");
    }

    // check range
    if (SDP_sp_block[k].NonZeroCount >= SDP_sp_block[k].NonZeroNumber) {
        rError("SparseLinearSpace::setElement NonZeroCount >= NonZeroNumber");
    }
    if ((i >= SDP_sp_block[k].nRow) || (j >= SDP_sp_block[k].nCol)) {
        rError("out of range in input data");
    }

    // set element
    int count = SDP_sp_block[k].NonZeroCount;
    SDP_sp_block[k].row_index[count] = i;
    SDP_sp_block[k].column_index[count] = j;
    SDP_sp_block[k].sp_ele[count] = ele;
    SDP_sp_block[k].NonZeroCount++;
    if (i == j) {
        SDP_sp_block[k].NonZeroEffect++;
    } else {
        SDP_sp_block[k].NonZeroEffect += 2;
    }
}

void SparseLinearSpace::setElement_LP(int block, mpfr_class ele)
{
    int k;

    for (k = 0; k < LP_sp_nBlock; k++) {
        if (LP_sp_index[k] == block) {
            break;
        }
    }
    if (k == LP_sp_nBlock) {
        rError("SparseLinearSpace::setElement no block");
    }
    LP_sp_block[k] = ele;
}

DenseLinearSpace::DenseLinearSpace()
{
    SDP_nBlock = 0;
    SDP_block = NULL;
    SOCP_nBlock = 0;
    SOCP_block = NULL;
    LP_nBlock = 0;
    LP_block = NULL;
}

DenseLinearSpace::~DenseLinearSpace()
{
    terminate();
}

void DenseLinearSpace::initialize(int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock,
                                  int* SOCP_blockStruct, int LP_nBlock)
{
    if (SDP_nBlock + SOCP_nBlock + LP_nBlock <= 0) {
        rError("DenseLinearSpace:: SDP + SOCP + LP Block is nonpositive");
    }

    // for SDP
    this->SDP_nBlock = SDP_nBlock;
    if (SDP_nBlock < 0) {
        rError("DenseLinearSpace:: SDP_nBlock is negative");
    }
    if ((SDP_nBlock > 0) && (SDP_block == NULL)) {
        rNewCheck();
        SDP_block = new DenseMatrix[SDP_nBlock];
        if (SDP_block == NULL) {
            rError("DenseLinearSpace:: memory exhausted");
        }
    }
    for (int l = 0; l < SDP_nBlock; ++l) {
        int size = SDP_blockStruct[l];
        if (size > 0) {
            SDP_block[l].initialize(size, size, DenseMatrix::DENSE);
        } else {
            rError("DenseLinearSpace:: SDP size is nonpositive");
        }
    }

    // for SOCP
    this->SOCP_nBlock = 0;
    this->SOCP_block = NULL;

    // for LP
    this->LP_nBlock = LP_nBlock;
    if (LP_nBlock < 0) {
        rError("DenseLinearSpace:: LP_nBlock is negative");
    }
    if ((LP_nBlock > 0) && (LP_block == NULL)) {
        rNewCheck();
        LP_block = new mpfr_class[LP_nBlock];
        if (LP_block == NULL) {
            rError("DenseLinearSpace:: memory exhausted");
        }
    }
    for (int l = 0; l < LP_nBlock; ++l) {
        LP_block[l] = 0.0;
    }
}

void DenseLinearSpace::terminate()
{
    // for SDP
    if (SDP_block && SDP_nBlock > 0) {
        for (int l = 0; l < SDP_nBlock; ++l) {
            SDP_block[l].terminate();
        }
        delete[] SDP_block;
        SDP_block = NULL;
    }

    // SOCP

    // LP
    if (LP_block && LP_nBlock > 0) {
        delete[] LP_block;
        LP_block = NULL;
    }
}

bool DenseLinearSpace::copyFrom(DenseLinearSpace& other)
{
    if (this == &other) {
        return _SUCCESS;
    }

    if (other.SDP_nBlock + other.SOCP_nBlock + other.LP_nBlock <= 0) {
        rError("DenseLinearSpace:: SDP + SOCP + LP Block is nonpositive");
    }
    bool total_judge = _SUCCESS;

    // for SDP
    if (other.SDP_nBlock < 0) {
        rError("DenseLinearSpace:: SDP_nBlock is negative");
    }
    if (SDP_nBlock != other.SDP_nBlock) {
        delete[] SDP_block;
        SDP_block = NULL;
    }
    SDP_nBlock = other.SDP_nBlock;
    if ((SDP_nBlock > 0) && (SDP_block == NULL)) {
        rNewCheck();
        SDP_block = new DenseMatrix[SDP_nBlock];
        if (SDP_block == NULL) {
            rError("DenseLinearSpace:: memory exhausted");
        }
    }
    for (int l = 0; l < SDP_nBlock; ++l) {
        total_judge = SDP_block[l].copyFrom(other.SDP_block[l]);
    }
    if (total_judge == FAILURE) {
        rError("DenseLinearSpace:: copy miss");
    }

    // for LP
    if (other.LP_nBlock < 0) {
        rError("DenseLinearSpace:: LP_nBlock is negative");
    }
    if (LP_nBlock != other.LP_nBlock) {
        delete[] LP_block;
        LP_block = NULL;
    }
    LP_nBlock = other.LP_nBlock;
    if ((LP_nBlock > 0) && (LP_block == NULL)) {
        LP_block = new mpfr_class[LP_nBlock];
        if (LP_block == NULL) {
            rError("DenseLinearSpace:: memory exhausted");
        }
    }
    for (int l = 0; l < LP_nBlock; ++l) {
        LP_block[l] = other.LP_block[l];
    }

    return total_judge;
}

void DenseLinearSpace::setElement_SDP(int block, int i, int j, mpfr_class ele)
{

    // check range
    if (block >= SDP_nBlock) {
        rError("out of range in input data");
    }
    if ((i >= SDP_block[block].nRow) || (j >= SDP_block[block].nCol)) {
        rError("out of range in input data");
    }

    int nCol = SDP_block[block].nCol;
    SDP_block[block].de_ele[i + j * nCol] = ele;
    SDP_block[block].de_ele[j + i * nCol] = ele;
}

void DenseLinearSpace::setElement_LP(int block, mpfr_class ele)
{
    // check range
    if (block >= LP_nBlock) {
        rError("out of range in input data");
    }
    LP_block[block] = ele;
}

void DenseLinearSpace::setZero()
{
    // for SDP
    if (SDP_nBlock > 0 && SDP_block) {
        for (int l = 0; l < SDP_nBlock; ++l) {
            SDP_block[l].setZero();
        }
    }

    // for LP
    if (LP_nBlock > 0 && LP_block) {
        for (int l = 0; l < LP_nBlock; ++l) {
            LP_block[l] = 0.0;
        }
    }
}

void DenseLinearSpace::setIdentity(mpfr_class scalar)
{
    // for SDP
    if (SDP_nBlock > 0 && SDP_block) {
        for (int l = 0; l < SDP_nBlock; ++l) {
            SDP_block[l].setIdentity(scalar);
        }
    }
    // for LP
    if (LP_nBlock > 0 && LP_block) {
        for (int l = 0; l < LP_nBlock; ++l) {
            LP_block[l] = scalar;
        }
    }
}

} // namespace sdpa

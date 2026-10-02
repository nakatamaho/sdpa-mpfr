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

// printing presicion of such as vector
#define P_FORMAT "%+18.12Re"

#ifndef __sdpa_struct_h__
#define __sdpa_struct_h__

#include <sdpa_include.h>

namespace sdpa {

class Vector {
public:
    int nDim;
    mpfr_class* ele;

    Vector();
    ~Vector();

    void initialize(int nDim, mpfr_class value = 0.0);
    void initialize(mpfr_class value);
    void terminate();

    void setZero();
    void display(FILE* fpout, mpfr_class scalar);
    bool copyFrom(Vector& other);
};

class BlockVector {
public:
    int nBlock;
    int* blockStruct;

    Vector* ele;

    BlockVector();
    ~BlockVector();

    void initialize(int nBlock, int* blockStruct, mpfr_class value = 0.0);
    void terminate();
};

class SparseMatrix {
public:
    int nRow, nCol;

    enum Type { SPARSE, DENSE };
    Type type;

    int NonZeroNumber;
    // for memory
    int NonZeroCount;
    // currentry stored
    int NonZeroEffect;
    // use for calculation of F1,F2,F3

    // for Dense
    mpfr_class* de_ele;

    // for Sparse
    int* row_index;
    int* column_index;
    mpfr_class* sp_ele;

    SparseMatrix();
    ~SparseMatrix();

    void initialize(int nRow, int nCol, Type type, int NonZeroNumber);
    void terminate();

    void changeToDense(bool forceChange = false);
};

class DenseMatrix {
public:
    int nRow, nCol;

    enum Type { DENSE, COMPLETION };
    Type type;

    mpfr_class* de_ele;

    DenseMatrix();
    ~DenseMatrix();

    void initialize(int nRow, int nCol, Type type);
    void terminate();

    void display(FILE* fpout = stdout);
    bool copyFrom(DenseMatrix& other);

    void setZero();
    void setIdentity(mpfr_class scalar = 1.0);
};

class SparseLinearSpace {
public:
    int SDP_sp_nBlock;
    int SOCP_sp_nBlock;
    int LP_sp_nBlock;

    int* SDP_sp_index;
    int* SOCP_sp_index;
    int* LP_sp_index;

    SparseMatrix* SDP_sp_block;
    SparseMatrix* SOCP_sp_block;
    mpfr_class* LP_sp_block;

    SparseLinearSpace();
    ~SparseLinearSpace();

    // sparse form of block index      2008/02/27 kazuhide nakata
    void initialize(int SDP_sp_nBlock, int* SDP_sp_index, int* SDP_sp_blockStruct,
                    int* SDP_sp_NonZeroNumber, int SOCP_sp_nBlock, int* SOCP_sp_index,
                    int* SOCP_sp_blockStruct, int* SOCP_sp_NonZeroNumber, int LP_sp_nBlock,
                    int* LP_sp_index);
    void terminate();

    void changeToDense(bool forceChange = false);

    void setElement_SDP(int block, int nCol, int nRow, mpfr_class ele);
    void setElement_LP(int block, mpfr_class ele);
};

class DenseLinearSpace {
public:
    int SDP_nBlock;
    int SOCP_nBlock;
    int LP_nBlock;

    DenseMatrix* SDP_block;
    DenseMatrix* SOCP_block;
    mpfr_class* LP_block;

    DenseLinearSpace();
    ~DenseLinearSpace();
    void initialize(int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock, int* SOCP_blockStruct,
                    int LP_nBlock);
    void terminate();

    bool copyFrom(DenseLinearSpace& other);
    void setElement_SDP(int block, int nCol, int nRow, mpfr_class ele);
    void setElement_LP(int block, mpfr_class ele);
    void setZero();
    void setIdentity(mpfr_class scalar = 1.0);
};

} // namespace sdpa

#endif // __sdpa_struct_h__

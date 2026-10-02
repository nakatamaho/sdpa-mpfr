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

#include <sdpa_io.h>
#include <vector>
#include <algorithm>

namespace sdpa {

// 2008/02/27  kazuhide nakata

void IO::read(FILE* fpData, FILE* fpout, int& m, char* str)
{
    while (true) {
        fgets(str, lengthOfString, fpData);
        if (str[0] == '*' || str[0] == '"') {
            fprintf(fpout, "%s", str);
        } else {
            sscanf(str, "%d", &m);
            break;
        }
    }
}

void IO::read(FILE* fpData, int& nBlock)
{
    fscanf(fpData, "%d", &nBlock);
}

void IO::read(FILE* fpData, int nBlock, int* blockStruct)
{
    for (int l = 0; l < nBlock; ++l) {
        fscanf(fpData, "%*[^0-9+-]%d", &blockStruct[l]);
    }
}

void IO::read(FILE* fpData, Vector& b)
{
    for (int k = 0; k < b.nDim; ++k) {
        sdpa_fscan_real(fpData, &b.ele[k]);
    }
}

void IO::read(FILE* fpData, DenseLinearSpace& xMat, Vector& yVec, DenseLinearSpace& zMat,
              bool inputSparse)
{
    // read initial point
    int SDP_nBlock = xMat.SDP_nBlock;
    int SOCP_nBlock = xMat.SOCP_nBlock;
    int LP_nBlock = xMat.LP_nBlock;

    // yVec is opposite sign
    for (int k = 0; k < yVec.nDim; ++k) {
        mpfr_class tmp;
        sdpa_fscan_real(fpData, &tmp);
        yVec.ele[k] = -tmp;
    }

    if (inputSparse) {
        // sparse case , zMat , xMat in this order
        int i, j, l, target;
        mpfr_class value;
        while (true) {
            if (fscanf(fpData, "%*[^0-9+-]%d", &target) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &l) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &i) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &j) <= 0) {
                break;
            }
            if (sdpa_fscan_real(fpData, &value) <= 0) {
                break;
            }

            if (l <= SDP_nBlock) {
                // SDP part
                if (target == 1) {
                    zMat.setElement_SDP(l - 1, i - 1, j - 1, value);
                } else {
                    xMat.setElement_SDP(l - 1, i - 1, j - 1, value);
                }
            } else if (l <= SDP_nBlock + SOCP_nBlock) {
                // SOCP part
                rError("io:: current version does not support SOCP");
            } else {
                // LP part
                if (i != j) {
                    rError("io:: LP part  3rd elemtn != 4th elemnt");
                }
                if (target == 1) {
                    zMat.setElement_LP(i - 1, value);
                } else {
                    xMat.setElement_LP(i - 1, value);
                }
            }
        } // end of 'while (true)'
    } else {
        // dense case , zMat , xMat in this order
        // for SDP
        for (int l = 0; l < SDP_nBlock; ++l) {
            int size = zMat.SDP_block[l].nRow;
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    mpfr_class tmp;
                    sdpa_fscan_real(fpData, &tmp);
                    if (i <= j && tmp != 0.0) {
                        zMat.setElement_SDP(l, i, j, tmp);
                    }
                }
            }
        }
        // for SOCP
        for (int l = 0; l < SOCP_nBlock; ++l) {
            rError("io:: current version does not support SOCP");
        }
        // for LP
        for (int j = 0; j < LP_nBlock; ++j) {
            mpfr_class tmp;
            sdpa_fscan_real(fpData, &tmp);
            if (tmp != 0.0) {
                zMat.setElement_LP(j, tmp);
            }
        }

        // for SDP
        for (int l = 0; l < SDP_nBlock; ++l) {
            int size = xMat.SDP_block[l].nRow;
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    mpfr_class tmp;
                    sdpa_fscan_real(fpData, &tmp);
                    if (i <= j && tmp != 0.0) {
                        xMat.setElement_SDP(l, i, j, tmp);
                    }
                }
            }
        }
        // for SOCP
        for (int l = 0; l < SOCP_nBlock; ++l) {
            rError("io:: current version does not support SOCP");
        }
        // for LP
        for (int j = 0; j < LP_nBlock; ++j) {
            mpfr_class tmp;
            sdpa_fscan_real(fpData, &tmp);
            if (tmp != 0.0) {
                xMat.setElement_LP(j, tmp);
            }
        }
    } // end of 'if (inputSparse)'
}

// 2008/02/27 kazuhide nakata
// without LP_ANonZeroCount
void IO::read(FILE* fpData, int m, int SDP_nBlock, int* SDP_blockStruct, int SOCP_nBlock,
              int* SOCP_blockStruct, int LP_nBlock, int nBlock, int* blockStruct, int* blockType,
              int* blockNumber, InputData& inputData, bool isDataSparse)
{
    inputData.initialize_bVec(m);
    read(fpData, inputData.b);
    long position = ftell(fpData);

    // C,A must be accessed "mpfr_class".

    //   initialize block struct of C and A
    setBlockStruct(fpData, inputData, m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct,
                   LP_nBlock, nBlock, blockStruct, blockType, blockNumber, position, isDataSparse);

    setElement(fpData, inputData, m, SDP_nBlock, SDP_blockStruct, SOCP_nBlock, SOCP_blockStruct,
               LP_nBlock, nBlock, blockStruct, blockType, blockNumber, position, isDataSparse);
}

// 2008/02/27 kazuhide nakata
// without LP_ANonZeroCount
void IO::setBlockStruct(FILE* fpData, InputData& inputData, int m, int SDP_nBlock,
                        int* SDP_blockStruct, int SOCP_nBlock, int* SOCP_blockStruct, int LP_nBlock,
                        int nBlock, int* blockStruct, int* blockType, int* blockNumber,
                        long position, bool isDataSparse)
{
    // seed the positon of C in the fpData
    fseek(fpData, position, 0);

    vector<int>* SDP_index = NULL;
    SDP_index = new vector<int>[m + 1];
    vector<int>* SOCP_index = NULL;
    SOCP_index = new vector<int>[m + 1];
    vector<int>* LP_index = NULL;
    LP_index = new vector<int>[m + 1];

    // for SDP
    int SDP_sp_nBlock;
    int* SDP_sp_index = NULL;
    int* SDP_sp_blockStruct = NULL;
    int* SDP_sp_NonZeroNumber = NULL;
    SDP_sp_index = new int[SDP_nBlock];
    SDP_sp_blockStruct = new int[SDP_nBlock];
    SDP_sp_NonZeroNumber = new int[SDP_nBlock];
    // for SOCP
    int SOCP_sp_nBlock = 0;
    int* SOCP_sp_blockStruct = NULL;
    int* SOCP_sp_index = NULL;
    int* SOCP_sp_NonZeroNumber = NULL;
    // for LP
    int LP_sp_nBlock;
    int* LP_sp_index;
    LP_sp_index = new int[LP_nBlock];

    if (isDataSparse) {
        int i, j, k, l;
        mpfr_class value;
        while (true) {
            if (fscanf(fpData, "%*[^0-9+-]%d", &k) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &l) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &i) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &j) <= 0) {
                break;
            }
            if (sdpa_fscan_real(fpData, &value) <= 0) {
                break;
            }

            if (blockType[l - 1] == 1) { // SDP part
                int l2 = blockNumber[l - 1];
                SDP_index[k].push_back(l2);
            } else if (blockType[l - 1] == 2) { // SOCP part
                rError("io:: current version does not support SOCP");
            } else if (blockType[l - 1] == 3) { // LP part
                if (i != j) {
                    printf("invalid data file k:%d, l:%d, i:%d, j:%d, ", k, l, i, j);
                    mpfr_printf("value:%9.1Re\n", value.get_mpfr_t());
                    rError("IO::initializeLinearSpace");
                }
                int l2 = blockNumber[l - 1];
                LP_index[k].push_back(l2 + i - 1);
            } else {
                rError("io::read not valid blockType");
            }
        } // end of 'while (true)'

    } else { // isDataSparse == false

        // k==0
        for (int l2 = 0; l2 < nBlock; ++l2) {
            if (blockType[l2] == 1) { // SDP part
                int l = blockNumber[l2];
                int size = SDP_blockStruct[l];
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < size; ++j) {
                        mpfr_class tmp;
                        sdpa_fscan_real(fpData, &tmp);
                        if (i <= j && tmp != 0.0) {
                            SDP_index[0].push_back(l);
                        }
                    }
                }
            } else if (blockType[l2] == 2) { // SOCP part
                rError("io:: current version does not support SOCP");
            } else if (blockType[l2] == 3) { // LP part
                for (int j = 0; j < blockStruct[l2]; ++j) {
                    mpfr_class tmp;
                    sdpa_fscan_real(fpData, &tmp);
                    if (tmp != 0.0) {
                        LP_index[0].push_back(blockNumber[l2] + j);
                    }
                }
            } else {
                rError("io::read not valid blockType");
            }
        }

        for (int k = 0; k < m; ++k) {
            // k>0
            for (int l2 = 0; l2 < nBlock; ++l2) {
                if (blockType[l2] == 1) { // SDP part
                    int l = blockNumber[l2];
                    int size = SDP_blockStruct[l];
                    for (int i = 0; i < size; ++i) {
                        for (int j = 0; j < size; ++j) {
                            mpfr_class tmp;
                            sdpa_fscan_real(fpData, &tmp);
                            if (i <= j && tmp != 0.0) {
                                SDP_index[k + 1].push_back(l);
                            }
                        }
                    }
                } else if (blockType[l2] == 2) { // SOCP part
                    rError("io:: current version does not support SOCP");
                } else if (blockType[l2] == 3) { // LP part
                    for (int j = 0; j < blockStruct[l2]; ++j) {
                        mpfr_class tmp;
                        sdpa_fscan_real(fpData, &tmp);
                        if (tmp != 0.0) {
                            LP_index[k + 1].push_back(blockNumber[l2] + j);
                        }
                    }
                } else {
                    rError("io::read not valid blockType");
                }
            }
        }

    } // end of 'if (isDataSparse)'

    inputData.A = new SparseLinearSpace[m];
    for (int k = 0; k < m + 1; k++) {
        sort(SDP_index[k].begin(), SDP_index[k].end());
        SDP_sp_nBlock = 0;
        int previous_index = -1;
        int index;
        for (int i = 0; i < SDP_index[k].size(); i++) {
            index = SDP_index[k][i];
            if (previous_index != index) {
                SDP_sp_index[SDP_sp_nBlock] = index;
                SDP_sp_blockStruct[SDP_sp_nBlock] = SDP_blockStruct[index];
                SDP_sp_NonZeroNumber[SDP_sp_nBlock] = 1;
                previous_index = index;
                SDP_sp_nBlock++;
            } else {
                SDP_sp_NonZeroNumber[SDP_sp_nBlock - 1]++;
            }
        }
        sort(LP_index[k].begin(), LP_index[k].end());
        LP_sp_nBlock = 0;
        previous_index = -1;
        for (int i = 0; i < LP_index[k].size(); i++) {
            index = LP_index[k][i];
            if (previous_index != index) {
                LP_sp_index[LP_sp_nBlock] = index;
                previous_index = index;
                LP_sp_nBlock++;
            }
        }

        if (k == 0) {
            inputData.C.initialize(SDP_sp_nBlock, SDP_sp_index, SDP_sp_blockStruct,
                                   SDP_sp_NonZeroNumber, SOCP_sp_nBlock, SOCP_sp_blockStruct,
                                   SOCP_sp_index, SOCP_sp_NonZeroNumber, LP_sp_nBlock, LP_sp_index);
        } else {
            inputData.A[k - 1].initialize(SDP_sp_nBlock, SDP_sp_index, SDP_sp_blockStruct,
                                          SDP_sp_NonZeroNumber, SOCP_sp_nBlock, SOCP_sp_blockStruct,
                                          SOCP_sp_index, SOCP_sp_NonZeroNumber, LP_sp_nBlock,
                                          LP_sp_index);
        }
    }

    delete[] SDP_index;
    SDP_index = NULL;
    delete[] SOCP_index;
    SOCP_index = NULL;
    delete[] LP_index;
    LP_index = NULL;

    delete[] SDP_sp_index;
    SDP_sp_index = NULL;
    delete[] SDP_sp_blockStruct;
    SDP_sp_blockStruct = NULL;
    delete[] SDP_sp_NonZeroNumber;
    SDP_sp_NonZeroNumber = NULL;
    delete[] LP_sp_index;
    LP_sp_index = NULL;
}

// 2008/02/27 kazuhide nakata
// without LP_ANonZeroCount
void IO::setElement(FILE* fpData, InputData& inputData, int m, int SDP_nBlock, int* SDP_blockStruct,
                    int SOCP_nBlock, int* SOCP_blockStruct, int LP_nBlock, int nBlock,
                    int* blockStruct, int* blockType, int* blockNumber, long position,
                    bool isDataSparse)
{
    // in Sparse, read C,A[k]

    // seed the positon of C in the fpData
    fseek(fpData, position, 0);

    if (isDataSparse) {
        int i, j, k, l;
        mpfr_class value;
        while (true) {
            if (fscanf(fpData, "%*[^0-9+-]%d", &k) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &l) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &i) <= 0) {
                break;
            }
            if (fscanf(fpData, "%*[^0-9+-]%d", &j) <= 0) {
                break;
            }
            if (sdpa_fscan_real(fpData, &value) <= 0) {
                break;
            }

            if (blockType[l - 1] == 1) { // SDP part
                int l2 = blockNumber[l - 1];
                if (k == 0) {
                    inputData.C.setElement_SDP(l2, i - 1, j - 1, -value);
                } else {
                    inputData.A[k - 1].setElement_SDP(l2, i - 1, j - 1, value);
                }
            } else if (blockType[l - 1] == 2) { // SOCP part
                rError("io:: current version does not support SOCP");
            } else if (blockType[l - 1] == 3) { // LP part
                if (i != j) {
                    rError("io:: LP part  3rd elemtn != 4th elemnt");
                }
                if (k == 0) {
                    inputData.C.setElement_LP(blockNumber[l - 1] + i - 1, -value);
                } else {
                    inputData.A[k - 1].setElement_LP(blockNumber[l - 1] + i - 1, value);
                }
            } else {
                rError("io::read not valid blockType");
            }
        }
    } else { // dense

        // k==0
        for (int l2 = 0; l2 < nBlock; ++l2) {
            if (blockType[l2] == 1) { // SDP part
                int l = blockNumber[l2];
                int size = SDP_blockStruct[l];
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < size; ++j) {
                        mpfr_class tmp;
                        sdpa_fscan_real(fpData, &tmp);
                        if (i <= j && tmp != 0.0) {
                            inputData.C.setElement_SDP(l, i, j, -tmp);
                        }
                    }
                }
            } else if (blockType[l2] == 2) { // SOCP part
                rError("io:: current version does not support SOCP");
            } else if (blockType[l2] == 3) { // LP part
                for (int j = 0; j < blockStruct[l2]; ++j) {
                    mpfr_class tmp;
                    sdpa_fscan_real(fpData, &tmp);
                    if (tmp != 0.0) {
                        inputData.C.setElement_LP(blockNumber[l2] + j, -tmp);
                    }
                }
            } else {
                rError("io::read not valid blockType");
            }
        }

        // k > 0
        for (int k = 0; k < m; ++k) {

            for (int l2 = 0; l2 < nBlock; ++l2) {
                if (blockType[l2] == 1) { // SDP part
                    int l = blockNumber[l2];
                    int size = SDP_blockStruct[l];
                    for (int i = 0; i < size; ++i) {
                        for (int j = 0; j < size; ++j) {
                            mpfr_class tmp;
                            sdpa_fscan_real(fpData, &tmp);
                            if (i <= j && tmp != 0.0) {
                                inputData.A[k].setElement_SDP(l, i, j, tmp);
                            }
                        }
                    }
                } else if (blockType[l2] == 2) { // SOCP part
                    rError("io:: current version does not support SOCP");
                } else if (blockType[l2] == 3) { // LP part
                    for (int j = 0; j < blockStruct[l2]; ++j) {
                        mpfr_class tmp;
                        sdpa_fscan_real(fpData, &tmp);
                        if (tmp != 0.0) {
                            inputData.A[k].setElement_LP(blockNumber[l2] + j, tmp);
                        }
                    }
                } else {
                    rError("io::read not valid blockType");
                }
            }

        } // for k

    } // end of 'if (isDataSparse)'
}

void IO::printHeader(FILE* fpout, FILE* Display)
{
    if (fpout) {
        fprintf(fpout, "   mu      thetaP  thetaD  objP      objD "
                       "     alphaP  alphaD  beta \n");
    }
    if (Display) {
        fprintf(Display, "   mu      thetaP  thetaD  objP      objD "
                         "     alphaP  alphaD  beta \n");
    }
}

void IO::printOneIteration(int pIteration, AverageComplementarity& mu,
                           RatioInitResCurrentRes& theta, SolveInfo& solveInfo, StepLength& alpha,
                           DirectionParameter& beta, FILE* fpout, FILE* Display)
{
    if (Display) {
        mpfr_class mtmp1 = -solveInfo.objValDual;
        mpfr_class mtmp2 = -solveInfo.objValPrimal;
        mpfr_fprintf(Display,
                     "%2d %4.1Re %4.1Re %4.1Re %+7.2Re %+7.2Re"
                     " %4.1Re %4.1Re %4.2Re\n",
                     pIteration, mu.current.get_mpfr_t(), theta.dual.get_mpfr_t(),
                     theta.primal.get_mpfr_t(), mtmp1.get_mpfr_t(), mtmp2.get_mpfr_t(),
                     alpha.dual.get_mpfr_t(), alpha.primal.get_mpfr_t(), beta.value.get_mpfr_t());
    }
    if (fpout) {
        mpfr_class mtmp1 = -solveInfo.objValDual;
        mpfr_class mtmp2 = -solveInfo.objValPrimal;
        mpfr_fprintf(fpout,
                     "%2d %4.1Re %4.1Re %4.1Re %+7.2Re %+7.2Re"
                     " %4.1Re %4.1Re %4.2Re\n",
                     pIteration, mu.current.get_mpfr_t(), theta.dual.get_mpfr_t(),
                     theta.primal.get_mpfr_t(), mtmp1.get_mpfr_t(), mtmp2.get_mpfr_t(),
                     alpha.dual.get_mpfr_t(), alpha.primal.get_mpfr_t(), beta.value.get_mpfr_t());
    }
}

void IO::printLastInfo(int pIteration, AverageComplementarity& mu, RatioInitResCurrentRes& theta,
                       SolveInfo& solveInfo, StepLength& alpha, DirectionParameter& beta,
                       Residuals& currentRes, Phase& phase, Solutions& currentPt, double cputime,
                       int nBlock, int* blockStruct, int* blockType, int* blockNumber,
                       InputData& inputData, WorkVariables& work, ComputeTime& com,
                       Parameter& param, FILE* fpout, FILE* Display, bool printTime)
{
    int nDim = currentPt.nDim;

    printOneIteration(pIteration, mu, theta, solveInfo, alpha, beta, fpout, Display);

    mpfr_class mean = (abs(solveInfo.objValPrimal) + abs(solveInfo.objValDual)) / 2.0;
    mpfr_class PDgap = abs(solveInfo.objValPrimal - solveInfo.objValDual);
    mpfr_class relgap;
    if (mean < 1.0) {
        relgap = PDgap;
    } else {
        relgap = PDgap / mean;
    }
    mpfr_class gap = mu.current * nDim;
    mpfr_class digits = 0.0;
    if (PDgap != 0.0 && mean != 0.0) {
        digits = -log10(abs(PDgap / mean));
    }

    if (Display) {
        fprintf(Display, "\n");
        phase.display(Display);
        mpfr_fprintf(Display, "   Iteration = %d\n", pIteration);
        mpfr_fprintf(Display, "          mu = %4.16Re\n", mu.current.get_mpfr_t());
        mpfr_fprintf(Display, "relative gap = %4.16Re\n", relgap.get_mpfr_t());
        mpfr_fprintf(Display, "         gap = %4.16Re\n", gap.get_mpfr_t());
        mpfr_fprintf(Display, "      digits = %4.16Re\n", digits.get_mpfr_t());

        mpfr_class mtmp1 = -solveInfo.objValDual;
        mpfr_class mtmp2 = -solveInfo.objValPrimal;
        mpfr_fprintf(Display, "objValPrimal = %10.16Re\n", mtmp1.get_mpfr_t());
        mpfr_fprintf(Display, "objValDual   = %10.16Re\n", mtmp2.get_mpfr_t());
        mpfr_fprintf(Display, "p.feas.error = %10.16Re\n", currentRes.normDualMat.get_mpfr_t());
        mpfr_fprintf(Display, "d.feas.error = %10.16Re\n", currentRes.normPrimalVec.get_mpfr_t());
        mpfr_fprintf(Display, "relative eps = %10.16Re\n", Rlamch_mpfr("E").get_mpfr_t());
        if (printTime == true) {
            fprintf(Display, "total time   = %.3f\n", cputime);
        }
    }
    if (fpout) {
        fprintf(fpout, "\n");
        phase.display(fpout);
        fprintf(fpout, "   Iteration = %d\n", pIteration);
        mpfr_fprintf(fpout, "          mu = %4.16Re\n", mu.current.get_mpfr_t());
        mpfr_fprintf(fpout, "relative gap = %4.16Re\n", relgap.get_mpfr_t());
        mpfr_fprintf(fpout, "         gap = %4.16Re\n", gap.get_mpfr_t());
        mpfr_fprintf(fpout, "      digits = %4.16Re\n", digits.get_mpfr_t());

        mpfr_class mtmp1 = -solveInfo.objValDual;
        mpfr_class mtmp2 = -solveInfo.objValPrimal;
        mpfr_fprintf(fpout, "objValPrimal = %10.16Re\n", mtmp1.get_mpfr_t());
        mpfr_fprintf(fpout, "objValDual   = %10.16Re\n", mtmp2.get_mpfr_t());
        mpfr_fprintf(fpout, "p.feas.error = %10.16Re\n", currentRes.normDualMat.get_mpfr_t());
        mpfr_fprintf(fpout, "d.feas.error = %10.16Re\n", currentRes.normPrimalVec.get_mpfr_t());
        mpfr_fprintf(fpout, "relative eps = %10.16Re\n", Rlamch_mpfr("E").get_mpfr_t());
        fprintf(fpout, "total time   = %.3f\n", cputime);

        fprintf(fpout, "\n\nParameters are\n");
        param.display(fpout);
        com.display(fpout);

        fprintf(fpout, "xVec = \n");
        currentPt.yVec.display(fpout, -1.0);
        fprintf(fpout, "xMat = \n");
        displayDenseLinarSpaceLast(currentPt.zMat, nBlock, blockStruct, blockType, blockNumber,
                                   fpout);
        fprintf(fpout, "yMat = \n");
        displayDenseLinarSpaceLast(currentPt.xMat, nBlock, blockStruct, blockType, blockNumber,
                                   fpout);
    }
}

void IO::displayDenseLinarSpaceLast(DenseLinearSpace& aMat, int nBlock, int* blockStruct,
                                    int* blockType, int* blockNumber, FILE* fpout)
{
    if (fpout == NULL) {
        return;
    }

    fprintf(fpout, "{\n");
    for (int i = 0; i < nBlock; i++) {
        if (blockType[i] == 1) {
            int l = blockNumber[i];
            aMat.SDP_block[l].display(fpout);
        } else if (blockType[i] == 2) {
            rError("io:: current version does not support SOCP");
        } else if (blockType[i] == 3) {
            fprintf(fpout, "{");
            for (int l = 0; l < blockStruct[i] - 1; ++l) {
                mpfr_fprintf(fpout, P_FORMAT ",", aMat.LP_block[blockNumber[i] + l].get_mpfr_t());
            }
            if (blockStruct[i] > 0) {
                mpfr_fprintf(fpout, P_FORMAT "}\n",
                             aMat.LP_block[blockNumber[i] + blockStruct[i] - 1].get_mpfr_t());
            } else {
                fprintf(fpout, "  }\n");
            }
        } else {
            rError("io::displayDenseLinarSpaceLast not valid blockType");
        }
    }
    fprintf(fpout, "}\n");
}

} // namespace sdpa

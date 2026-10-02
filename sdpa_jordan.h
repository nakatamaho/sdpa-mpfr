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

#ifndef __sdpa_jordan_h__
#define __sdpa_jordan_h__

#include <sdpa_linear.h>

namespace sdpa {
class WorkVariables;
}

namespace sdpa {

class Jal {
public:
    // calculate the minimum eigen value of lMat*xMat*(lMat^T)
    // by Lanczos methods.
    // lMat is lower triangular¡¢xMat is symmetric
    // block size > 20   : Lanczos method
    // block size <= 20  : QR method
    static mpfr_class getMinEigen(DenseLinearSpace& lMat, DenseLinearSpace& xMat,
                                  WorkVariables& work);

    static bool getInvChol(DenseLinearSpace& invCholMat, DenseLinearSpace& aMat,
                           DenseLinearSpace& workMat);

    static bool getInvCholAndInv(DenseLinearSpace& invCholMat, DenseLinearSpace& inverseMat,
                                 DenseLinearSpace& aMat, DenseLinearSpace& workMat);

    static bool multiply(DenseLinearSpace& retMat, DenseLinearSpace& aMat, DenseLinearSpace& bMat,
                         mpfr_class* scalar = NULL);

    //  retMat = A * B * C
    static bool ns_jordan_triple_product(DenseLinearSpace& retMat, DenseLinearSpace& aMat,
                                         DenseLinearSpace& bMat, DenseLinearSpace& cMat,
                                         DenseLinearSpace& work);
};

} // namespace sdpa

#endif // __sdpa_jordan_h__

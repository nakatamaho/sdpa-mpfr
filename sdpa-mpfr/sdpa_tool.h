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
/*--------------------------------------------------
  rsdpa_tool.h
  $Id: rsdpa_tool.h,v 1.2 2004/09/01 06:34:12 makoto Exp $
--------------------------------------------------*/

#ifndef __sdpa_tool_h__
#define __sdpa_tool_h__

#include <sdpa_right.h>

#include <iostream>
#include <string>

#include <cstdio>
#include <mpfrxx_mkII.h>

using mpfrxx::mpfr_class;

namespace sdpa {

#define rMessage(message) cout << message << " :: line " << __LINE__ << " in " << __FILE__ << endl

#define rError(message)                                                                            \
    cout << message << " :: line " << __LINE__ << " in " << __FILE__ << endl;                      \
    exit(false)

#define rNewCheck() ;

// These are constant. Do NOT change
extern int IONE;         // =  1;
extern mpfr_class MZERO; // =  0.0;
extern mpfr_class MONE;  // =  1.0;
extern mpfr_class MMONE; // = -1.0;

// Sets the MPFR default precision (in bits) used by subsequently
// constructed mpfr_class objects, and re-initializes MZERO, MONE, MMONE
// at that precision.
void setDefaultPrecision(int precision);

// Replacement for gmp_fscanf(fp, "%*[^0-9+-]%Fe", value).
// Skips characters other than [0-9+-.], then reads one decimal
// floating-point number [+-]digits[.digits][(e|E)[+-]digits] into *value
// at its current precision.  The character following the number is left
// in the stream.
// Returns 1 on success, 0 on a malformed number, EOF at end of file.
int sdpa_fscan_real(FILE* fp, mpfr_class* value);

class Time {
public:
    static double rGetUseTime();
};

#define TimeStart(START__)                                                                         \
    static double START__;                                                                         \
    START__ = Time::rGetUseTime()
#define TimeEnd(END__)                                                                             \
    static double END__;                                                                           \
    END__ = Time::rGetUseTime()
#define TimeCal(START__, END__) (END__ - START__)

} // namespace sdpa

#endif // __sdpa_tool_h__

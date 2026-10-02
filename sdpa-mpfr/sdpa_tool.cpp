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
/*-----------------------------------------
  rsdpa_tool.cpp
  $Id: rsdpa_tool.cpp,v 1.2 2004/09/01 06:34:12 makoto Exp $
-----------------------------------------*/

#include <sdpa_tool.h>
#include <sys/times.h>
#include <time.h>

#include <unistd.h>
#ifndef CLK_TCK
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif

using std::cout;
using std::endl;

namespace sdpa {

// These are constant.
// Do Not Change .
int IONE = 1;
mpfr_class MZERO = 0.0;
mpfr_class MONE = 1.0;
mpfr_class MMONE = -1.0;

void setDefaultPrecision(int precision)
{
    if (precision < (int)MPFR_PREC_MIN) {
        rError("precision must be at least " << (int)MPFR_PREC_MIN << " bits");
    }
    mpfrxx::set_default_precision_bits((mpfr_prec_t)precision);
    MZERO.set_prec(precision);
    MZERO = 0.0;
    MONE.set_prec(precision);
    MONE = 1.0;
    MMONE.set_prec(precision);
    MMONE = -1.0;
}

static inline bool isDigit(int c)
{
    return '0' <= c && c <= '9';
}

// Appends a run of decimal digits to token; returns the number of digits.
static int scanDigits(FILE* fp, std::string& token, int& c)
{
    int n = 0;
    while (isDigit(c)) {
        token += (char)c;
        c = getc(fp);
        ++n;
    }
    return n;
}

int sdpa_fscan_real(FILE* fp, mpfr_class* value)
{
    int c;
    // skip separators such as spaces, commas, braces and parentheses
    while ((c = getc(fp)) != EOF) {
        if (isDigit(c) || c == '+' || c == '-' || c == '.') {
            break;
        }
    }
    if (c == EOF) {
        return EOF;
    }
    // [+-] digits [. digits] [(e|E) [+-] digits]
    std::string token;
    if (c == '+' || c == '-') {
        token += (char)c;
        c = getc(fp);
    }
    int nDigits = scanDigits(fp, token, c);
    if (c == '.') {
        token += (char)c;
        c = getc(fp);
        nDigits += scanDigits(fp, token, c);
    }
    bool valid = nDigits > 0;
    if (valid && (c == 'e' || c == 'E')) {
        token += (char)c;
        c = getc(fp);
        if (c == '+' || c == '-') {
            token += (char)c;
            c = getc(fp);
        }
        valid = scanDigits(fp, token, c) > 0;
    }
    // the first character after the number is left in the stream
    if (c != EOF) {
        ungetc(c, fp);
    }
    if (!valid) {
        return 0;
    }
    mpfr_set_str(value->get_mpfr_t(), token.c_str(), 10, MPFR_RNDN);
    return 1;
}

double Time::rGetUseTime()
{
    struct tms TIME;
    times(&TIME);
    return (double)TIME.tms_utime / (double)CLK_TCK;
}

} // namespace sdpa

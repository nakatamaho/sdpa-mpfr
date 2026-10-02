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
#include <sys/time.h>
#include <time.h>

#include <unistd.h>
#ifndef CLK_TCK
#define  CLK_TCK  sysconf(_SC_CLK_TCK)
#endif

using std::cout;
using std::endl;

namespace sdpa {

// These are constant.
// Do Not Change .
int IZERO =  0;
int IONE  =  1;
int IMONE = -1;
mpfr_class MZERO =  0.0;
mpfr_class MONE  =  1.0;
mpfr_class MMONE = -1.0;

void setDefaultPrecision(int precision)
{
  if (precision < (int)MPFR_PREC_MIN) {
    rError("precision must be at least " << (int)MPFR_PREC_MIN << " bits");
  }
  mpfrxx::set_default_precision_bits((mpfr_prec_t)precision);
  MZERO.set_prec(precision); MZERO =  0.0;
  MONE .set_prec(precision); MONE  =  1.0;
  MMONE.set_prec(precision); MMONE = -1.0;
}

int sdpa_fscan_real(FILE* fp, mpfr_class* value)
{
  int c;
  // skip separators such as spaces, commas, braces and parentheses
  while ((c = getc(fp)) != EOF) {
    if (('0' <= c && c <= '9') || c == '+' || c == '-' || c == '.') {
      break;
    }
  }
  if (c == EOF) {
    return EOF;
  }
  std::string token;
  // accept the characters of a decimal floating-point literal
  while (c != EOF
	 && (('0' <= c && c <= '9') || c == '+' || c == '-'
	     || c == '.' || c == 'e' || c == 'E')) {
    token += (char)c;
    c = getc(fp);
  }
  if (c != EOF) {
    ungetc(c, fp);
  }
  char* end = NULL;
  mpfr_strtofr(value->get_mpfr_t(), token.c_str(), &end, 10, MPFR_RNDN);
  if (end == token.c_str()) {
    return 0;
  }
  return 1;
}

double Time::rGetUseTime()
{
  struct tms TIME;
  times(&TIME);
  return (double)TIME.tms_utime/(double)CLK_TCK; 
}

void Time::rSetTimeVal(struct timeval& targetVal)
{
  static struct timezone tz;
  gettimeofday(&targetVal,&tz);
}

double Time::rGetRealTime(const struct timeval& start,
			   const struct timeval& end)
{
  const long int second = end.tv_sec - start.tv_sec;
  const long int usecond = end.tv_usec - start.tv_usec;
  return ((double)second) + ((double)usecond)*(1.0e-6);
}

}


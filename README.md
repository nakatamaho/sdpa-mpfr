# sdpa-mpfr

SDPA-MPFR is a multiple precision version of SDPA (SemiDefinite Programming
Algorithm) based on [MPFR](https://www.mpfr.org/).
It is derived from SDPA-GMP 7.1.2; the GMP `mpf_class` arithmetic has been
replaced by `mpfrxx::mpfr_class` of
[gmpfrxx_mkII](https://github.com/nakatamaho/gmpfrxx_mkII), a header-only
C++17 wrapper of GMP/MPFR.

## Requirements

- C++17 compiler
- GMP
- MPFR (built with thread-local storage support, which is the default)
- gmpfrxx_mkII (header only; `mpfrxx_mkII.h` must be reachable)

## Build

Install gmpfrxx_mkII (its generated headers are produced by CMake):

```sh
git clone https://github.com/nakatamaho/gmpfrxx_mkII
cmake -S gmpfrxx_mkII -B gmpfrxx_mkII/build -DGMPFRXX_MKII_COMPONENTS=GMP,MPFR \
      -DGMPFRXX_MKII_BUILD_EXAMPLES=OFF -DGMPFRXX_MKII_BUILD_BENCHMARKS=OFF
cmake --build gmpfrxx_mkII/build
cmake --install gmpfrxx_mkII/build --prefix $HOME/opt/gmpfrxx_mkII
```

Then build SDPA-MPFR:

```sh
./configure --with-gmpfrxx-mkii-includedir=$HOME/opt/gmpfrxx_mkII/include
make
```

Other options: `--with-gmp-includedir`, `--with-gmp-libdir`,
`--with-mpfr-includedir`, `--with-mpfr-libdir`, `--with-system-spooles`,
`--enable-shared`, `--disable-mpfr-fast`.

By default the gmpfrxx_mkII fast paths `GMPFRXX_MKII_FAST_FIXED_PREC` and
`GMPFRXX_MKII_FAST_STABLE_RND` are enabled.  Their contracts hold for
SDPA-MPFR: every `mpfr_class` is kept at the precision given in
`param.sdpa`, and the MPFR rounding mode is never changed.  Together with
the `X += a * b` form used in the MBLAS kernels they let gmpfrxx_mkII reuse
scratch storage instead of allocating a temporary per operation (about
15-20% faster on SDPLIB theta1, control2 and mcp100).  Use
`--disable-mpfr-fast` to build without them.

## Usage

```sh
./sdpa_mpfr example1.dat example1.result
./sdpa_mpfr -ds example1.dat-s -o example1.result -p param.sdpa
```

## Precision

The last line of `param.sdpa` (default 200) gives the MPFR precision in bits.
MPFR uses exactly the requested number of bits, whereas GMP `mpf_t` rounded
the precision up to whole limbs (and added a guard limb).  With the same
`precision` value SDPA-MPFR therefore works with fewer bits than SDPA-GMP;
increase the value if you need the same accuracy.

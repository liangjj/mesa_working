*deck i1mach
      integer function i1mach (i)
c***begin prologue  i1mach
c***purpose  return integer machine dependent constants.
c***library   slatec
c***category  r1
c***type      integer (i1mach-i)
c***keywords  machine constants
c***author  fox, p. a., (bell labs)
c           hall, a. d., (bell labs)
c           schryer, n. l., (bell labs)
c***description
c
c   i1mach can be used to obtain machine-dependent parameters for the
c   local machine environment.  it is a function subprogram with one
c   (input) argument and can be referenced as follows:
c
c        k = i1mach(i)
c
c   where i=1,...,16.  the (output) value of k above is determined by
c   the (input) value of i.  the results for various values of i are
c   discussed below.
c
c   i/o unit numbers:
c     i1mach( 1) = the standard input unit.
c     i1mach( 2) = the standard output unit.
c     i1mach( 3) = the standard punch unit.
c     i1mach( 4) = the standard error message unit.
c
c   words:
c     i1mach( 5) = the number of bits per integer storage unit.
c     i1mach( 6) = the number of characters per integer storage unit.
c
c   integers:
c     assume integers are represented in the s-digit, base-a form
c
c                sign ( x(s-1)*a**(s-1) + ... + x(1)*a + x(0) )
c
c                where 0 .le. x(i) .lt. a for i=0,...,s-1.
c     i1mach( 7) = a, the base.
c     i1mach( 8) = s, the number of base-a digits.
c     i1mach( 9) = a**s - 1, the largest magnitude.
c
c   floating-point numbers:
c     assume floating-point numbers are represented in the t-digit,
c     base-b form
c                sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c                where 0 .le. x(i) .lt. b for i=1,...,t,
c                0 .lt. x(1), and emin .le. e .le. emax.
c     i1mach(10) = b, the base.
c
c   single-precision:
c     i1mach(11) = t, the number of base-b digits.
c     i1mach(12) = emin, the smallest exponent e.
c     i1mach(13) = emax, the largest exponent e.
c
c   double-precision:
c     i1mach(14) = t, the number of base-b digits.
c     i1mach(15) = emin, the smallest exponent e.
c     i1mach(16) = emax, the largest exponent e.
c
c   to alter this function for a particular environment, the desired
c   set of data statements should be activated by removing the c from
c   column 1.  also, the values of i1mach(1) - i1mach(4) should be
c   checked for consistency with the local operating system.
c
c***references  p. a. fox, a. d. hall and n. l. schryer, framework for
c                 a portable library, acm transactions on mathematical
c                 software 4, 2 (june 1978), pp. 177-188.
c***routines called  (none)
c***revision history  (yymmdd)
c   750101  date written
c   891012  added vax g-floating constants.  (wrb)
c   891012  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900618  added dec risc constants.  (wrb)
c   900723  added ibm rs 6000 constants.  (wrb)
c   901009  correct i1mach(7) for ibm mainframes. should be 2 not 16.
c           (rwc)
c   910710  added hp 730 constants.  (smr)
c   911114  added convex ieee constants.  (wrb)
c   920121  added sun -r8 compiler option constants.  (wrb)
c   920229  added touchstone delta i860 constants.  (wrb)
c   920501  reformatted the references section.  (wrb)
c   920625  added convex -p8 and -pd8 compiler option constants.
c           (bks, wrb)
c   930201  added dec alpha and sgi constants.  (rwc and wrb)
c   930618  corrected i1mach(5) for convex -p8 and -pd8 compiler
c           options.  (dwl, rwc and wrb).
c***end prologue  i1mach
c
      integer imach(16),output
      save imach
      equivalence (imach(4),output)
c
c     machine constants for the amiga
c     absoft compiler
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1022 /
c     data imach(16) /       1023 /
c
c     machine constants for the apollo
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        129 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1025 /
c
c     machine constants for the burroughs 1700 system
c
c     data imach( 1) /          7 /
c     data imach( 2) /          2 /
c     data imach( 3) /          2 /
c     data imach( 4) /          2 /
c     data imach( 5) /         36 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         33 /
c     data imach( 9) / z1ffffffff /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -256 /
c     data imach(13) /        255 /
c     data imach(14) /         60 /
c     data imach(15) /       -256 /
c     data imach(16) /        255 /
c
c     machine constants for the burroughs 5700 system
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         48 /
c     data imach( 6) /          6 /
c     data imach( 7) /          2 /
c     data imach( 8) /         39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /          8 /
c     data imach(11) /         13 /
c     data imach(12) /        -50 /
c     data imach(13) /         76 /
c     data imach(14) /         26 /
c     data imach(15) /        -50 /
c     data imach(16) /         76 /
c
c     machine constants for the burroughs 6700/7700 systems
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         48 /
c     data imach( 6) /          6 /
c     data imach( 7) /          2 /
c     data imach( 8) /         39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /          8 /
c     data imach(11) /         13 /
c     data imach(12) /        -50 /
c     data imach(13) /         76 /
c     data imach(14) /         26 /
c     data imach(15) /     -32754 /
c     data imach(16) /      32780 /
c
c     machine constants for the cdc 170/180 series using nos/ve
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         64 /
c     data imach( 6) /          8 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 9223372036854775807 /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /      -4095 /
c     data imach(13) /       4094 /
c     data imach(14) /         94 /
c     data imach(15) /      -4095 /
c     data imach(16) /       4094 /
c
c     machine constants for the cdc 6000/7000 series
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /    6loutput/
c     data imach( 5) /         60 /
c     data imach( 6) /         10 /
c     data imach( 7) /          2 /
c     data imach( 8) /         48 /
c     data imach( 9) / 00007777777777777777b /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /       -929 /
c     data imach(13) /       1070 /
c     data imach(14) /         94 /
c     data imach(15) /       -929 /
c     data imach(16) /       1069 /
c
c     machine constants for the celerity c1260
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          0 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / z'7fffffff' /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1022 /
c     data imach(16) /       1023 /
c
c     machine constants for the convex
c     using the -fn compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the convex
c     using the -fi compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the convex
c     using the -p8 compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         64 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 9223372036854775807 /
c     data imach(10) /          2 /
c     data imach(11) /         53 /
c     data imach(12) /      -1023 /
c     data imach(13) /       1023 /
c     data imach(14) /        113 /
c     data imach(15) /     -16383 /
c     data imach(16) /      16383 /
c
c     machine constants for the convex
c     using the -pd8 compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         64 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 9223372036854775807 /
c     data imach(10) /          2 /
c     data imach(11) /         53 /
c     data imach(12) /      -1023 /
c     data imach(13) /       1023 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the cray
c     using the 46 bit integer compiler option
c
c     data imach( 1) /        100 /
c     data imach( 2) /        101 /
c     data imach( 3) /        102 /
c     data imach( 4) /        101 /
c     data imach( 5) /         64 /
c     data imach( 6) /          8 /
c     data imach( 7) /          2 /
c     data imach( 8) /         46 /
c     data imach( 9) / 1777777777777777b /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /      -8189 /
c     data imach(13) /       8190 /
c     data imach(14) /         94 /
c     data imach(15) /      -8099 /
c     data imach(16) /       8190 /
c
c     machine constants for the cray
c     using the 64 bit integer compiler option
c
c     data imach( 1) /        100 /
c     data imach( 2) /        101 /
c     data imach( 3) /        102 /
c     data imach( 4) /        101 /
c     data imach( 5) /         64 /
c     data imach( 6) /          8 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 777777777777777777777b /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /      -8189 /
c     data imach(13) /       8190 /
c     data imach(14) /         94 /
c     data imach(15) /      -8099 /
c     data imach(16) /       8190 /
c
c     machine constants for the data general eclipse s/200
c
c     data imach( 1) /         11 /
c     data imach( 2) /         12 /
c     data imach( 3) /          8 /
c     data imach( 4) /         10 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /         16 /
c     data imach(11) /          6 /
c     data imach(12) /        -64 /
c     data imach(13) /         63 /
c     data imach(14) /         14 /
c     data imach(15) /        -64 /
c     data imach(16) /         63 /
c
c     machine constants for the dec alpha
c     using g_float
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the dec alpha
c     using ieee_float
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the dec risc
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the dec vax
c     using d_floating
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the dec vax
c     using g_floating
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the elxsi 6400
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         32 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1022 /
c     data imach(16) /       1023 /
c
c     machine constants for the harris 220
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          0 /
c     data imach( 4) /          6 /
c     data imach( 5) /         24 /
c     data imach( 6) /          3 /
c     data imach( 7) /          2 /
c     data imach( 8) /         23 /
c     data imach( 9) /    8388607 /
c     data imach(10) /          2 /
c     data imach(11) /         23 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         38 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the honeywell 600/6000 series
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /         43 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          6 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         63 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the hp 730
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the hp 2100
c     3 word double precision option with ftn4
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          4 /
c     data imach( 4) /          1 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         23 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         39 /
c     data imach(15) /       -128 /
c     data imach(16) /        127 /
c
c     machine constants for the hp 2100
c     4 word double precision option with ftn4
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          4 /
c     data imach( 4) /          1 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         23 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         55 /
c     data imach(15) /       -128 /
c     data imach(16) /        127 /
c
c     machine constants for the hp 9000
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          7 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         32 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1015 /
c     data imach(16) /       1017 /
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9, the sel systems 85/86, and
c     the perkin elmer (interdata) 7/32.
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) /  z7fffffff /
c     data imach(10) /         16 /
c     data imach(11) /          6 /
c     data imach(12) /        -64 /
c     data imach(13) /         63 /
c     data imach(14) /         14 /
c     data imach(15) /        -64 /
c     data imach(16) /         63 /
c
c     machine constants for the ibm pc
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          0 /
c     data imach( 4) /          0 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1023 /
c
c     machine constants for the ibm rs 6000
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          0 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the intel i860
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the pdp-10 (ka processor)
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          5 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / "377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         54 /
c     data imach(15) /       -101 /
c     data imach(16) /        127 /
c
c     machine constants for the pdp-10 (ki processor)
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          5 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / "377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         62 /
c     data imach(15) /       -128 /
c     data imach(16) /        127 /
c
c     machine constants for pdp-11 fortran supporting
c     32-bit integer arithmetic.
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for pdp-11 fortran supporting
c     16-bit integer arithmetic.
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the silicon graphics
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the sun
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the sun
c     using the -r8 compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         53 /
c     data imach(12) /      -1021 /
c     data imach(13) /       1024 /
c     data imach(14) /        113 /
c     data imach(15) /     -16381 /
c     data imach(16) /      16384 /
c
c     machine constants for the univac 1100 series ftn compiler
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          1 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         60 /
c     data imach(15) /      -1024 /
c     data imach(16) /       1023 /
c
c     machine constants for the z80 microprocessor
c
c     data imach( 1) /          1 /
c     data imach( 2) /          1 /
c     data imach( 3) /          0 /
c     data imach( 4) /          1 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c***first executable statement  i1mach
      if (i .lt. 1  .or.  i .gt. 16) go to 10
c
      i1mach = imach(i)
      return
c
   10 continue
      write (unit = output, fmt = 9000)
 9000 format ('1error    1 in i1mach - i out of bounds')
c
c     call fdump
c
      stop
      end

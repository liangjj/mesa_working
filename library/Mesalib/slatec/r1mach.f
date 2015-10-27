*deck r1mach
      real function r1mach (i)
c***begin prologue  r1mach
c***purpose  return floating point machine dependent constants.
c***library   slatec
c***category  r1
c***type      single precision (r1mach-s, d1mach-d)
c***keywords  machine constants
c***author  fox, p. a., (bell labs)
c           hall, a. d., (bell labs)
c           schryer, n. l., (bell labs)
c***description
c
c   r1mach can be used to obtain machine-dependent parameters for the
c   local machine environment.  it is a function subprogram with one
c   (input) argument, and can be referenced as follows:
c
c        a = r1mach(i)
c
c   where i=1,...,5.  the (output) value of a above is determined by
c   the (input) value of i.  the results for various values of i are
c   discussed below.
c
c   r1mach(1) = b**(emin-1), the smallest positive magnitude.
c   r1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
c   r1mach(3) = b**(-t), the smallest relative spacing.
c   r1mach(4) = b**(1-t), the largest relative spacing.
c   r1mach(5) = log10(b)
c
c   assume single precision numbers are represented in the t-digit,
c   base-b form
c
c              sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c   where 0 .le. x(i) .lt. b for i=1,...,t, 0 .lt. x(1), and
c   emin .le. e .le. emax.
c
c   the values of b, t, emin and emax are provided in i1mach as
c   follows:
c   i1mach(10) = b, the base.
c   i1mach(11) = t, the number of base-b digits.
c   i1mach(12) = emin, the smallest exponent e.
c   i1mach(13) = emax, the largest exponent e.
c
c   to alter this function for a particular environment, the desired
c   set of data statements should be activated by removing the c from
c   column 1.  also, the values of r1mach(1) - r1mach(4) should be
c   checked for consistency with the local operating system.
c
c***references  p. a. fox, a. d. hall and n. l. schryer, framework for
c                 a portable library, acm transactions on mathematical
c                 software 4, 2 (june 1978), pp. 177-188.
c***routines called  xermsg
c***revision history  (yymmdd)
c   790101  date written
c   890213  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900618  added dec risc constants.  (wrb)
c   900723  added ibm rs 6000 constants.  (wrb)
c   910710  added hp 730 constants.  (smr)
c   911114  added convex ieee constants.  (wrb)
c   920121  added sun -r8 compiler option constants.  (wrb)
c   920229  added touchstone delta i860 constants.  (wrb)
c   920501  reformatted the references section.  (wrb)
c   920625  added convex -p8 and -pd8 compiler option constants.
c           (bks, wrb)
c   930201  added dec alpha and sgi constants.  (rwc and wrb)
c***end prologue  r1mach
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
c
      real rmach(5)
      save rmach
c
      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))
c
c     machine constants for the amiga
c     absoft fortran compiler using the 68020/68881 compiler option
c
c     data small(1) / z'00800000' /
c     data large(1) / z'7f7fffff' /
c     data right(1) / z'33800000' /
c     data diver(1) / z'34000000' /
c     data log10(1) / z'3e9a209b' /
c
c     machine constants for the amiga
c     absoft fortran compiler using software floating point
c
c     data small(1) / z'00800000' /
c     data large(1) / z'7effffff' /
c     data right(1) / z'33800000' /
c     data diver(1) / z'34000000' /
c     data log10(1) / z'3e9a209b' /
c
c     machine constants for the apollo
c
c     data small(1) / 16#00800000 /
c     data large(1) / 16#7fffffff /
c     data right(1) / 16#33800000 /
c     data diver(1) / 16#34000000 /
c     data log10(1) / 16#3e9a209b /
c
c     machine constants for the burroughs 1700 system
c
c     data rmach(1) / z400800000 /
c     data rmach(2) / z5ffffffff /
c     data rmach(3) / z4e9800000 /
c     data rmach(4) / z4ea800000 /
c     data rmach(5) / z500e730e8 /
c
c     machine constants for the burroughs 5700/6700/7700 systems
c
c     data rmach(1) / o1771000000000000 /
c     data rmach(2) / o0777777777777777 /
c     data rmach(3) / o1311000000000000 /
c     data rmach(4) / o1301000000000000 /
c     data rmach(5) / o1157163034761675 /
c
c     machine constants for the cdc 170/180 series using nos/ve
c
c     data rmach(1) / z"3001800000000000" /
c     data rmach(2) / z"4ffefffffffffffe" /
c     data rmach(3) / z"3fd2800000000000" /
c     data rmach(4) / z"3fd3800000000000" /
c     data rmach(5) / z"3fff9a209a84fbcf" /
c
c     machine constants for the cdc 6000/7000 series
c
c     data rmach(1) / 00564000000000000000b /
c     data rmach(2) / 37767777777777777776b /
c     data rmach(3) / 16414000000000000000b /
c     data rmach(4) / 16424000000000000000b /
c     data rmach(5) / 17164642023241175720b /
c
c     machine constants for the celerity c1260
c
c     data small(1) / z'00800000' /
c     data large(1) / z'7f7fffff' /
c     data right(1) / z'33800000' /
c     data diver(1) / z'34000000' /
c     data log10(1) / z'3e9a209b' /
c
c     machine constants for the convex
c     using the -fn compiler option
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7fffffff' /
c     data rmach(3) / z'34800000' /
c     data rmach(4) / z'35000000' /
c     data rmach(5) / z'3f9a209b' /
c
c     machine constants for the convex
c     using the -fi compiler option
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7f7fffff' /
c     data rmach(3) / z'33800000' /
c     data rmach(4) / z'34000000' /
c     data rmach(5) / z'3e9a209b' /
c
c     machine constants for the convex
c     using the -p8 or -pd8 compiler option
c
c     data rmach(1) / z'0010000000000000' /
c     data rmach(2) / z'7fffffffffffffff' /
c     data rmach(3) / z'3cc0000000000000' /
c     data rmach(4) / z'3cd0000000000000' /
c     data rmach(5) / z'3ff34413509f79ff' /
c
c     machine constants for the cray
c
c     data rmach(1) / 200034000000000000000b /
c     data rmach(2) / 577767777777777777776b /
c     data rmach(3) / 377224000000000000000b /
c     data rmach(4) / 377234000000000000000b /
c     data rmach(5) / 377774642023241175720b /
c
c     machine constants for the data general eclipse s/200
c     note - it may be appropriate to include the following card -
c     static rmach(5)
c
c     data small /    20k,       0 /
c     data large / 77777k, 177777k /
c     data right / 35420k,       0 /
c     data diver / 36020k,       0 /
c     data log10 / 40423k,  42023k /
c
c     machine constants for the dec alpha
c     using g_float
c
c     data rmach(1) / '00000080'x /
c     data rmach(2) / 'ffff7fff'x /
c     data rmach(3) / '00003480'x /
c     data rmach(4) / '00003500'x /
c     data rmach(5) / '209b3f9a'x /
c
c     machine constants for the dec alpha
c     using ieee_float
c
c     data rmach(1) / '00800000'x /
c     data rmach(2) / '7f7fffff'x /
c     data rmach(3) / '33800000'x /
c     data rmach(4) / '34000000'x /
c     data rmach(5) / '3e9a209b'x /
c
c     machine constants for the dec risc
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7f7fffff' /
c     data rmach(3) / z'33800000' /
c     data rmach(4) / z'34000000' /
c     data rmach(5) / z'3e9a209b' /
c
c     machine constants for the dec vax
c     (expressed in integer and hexadecimal)
c     the hex format below may not be suitable for unix systems
c     the integer format should be ok for unix systems
c
c     data small(1) /       128 /
c     data large(1) /    -32769 /
c     data right(1) /     13440 /
c     data diver(1) /     13568 /
c     data log10(1) / 547045274 /
c
c     data small(1) / z00000080 /
c     data large(1) / zffff7fff /
c     data right(1) / z00003480 /
c     data diver(1) / z00003500 /
c     data log10(1) / z209b3f9a /
c
c     machine constants for the elxsi 6400
c     (assuming real*4 is the default real)
c
c     data small(1) / '00800000'x /
c     data large(1) / '7f7fffff'x /
c     data right(1) / '33800000'x /
c     data diver(1) / '34000000'x /
c     data log10(1) / '3e9a209b'x /
c
c     machine constants for the harris 220
c
c     data small(1), small(2) / '20000000, '00000201 /
c     data large(1), large(2) / '37777777, '00000177 /
c     data right(1), right(2) / '20000000, '00000352 /
c     data diver(1), diver(2) / '20000000, '00000353 /
c     data log10(1), log10(2) / '23210115, '00000377 /
c
c     machine constants for the honeywell 600/6000 series
c
c     data rmach(1) / o402400000000 /
c     data rmach(2) / o376777777777 /
c     data rmach(3) / o714400000000 /
c     data rmach(4) / o716400000000 /
c     data rmach(5) / o776464202324 /
c
c     machine constants for the hp 730
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7f7fffff' /
c     data rmach(3) / z'33800000' /
c     data rmach(4) / z'34000000' /
c     data rmach(5) / z'3e9a209b' /
c
c     machine constants for the hp 2100
c     3 word double precision with ftn4
c
c     data small(1), small(2) / 40000b,       1 /
c     data large(1), large(2) / 77777b, 177776b /
c     data right(1), right(2) / 40000b,    325b /
c     data diver(1), diver(2) / 40000b,    327b /
c     data log10(1), log10(2) / 46420b,  46777b /
c
c     machine constants for the hp 2100
c     4 word double precision with ftn4
c
c     data small(1), small(2) / 40000b,       1 /
c     data large(1), large(2) / 77777b, 177776b /
c     data right(1), right(2) / 40000b,    325b /
c     data diver(1), diver(2) / 40000b,    327b /
c     data log10(1), log10(2) / 46420b,  46777b /
c
c     machine constants for the hp 9000
c
c     data small(1) / 00004000000b /
c     data large(1) / 17677777777b /
c     data right(1) / 06340000000b /
c     data diver(1) / 06400000000b /
c     data log10(1) / 07646420233b /
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9, the sel systems 85/86  and
c     the perkin elmer (interdata) 7/32.
c
c     data rmach(1) / z00100000 /
c     data rmach(2) / z7fffffff /
c     data rmach(3) / z3b100000 /
c     data rmach(4) / z3c100000 /
c     data rmach(5) / z41134413 /
c
c     machine constants for the ibm pc
c
c     data small(1) / 1.18e-38      /
c     data large(1) / 3.40e+38      /
c     data right(1) / 0.595e-07     /
c     data diver(1) / 1.19e-07      /
c     data log10(1) / 0.30102999566 /
c
c     machine constants for the ibm rs 6000
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7f7fffff' /
c     data rmach(3) / z'33800000' /
c     data rmach(4) / z'34000000' /
c     data rmach(5) / z'3e9a209b' /
c
c     machine constants for the intel i860
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7f7fffff' /
c     data rmach(3) / z'33800000' /
c     data rmach(4) / z'34000000' /
c     data rmach(5) / z'3e9a209b' /
c
c     machine constants for the pdp-10 (ka or ki processor)
c
c     data rmach(1) / "000400000000 /
c     data rmach(2) / "377777777777 /
c     data rmach(3) / "146400000000 /
c     data rmach(4) / "147400000000 /
c     data rmach(5) / "177464202324 /
c
c     machine constants for pdp-11 fortran supporting
c     32-bit integers (expressed in integer and octal).
c
c     data small(1) /    8388608 /
c     data large(1) / 2147483647 /
c     data right(1) /  880803840 /
c     data diver(1) /  889192448 /
c     data log10(1) / 1067065499 /
c
c     data rmach(1) / o00040000000 /
c     data rmach(2) / o17777777777 /
c     data rmach(3) / o06440000000 /
c     data rmach(4) / o06500000000 /
c     data rmach(5) / o07746420233 /
c
c     machine constants for pdp-11 fortran supporting
c     16-bit integers  (expressed in integer and octal).
c
c     data small(1), small(2) /   128,     0 /
c     data large(1), large(2) / 32767,    -1 /
c     data right(1), right(2) / 13440,     0 /
c     data diver(1), diver(2) / 13568,     0 /
c     data log10(1), log10(2) / 16282,  8347 /
c
c     data small(1), small(2) / o000200, o000000 /
c     data large(1), large(2) / o077777, o177777 /
c     data right(1), right(2) / o032200, o000000 /
c     data diver(1), diver(2) / o032400, o000000 /
c     data log10(1), log10(2) / o037632, o020233 /
c
c     machine constants for the silicon graphics
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7f7fffff' /
c     data rmach(3) / z'33800000' /
c     data rmach(4) / z'34000000' /
c     data rmach(5) / z'3e9a209b' /
c
c     machine constants for the sun
c
c     data rmach(1) / z'00800000' /
c     data rmach(2) / z'7f7fffff' /
c     data rmach(3) / z'33800000' /
c     data rmach(4) / z'34000000' /
c     data rmach(5) / z'3e9a209b' /
c
c     machine constants for the sun
c     using the -r8 compiler option
c
c     data rmach(1) / z'0010000000000000' /
c     data rmach(2) / z'7fefffffffffffff' /
c     data rmach(3) / z'3ca0000000000000' /
c     data rmach(4) / z'3cb0000000000000' /
c     data rmach(5) / z'3fd34413509f79ff' /
c
c     machine constants for the univac 1100 series
c
c     data rmach(1) / o000400000000 /
c     data rmach(2) / o377777777777 /
c     data rmach(3) / o146400000000 /
c     data rmach(4) / o147400000000 /
c     data rmach(5) / o177464202324 /
c
c     machine constants for the z80 microprocessor
c
c     data small(1), small(2) /     0,    256/
c     data large(1), large(2) /    -1,   -129/
c     data right(1), right(2) /     0,  26880/
c     data diver(1), diver(2) /     0,  27136/
c     data log10(1), log10(2) /  8347,  32538/
c
c***first executable statement  r1mach
      if (i .lt. 1 .or. i .gt. 5) call xermsg ('slatec', 'r1mach',
     +   'i out of bounds', 1, 2)
c
      r1mach = rmach(i)
      return
c
      end

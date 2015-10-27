*deck fundoc
      subroutine fundoc
c***begin prologue  fundoc
c***purpose  documentation for fnlib, a collection of routines for
c            evaluating elementary and special functions.
c***library   slatec
c***category  c, z
c***type      all (fundoc-a)
c***keywords  documentation, elementary functions, special functions
c***author  kahaner, d. k., (nbs)
c***description
c
c the slatec library --  elementary and special functions
c
c this describes the elementary and special function routines available
c in the slatec library.  most of the these routines were written by
c wayne fullerton while at lanl.  some were written by don amos of snla.
c there are approximately 63 single precision, 63 double precision and
c 25 complex user callable elementary and special function routines.
c
c the table below gives a breakdown of routines according to their
c function.  unless otherwise indicated all routines are function
c subprograms.
c                                             sngl.      dble.
c description              notation           prec.      prec.   complex
c
c         ***intrinsic functions and fundamental functions***
c unpack floating point              call r9upak(x,y,n)  d9upak    --
c  number
c pack floating point                        r9pak(y,n)  d9pak     --
c  number
c initialize orthogonal               inits(os,nos,eta)  initds    --
c  polynomial series
c evaluate chebyshev       summation for  csevl(x,cs,n)  dcsevl    --
c series                  i = 1 to n of
c                          cs(i)*(2*x)**(i-1)
c
c                  ***elementary functions***
c argument = theta in      z = \ z \ *          --         --    carg(z)
c  radians                 e**(i * theta)
c cube root                                   cbrt(x)    dcbrt   ccbrt
c relative error exponen-  ((e**x) -1) / x    exprel(x)  dexprl  cexprl
c  tial from first order
c common logarithm         log to the base 10   --         --  clog10(z)
c                          of z
c relative error logarithm ln(1 + x)          alnrel(x)  dlnrel  clnrel
c relative error logarithm (ln(1 + x) - x     r9ln2r(x)  d9ln2r  c9ln2r
c from second order        + x**2/2) / x**3
c               ***trigonometric and hyperbolic functions***
c tangent                  tan z                --         --    ctan(z)
c cotangent                cot x              cot(x)     dcot    ccot
c sine x in degrees        sin((2*pi*x)/360)  sindg(x)   dsindg    --
c cosine x in degrees      cos((2*pi*x)/360)  cosdg(x)   dcosdg    --
c arc sine                 arcsin (z)           --         --   casin(z)
c arc cosine               arccos (z)           --         --   cacos(z)
c arc tangent              arctan (z)           --         --   catan(z)
c quadrant correct         arctan (z1/z2)       --         -- catan2(z1,
c  arc tangent                                                       z2)
c hyperbolic sine          sinh z               --         --   csinh(z)
c hyperbolic cosine        cosh z               --         --   ccosh(z)
c hyperbolic tangent       tanh z               --         --   ctanh(z)
c arc hyperbolic sine      arcsinh (x)        asinh(x)   dasinh  casinh
c arc hyperbolic cosine    arccosh (x)        acosh(x)   dacosh  cacosh
c arc hyperbolic tangent   arctanh (x)        atanh(x)   datanh  catanh
c relative error arc       (arctan (x) - x)   r9atn1(x)  d9atn1    --
c  tangent from first order   / x**3
c              ***exponential integrals and related functions***
c exponential integral     ei(x) = (minus)    ei(x)      dei       --
c                          the integral from
c                          -x to infinity of
c                            (e**-t / t)dt
c exponential integral     e sub 1 (x) =      e1(x)      de1       --
c                          the integral from x
c                            to infinity of
c                          (e**-t / t) dt
c logarithmic integral     li(x) = the        ali(x)     dli       --
c                          integral from 0 to
c                          x of (1 / ln t) dt
c   sequences of exponential integrals.
c   m values are computed where
c   k=0,1,...m-1 and n>=1
c exponential integral     e sub n+k (x) call exint(x,   dexint    --
c                        =the integral from   n,kode,m,tol,
c                         1 to infinity of    en,ierr)
c                       (e**(-x*t)/t**(n+k))dt
c                 ***gamma functions and related functions***
c factorial                n!                 fac(n)     dfac      --
c binomial                 n!/(m!*(n-m)!)     binom(n,m) dbinom    --
c gamma                    gamma(x)           gamma(x)   dgamma  cgamma
c gamma(x) under and                     call gamlim(    dgamlm    --
c  overflow limits                           xmin,xmax)
c reciprocal gamma         1 / gamma(x)       gamr(x)    dgamr   cgamr
c log abs gamma            ln \gamma(x)\      alngam(x)  dlngam    --
c log gamma                ln gamma(z)          --         --    clngam
c log abs gamma       g = ln \gamma(x)\  call algams(x,  dlgams    --
c with sign           s = sign gamma(x)      g,s)
c incomplete gamma         gamma(a,x) =       gami(a,x)  dgami     --
c                          the integral from
c                          0 to x of
c                         (t**(a-1) * e**-t)dt
c complementary            gamma(a,x) =       gamic(a,x) dgamic    --
c  incomplete gamma        the integral from
c                          x to infinity of
c                         (t**(a-1) * e**-t)dt
c tricomi's             gamma super star(a,x) gamit(a,x) dgamit    --
c  incomplete gamma        = x**-a *
c                         incomplete gamma(a,x)
c                          / gamma(a)
c psi (digamma)            psi(x) = gamma'(x) psi(x)     dpsi    cpsi
c                          / gamma(x)
c pochhammer's         (a) sub x = gamma(a+x) poch(a,x)  dpoch     --
c  generalized symbol      / gamma(a)
c pochhammer's symbol    ((a) sub x -1) / x   poch1(a,x) dpoch1    --
c  from first order
c beta                     b(a,b) = (gamma(a) beta(a,b)  dbeta   cbeta
c                          * gamma(b))
c                          / gamma(a+b)
c                           = the integral
c                           from 0 to 1 of
c                           (t**(a-1) *
c                           (1-t)**(b-1))dt
c log beta                 ln b(a,b)         albeta(a,b) dlbeta  clbeta
c incomplete beta          i sub x (a,b) =  betai(x,a,b) dbetai    __
c                          b sub x (a,b) / b(a,b)
c                           = 1 / b(a,b) *
c                          the integral
c                          from 0 to x of
c                          (t**(a-1) *
c                          (1-t)**(b-1))dt
c log gamma correction     ln gamma(x) -      r9lgmc(x)  d9lgmc  c9lgmc
c  term when stirling's    (ln(2 * pi))/2 -
c  approximation is valid  (x - 1/2) * ln(x) + x
c                ***error functions and fresnel integrals***
c error function           erf x = (2 /       erf(x)     derf      --
c                          square root of pi) *
c                          the integral from
c                          0 to x of
c                          e**(-t**2)dt
c complementary            erfc x = (2 /      erfc(x)    derfc     --
c  error function          square root of pi) *
c                          the integral from
c                          x to infinity of
c                          e**(-t**2)dt
c dawson's function        f(x) = e**(-x**2)  daws(x)    ddaws     --
c                          * the integral from
c                          from 0 to x of
c                          e**(t**2)dt
c                         ***bessel functions***
c   bessel functions of special integer order
c first kind, order zero   j sub 0 (x)        besj0(x)   dbesj0    --
c first kind, order one    j sub 1 (x)        besj1(x)   dbesj1    --
c second kind, order zero  y sub 0 (x)        besy0(x)   dbesy0    --
c second kind, order one   y sub 1 (x)        besy1(x)   dbesy1    --
c   modified (hyperbolic) bessel functions of special integer order
c first kind, order zero   i sub 0 (x)        besi0(x)   dbesi0    --
c first kind, order one    i sub 1 (x)        besi1(x)   dbesi1    --
c third kind, order zero   k sub 0 (x)        besk0(x)   dbesk0    --
c third kind, order one    k sub 1 (x)        besk1(x)   dbesk1    --
c   modified (hyperbolic) bessel functions of special integer order
c   scaled by an exponential
c first kind, order zero   e**-\x\ * i sub 0(x) besi0e(x) dbsi0e   --
c first kind, order one    e**-\x\ * i sub 1(x) besi1e(x) dbsi1e   --
c third kind, order zero   e**x * k sub 0 (x)   besk0e(x) dbsk0e   --
c third kind, order one    e**x * k sub 1 (x)   besk1e(x) dbsk1e   --
c   sequences of bessel functions of general order.
c   n values are computed where  k = 1,2,...n and v .ge. 0.
c modified first kind      i sub v+k-1 (x) call besi(x,   dbesi    --
c                          optional scaling  alpha,kode,n,
c                          by e**(-x)        y,nz)
c first kind               j sub v+k-1 (x) call besj(x,   dbesj    --
c                                            alpha,n,y,nz)
c second kind              y sub v+k-1 (x) call besy(x,   dbesy    --
c                                            fnu,n,y)
c modified third kind      k sub v+k-1 (x) call besk(x,   dbesk    --
c                          optional scaling  fnu,kode,n,y,
c                          by e**(x)         nz)
c   sequences of bessel functions.  \n\ values are computed where
c   i = 0, 1, 2, ..., n-1  for n > 0  or i = 0, -1, -2, ..., n+1
c   for n < 0.
c modified third kind      k sub v+i (x)   call besks(    dbesks   --
c                                           xnu,x,n,bk)
c   sequences of bessel functions scaled by an exponential.
c   \n\ values are computed where  i = 0, 1, 2, ..., n-1
c   for n > 0  or  i = 0, -1, -2, ..., n+1  for n < 0.
c modified third kind      e**x *         call beskes(    dbskes   --
c                          k sub v+i (x)     xnu,x,n,bk)
c                ***bessel functions of fractional order***
c   airy functions
c airy                     ai(x)              ai(x)      dai       --
c bairy                    bi(x)              bi(x)      dbi       --
c   exponentially scaled airy functions
c airy                     ai(x), x <= 0      aie(x)     daie      --
c                          exp(2/3 * x**(3/2))
c                          * ai(x), x >= 0
c bairy                    bi(x), x <= 0      bie(x)     dbie      --
c                          exp(-2/3 * x**(3/2))
c                          * bi(x), x >= 0
c                 ***confluent hypergeometric functions***
c confluent                u(a,b,x)           chu(a,b,x) dchu      --
c  hypergeometric
c                     ***miscellaneous functions***
c spence                   s(x) = - the       spenc(x)   dspenc    --
c  dilogarithm             integral from
c                          0 to x of
c                          ((ln \1-y\) / y)dy
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   801015  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  routine name changed from fnlibd to fundoc.  (wrb)
c   900723  purpose section revised.  (wrb)
c***end prologue  fundoc
c***first executable statement  fundoc
      return
      end

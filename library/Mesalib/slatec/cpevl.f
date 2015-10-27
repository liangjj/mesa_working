*deck cpevl
      subroutine cpevl (n, m, a, z, c, b, kbd)
c***begin prologue  cpevl
c***subsidiary
c***purpose  subsidiary to cpzero
c***library   slatec
c***type      single precision (cpevl-s)
c***author  (unknown)
c***description
c
c        evaluate a complex polynomial and its derivatives.
c        optionally compute error bounds for these values.
c
c   input...
c        n = degree of the polynomial
c        m = number of derivatives to be calculated,
c            m=0 evaluates only the function
c            m=1 evaluates the function and first derivative, etc.
c             if m .gt. n+1 function and all n derivatives will be
c                calculated.
c       a = complex vector containing the n+1 coefficients of polynomial
c               a(i)= coefficient of z**(n+1-i)
c        z = complex point at which the evaluation is to take place.
c        c = array of 2(m+1) words into which values are placed.
c        b = array of 2(m+1) words only needed if bounds are to be
c              calculated.  it is not used otherwise.
c        kbd = a logical variable, e.g. .true. or .false. which is
c              to be set .true. if bounds are to be computed.
c
c  output...
c        c =  c(i+1) contains the complex value of the i-th
c              derivative at z, i=0,...,m
c        b =  b(i) contains the bounds on the real and imaginary parts
c              of c(i) if they were requested.
c
c***see also  cpzero
c***routines called  i1mach
c***revision history  (yymmdd)
c   810223  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cpevl
c
      complex a(*),c(*),z,ci,cim1,b(*),bi,bim1,t,za,q
      logical kbd
      save d1
      data d1 /0.0/
      za(q)=cmplx(abs(real(q)),abs(aimag(q)))
c***first executable statement  cpevl
      if (d1 .eq. 0.0) d1 = real(i1mach(10))**(1-i1mach(11))
      np1=n+1
      do 1 j=1,np1
         ci=0.0
         cim1=a(j)
         bi=0.0
         bim1=0.0
         mini=min(m+1,n+2-j)
            do 1 i=1,mini
               if(j .ne. 1) ci=c(i)
               if(i .ne. 1) cim1=c(i-1)
               c(i)=cim1+z*ci
               if(.not. kbd) go to 1
               if(j .ne. 1) bi=b(i)
               if(i .ne. 1) bim1=b(i-1)
               t=bi+(3.*d1+4.*d1*d1)*za(ci)
               r=real(za(z)*cmplx(real(t),-aimag(t)))
               s=aimag(za(z)*t)
               b(i)=(1.+8.*d1)*(bim1+d1*za(cim1)+cmplx(r,s))
               if(j .eq. 1) b(i)=0.0
    1 continue
      return
      end

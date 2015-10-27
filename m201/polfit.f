*deck @(#)polfit.f	5.1  11/6/94
      subroutine polfit(n,ndiff,x,values,xmat,iscr,scr,a,ok)
c***begin prologue     polfit.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)polfit.f	5.1   11/6/94
c***purpose            
c***description
c     this subroutine fits a polynomial to n values of a function
c     and its first (ndiff-1, ndiff.ge.1) derivatives.
c     the polynomial is of degree (n*ndiff=nprod) -1.  x(n) has the
c     ordinates and values(n,ndiff) has the values of the function and
c     its derivatives.  xmat(nprod,nprod) is a scratch matrix, as are
c     iscr(nprod) and scr(nprod).  a is returned with the coefficients
c     of the polynomial
c        f=a(1)*x**(nprod-1) +... + a(nprod)
c     ok is returned .true. if the fit was successful.
c     
c***references
c
c***routines called
c
c***end prologue       polfit.f
      implicit none
c     --- input variables -----
      integer n,ndiff
c     --- input arrays (unmodified) ---
      real*8 x(n),values(n,ndiff)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 a(n,ndiff)
c     --- output variables ---
      logical ok
c     --- scratch arrays ---
      integer iscr(n*ndiff)
      real*8 xmat(n,ndiff,n*ndiff),scr(n*ndiff)
c     --- local variables ---
      integer nprod,i,j,k,exp
      real*8 prod,zero,one,det
      parameter (zero=0.0d+00,one=1.0d+00)
c
c
c     --- load xmat with the powers of x at the specified points.
c         the first n rows reflect the function values, and the next
c         n the first derivatives, etc.
      nprod=n*ndiff
      do 30 i=1,n
         do 20 j=1,nprod
            prod=one
            exp=nprod-j
            do 10 k=1,ndiff
               if(exp.lt.0) then
                  xmat(i,k,j)=zero
               else if(exp.eq.0) then
                  xmat(i,k,j)=prod
               else
                  xmat(i,k,j)=prod*x(i)**exp
                  prod=prod*float(exp)
               endif
               exp=exp-1
   10       continue
   20    continue
   30 continue
c
c     --- invert and multiply to form the coefficients.
      call minvrt(xmat,nprod,nprod,det,iscr,scr)
      call ebc(a,xmat,values,nprod,nprod,1)
      ok=.true.
c
c
      return
      end

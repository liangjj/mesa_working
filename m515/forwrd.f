*deck @(#)forwrd.f	5.1  11/28/95
deck @(#)forwrd.f	5.1   11/28/95
      subroutine forwrd(flm,j,y,wt,int0,psilm,n)
c***begin prologue     forwrd.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             schneider, barry(nsf)
c***source             @(#)forwrd.f	5.1   11/28/95
c***purpose            forward integration of indefinite integral
c***description
c   performs the set of indefinite integrals which result
c   from integrating a function on [a,b] as [a,r(1)], [r(1),r(2)]
c   [r(2),r(3)]....[a,b] where r(i) are the quadrature points in
c   the interval. the integral is initialized as int0
c   which is either its last value or zero depending on the
c   interval. 
c
c   on exit, 
c      psilm(1)=[a,r1] and psilm(n-1)=[a,b]   
c
c***references
c
c***routines called
c
c***end prologue       forwrd.f
      implicit none
c     --- input variables -----
      integer n
      real*8 int0
c     --- input arrays (unmodified) ---
      real*8 flm(n),j(n),y(n),wt(n,n-1)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 psilm(n-1)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer pti,ptj
      real*8 sumf
c
      common /io/ inp, iout
c
      sumf=int0
      do 20 pti=1,n-1
         do 10 ptj=1,n
             sumf=sumf+wt(ptj,pti)*flm(ptj)*j(ptj)
   10    continue
         psilm(pti)=sumf
   20 continue
c
      int0=sumf
      do 30 pti=1,n-1
         psilm(pti)=psilm(pti)*y(pti+1)
   30 continue
c
c
      return
      end

*deck @(#)vmkghi.f	5.1 11/6/94
      subroutine vmkghi(g,h,i2,f00,b00,b10,bp01,c00,cp00,nmax,mmax,
     $                  imax,jmax,kmax,lmax,c,nat,symcen,lenv)
c***begin prologue     vmkghi.f
c***date written       851017  (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul
c***source             @(#)vmkghi.f	5.1 11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       vmkghi.f
      implicit none
c     --- input variables -----
      integer nmax,mmax,imax,jmax,kmax,lmax,nat,lenv
c     --- input arrays (unmodified) ---
      integer symcen(4)
      real*8 f00(lenv),b00(lenv),b10(lenv),bp01(lenv),c00(lenv,3)
      real*8 cp00(lenv,3)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
      real*8 g(lenv,3,0:nmax,0:mmax),h(lenv,3,0:nmax,0:jmax,0:mmax)
c     --- output arrays ---
      real*8 i2(lenv,3,0:imax,0:jmax,0:mmax,0:lmax)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
c
      call vmakg(g,f00,b00,b10,bp01,c00,cp00,nmax,mmax,lenv,lenv)
      call vmakh(h,g,c,nat,symcen,nmax,mmax,imax,jmax,lenv)
      call vmaki2(i2,h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,lenv)
c
c
      return
      end

*deck @(#)vmkghi.f	5.1  11/28/95
      subroutine vmkghi(g,h,i2,f00,b00,b10,bp01,c00,cp00,nmax,mmax,
     $                  imax,jmax,kmax,lmax,c,nat,symcen,lenv)
c***begin prologue     vmkghi.f
c***date written       851017  (yymmdd)  
c***revision date      4/18/95      
c
c***keywords           
c***author             saxe, paul
c***source             @(#)vmkghi.f	5.1   11/28/95
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
      integer inp,iout
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timg,timh,timi2
      logical timeit
      common/io/inp,iout
      parameter (timeit=.false.)
      data timg,timh,timi2/3*0.0d0/

      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
      call vmakg(g,f00,b00,b10,bp01,c00,cp00,nmax,mmax,lenv,lenv)
      if(timeit) then
         call timing(dum4,dum5,dum6)
         timg=timg+dum4-dum1
      endif
      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
      call vmakh(h,g,c,nat,symcen,nmax,mmax,imax,jmax,lenv)
      if(timeit) then
         call timing(dum4,dum5,dum6)
         timh=timh+dum4-dum1
      endif
      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
      call vmaki2(i2,h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,lenv)
      if(timeit) then
         call timing(dum4,dum5,dum6)
         timi2=timi2+dum4-dum1
      endif
c
c
      if(timeit) then
         write(iout,*) 'timg,timh,timi2',timg,timh,timi2
      endif
c
c
      return
      end

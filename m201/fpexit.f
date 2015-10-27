*deck @(#)fpexit.f	5.1  11/6/94
      subroutine fpexit(nvar,nz,toang,lbl,pool0,delvar,delloc,varloc,
     $                  vname,fzero)
c***begin prologue     fpexit.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkely, et al. (g82)
c***source             @(#)fpexit.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fpexit.f
      implicit none
c     --- input variables -----
      integer nvar,nz
      real*8 toang,fzero
c     --- input arrays (unmodified) ---
      integer lbl(nz)
      character*(*) vname(nvar)
      real*8 pool0(nvar),delvar(nvar),delloc(nvar),varloc(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j
      logical bond,angle
      real*8 todeg,one,f45
c
      parameter (f45=4.5d+01,one=1.0d+00)
c
      common/io/inp,iout
c
 1000 format(5x,'final pool of variables:'
     $      /10x,'name',13x,'value',8x,'increment',
     $      /10x,'----',13x,'-----',8x,'---------')
 1010 format(5x,2x,a16,f10.5,5x,f10.5)
 1020 format(5x,'final optimized value = ',e20.13)
c
c     --- convert pool0 and delvar from bohr/radian units to
c         angstrom/degree units for printing.
      call iosys('read integer zlbl from rwf',nz,lbl,0,' ')
      todeg=f45/atan(one)
      do 100 i=1,nvar
         bond=.false.
         angle=.false.
         do 80 j=2,nz
            if(abs(lbl(j)).eq.i) bond=.true.
   80    continue
         if(bond) then
            delloc(i)=delvar(i)*toang
            varloc(i)=pool0(i)*toang
         else
            delloc(i)=delvar(i)*todeg
            varloc(i)=pool0(i)*todeg
         endif
  100 continue
c
c     --- print final pool of variables and energy.
      write(iout,1000)
      do 110 i=1,nvar
         write(iout,1010) vname(i),varloc(i),delloc(i)
  110 continue
      write(iout,1020) fzero
c
c
      return
      end

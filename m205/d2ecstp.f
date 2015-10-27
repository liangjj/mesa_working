*deck @(#)d2ecstp.f	5.1  11/6/94
      subroutine d2ecstp(maxpt,xc,fc,xxc,ffc,d2edone,xnew,
     $                   natoms)
c***begin prologue     d2ecstp.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)d2ecstp.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       d2ecstp.f
      implicit none
c     --- input variables -----
      integer maxpt,natoms
c     --- input arrays (unmodified) ---
      real*8 fc(natoms*3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 xc(natoms*3),xxc(natoms*3,maxpt)
      real*8 ffc(natoms*3,maxpt)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 xnew(natoms*3)
c     --- local variables ---
      integer inp,iout
      integer d2ecycl
      integer i,jj,nat3
      logical prnt,chkpt,singpt,cartfx
      logical debug,d2edone
      real*8 energy,rmax,rmin,rlim,stpsize
c
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common/io/inp,iout
c
c
 1020 format(5x,'cycle number ',i3,' out of a maximum of ',i3,
     $      /5x,'all quantities printed in hartrees-bohr-radians.')
 1030 format(5x,'step:',i3,'      gradient   step    coordinates' )
 1060 format(3x,a16,2x,6f10.5)
 1070 format(' coordinates and gradient')
 1080 format(2f20.10)
 1090 format(/'moving variable ',i2,' by ', f5.3)
c
c     --- increment the cycle counter
      nat3=natoms*3
      d2ecycl=d2ecycl+1
      if(debug) then
         write(iout,1020) d2ecycl,maxpt
      endif
c
c     --- put current coordinates and gradient into xxc and ffc
      do 50 i=1,nat3
         xxc(i,d2ecycl)=xc(i)
         ffc(i,d2ecycl)=fc(i)
   50 continue
      if(debug) then
         write(iout,1070)
         do 60 jj=1,nat3
            write(iout,1080) xxc(jj,d2ecycl),ffc(jj,d2ecycl)
   60    continue
      endif
      call rzero(xnew,nat3)
c
c     --- determine the step vector
      if(singpt) then
c        --- single point difference formula
c            increment a coordinate
         if(d2ecycl.le.nat3) then
            xnew(d2ecycl) = stpsize
            if(debug) then
              write(iout,1090) d2ecycl, stpsize
            endif
         else
            d2edone=.true.
         endif
      else
c        --- double point difference formula
         if(d2ecycl.le.(2*nat3)) then
            if(mod(d2ecycl,2).eq.1) then
               xnew((d2ecycl+1)/2)=stpsize
            else
               xnew(d2ecycl/2)=-stpsize
            endif
         else
            d2edone=.true.
         endif
      endif
c
      call vadd(xc,xxc(1,1),xnew,nat3)
      if(d2edone) then
c        --- put the unperturbed coordinates back in x
         do 70 i=1,nat3
            xc(i)=xxc(i,1)
   70    continue
      endif
c
c
      return
      end

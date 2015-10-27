*deck @(#)d2estep.f	5.1 11/6/94
      subroutine d2estep(nz,nvar,nvv,maxpt,toang,vname,x,f,xx,
     $                  ff,d2edone,xnew)
c***begin prologue     d2estep.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)d2estep.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       d2estep.f
      implicit none
c     --- input variables -----
      integer nz,nvar,nvv,maxpt
      real*8 toang
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
      real*8 x(nvar),f(nvar)
c     --- input arrays (scratch) ---
      real*8 xnew(nvar)
c     --- output arrays ---
      real*8 xx(nvar,maxpt),ff(nvar,maxpt)
c     --- output variables ---
      logical d2edone
c     --- scratch arrays ---
c     --- local variables ---
      integer d2ecycl
      integer inp,iout
      integer i
      character tname*16
      logical prnt,chkpt,singpt,cartfx
      real*8 energy,rmax,rmin,rlim,stpsize
      logical debug
c
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common/io/inp,iout
c
 1020 format(5x,'cycle number ',i3,' out of a maximum of ',i3,
     $      /5x,'all quantities printed in hartrees-bohr-radians.')
 1030 format(5x,'step:',i3,'      gradient   step    coordinates' )
 1060 format(3x,a16,2x,6f10.5)
c
c     --- increment the cycle counter
      d2ecycl=d2ecycl+1
      if(debug) then
         write(iout,1020) d2ecycl,maxpt
      endif
c     --- put current coordinates and gradient into xx and ff
      do 55 i=1,nvar
         xx(i,d2ecycl)=x(i)
         ff(i,d2ecycl)=f(i)
   55 continue
      call rzero(xnew,nvar)
c
c     --- determine the step vector
      if(singpt) then
c        ---single point difference formula
c        increment a coordinate
         if(d2ecycl.le.nvar) then
            xnew(d2ecycl) = stpsize
c
            if(debug) then
               write(iout,1300) d2ecycl, stpsize
 1300          format(/'moving variable ',i2,' by ', f5.3)
            endif
        else
           d2edone=.true.
        endif
      else
c        --- double point difference formula
         if(d2ecycl.le.(2*nvar)) then
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
      call vadd(x,xx(1,1),xnew,nvar)
      if(debug) then
         write(iout,1030)d2ecycl
         do 40 i=1,nvar
            call crjust(vname(i),tname)
            write(iout,1060) tname,f(i),xnew(i),x(i)
   40    continue
      endif
c
      if(d2edone) then
c        --- put the unperturbed coordinates back in x
         do 90 i=1,nvar
            x(i)=xx(i,1)
   90    continue
      endif
c
c
      return
      end

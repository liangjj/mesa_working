*deck @(#)d2eform.f	5.1  11/6/94
      subroutine d2eform(nvar,nvv,maxpt,frcnst,xx,ff,zsq,vname)
c***begin prologue     d2eform.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)d2eform.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       d2eform.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,maxpt
c     --- input arrays (unmodified) ---
      real*8 xx(nvar,maxpt),ff(nvar,maxpt)
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 frcnst(nvv),zsq(nvar,nvar)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer d2ecycl
      integer inp,iout
      integer i,j,k,ij
      logical prnt,chkpt,singpt,cartfx
      real*8 energy,rmax,rmin,rlim,stpsize
      logical debug
c
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common /io/ inp,iout
c
 1000    format(//' ready to make force constant matrix '/
     $            'coordinates and gradients'//)
 1005 format(5x,'finite difference second derivative calculation:')
 1010 format(5x,'single point difference formula: stepsize = ',f5.3)
 1015 format(5x,'double point difference formula: stepsize = ',f5.3)
 1020 format(8f8.4)
 1030 format(//)
 1040 format(3f10.5)
 1050 format(//'     symmetrized force constants ',/,
     $         '  i    j   ij    frcnst(ij) '/)
 1060 format(3i5,f12.6)
 1070 format(5x,'the internal coordinate force constant matrix')
c
c
      if(debug) then
         write(iout,1000)
         call matout(xx,nvar,maxpt,nvar,maxpt,iout)
         call matout(ff,nvar,maxpt,nvar,maxpt,iout)
         do 10 j=1,maxpt
            write(iout,1020)(xx(j,i),i=1,nvar+1),(ff(j,i),i=1,nvar+1)
   10    continue
         write(iout,1030)
      endif
c
c     --- branch on single vs. double point differencing.
      if(singpt) then
         do 30 j=1,nvar
            do 20 i=1,nvar
c              --- the minus sign in the difference formula is
c                  because these are forces, not gradients
               zsq(i,j)=-(ff(i,j+1)-ff(i,1))/stpsize
   20       continue
            if(debug) then
               write(iout,1040) (zsq(k,j),k=1,nvar)
            endif
   30    continue
         if(debug) then
            write(iout,1050)
         endif
c        --- symmetrize the force constant matrix
         ij=0
         do 40 j=1,nvar
            do 40 i=1,j
               ij=ij+1
               frcnst(ij)=(zsq(i,j)+zsq(j,i))/2.0
               if(debug) write(iout,1060) i,j,ij,frcnst(ij)
   40    continue
      else
c        --- double point
         do 60 j=1,nvar
            do 50 i=1,nvar
c              --- the minus sign in the difference formula is
c                  because these are forces, not gradients
               zsq(i,j)=-(ff(i,2*j)-ff(i,2*j+1))/(2*stpsize)
   50    continue
             if(debug) write(iout,1040) (zsq(k,j),k=1,nvar)
   60    continue
         if(debug) write(iout,1050)
c        --- symmetrize the force constant matrix
         ij=0
         do 70 j=1,nvar
            do 70 i=1,j
               ij=ij+1
               frcnst(ij)=(zsq(i,j)+zsq(j,i))/2.0
               if(debug)  write(iout,1060) i,j,ij,frcnst(ij)
   70    continue
      endif
c
c     --- reprot the results
      write(iout,1005)
      if(singpt) then
         write(iout,1010) stpsize
      else
           write(iout,1015) stpsize
      endif
      write(iout,1070)
      call frqprt(zsq,nvar,nvar,nvar,nvar,1,1,vname,vname,
     $            0,frcnst,.false.)
c
c
      return
      end

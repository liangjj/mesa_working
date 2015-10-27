*deck @(#)d2ecfrm.f	5.1  11/6/94
      subroutine d2ecfrm(nat3,nnsq,maxpt,xxc,ffc,nvv,frcnst,zsq)
c***begin prologue     d2ecfrm.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)d2ecfrm.f	5.1   11/6/94
c***purpose            forms cartesian force constants by finite difference.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       d2ecfrm.f
      implicit none
c     --- input variables -----
      integer nat3,nnsq,maxpt,nvv
c     --- input arrays (unmodified) ---
      real*8 ffc(nat3,maxpt),xxc(nat3,maxpt)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 frcnst(nvv),zsq(nat3,nat3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer d2ecycl
      integer inp,iout
      integer i,j,k,ij,nnvc
      logical prnt,chkpt,singpt,cartfx
      logical debug
      real*8 energy,rmax,rmin,rlim,stpsize
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common /io/ inp,iout
c
 1000 format(//' ready to make force constant matrix '/
     $            'coordinates and gradients'//)
 1005 format(5x,'finite difference second derivative calculation:')
 1010 format(5x,'single point difference formula: stepsize = ',f5.3)
 1015 format(5x,'double point difference formula: stepsize = ',f5.3)
 1018 format(8f8.4)
 1020 format(3f10.5)
 1030 format(//'     symmetrized force constants ',/,
     $         '  i    j   ij    frcnst(ij) '/)
 1040 format(3i5,f12.6)
 1050 format(/,5x,'the force constant matrix')
c
c
      if(debug) then
         write(iout,1000)
         call matout(xxc,nat3,maxpt,nat3,maxpt,iout)
         call matout(ffc,nat3,maxpt,nat3,maxpt,iout)
         do 10 j=1,maxpt
            write(iout,1018)(xxc(j,i),i=1,nat3+1),(ffc(j,i),i=1,nat3+1)
   10    continue
      endif
c
c     --- divide depending on single vs. double point differencing
      if(singpt) then
         do 30 j=1,nat3
            do 20 i=1,nat3
               zsq(i,j)=(ffc(i,j+1)-ffc(i,1))/stpsize
   20       continue
            if(debug) then
               write(iout,1020) (zsq(k,j),k=1,nat3)
            endif
   30    continue
         if(debug) then
            write(iout,1030)
         endif
c
c        --- symmetrize the force constant matrix
         ij=0
         do 40 j=1,nat3
            do 40 i=1,j
               ij=ij+1
               frcnst(ij)=(zsq(i,j)+zsq(j,i))/2.0
               zsq(i,j)=frcnst(ij)
               zsq(j,i)=frcnst(ij)
c              if(debug)  write(iout,1040) i,j,ij,frcnst(ij)
   40    continue
      else
c
c        --- double point
         do 60 j=1,nat3
            do 50 i=1,nat3
               zsq(i,j)=(ffc(i,2*j)-ffc(i,2*j+1))/(2*stpsize)
   50    continue
            if(debug) write(iout,1020) (zsq(k,j),k=1,nat3)
   60    continue
c
c        --- symmetrize the force constant matrix
         if(debug) write(iout,1030)
         ij=0
         do 70 j=1,nat3
            do 70 i=1,j
               ij=ij+1
               frcnst(ij)=(zsq(i,j)+zsq(j,i))/2.0
               zsq(i,j)=frcnst(ij)
               zsq(j,i)=frcnst(ij)
               if(debug) write(iout,1040) i,j,ij,frcnst(ij)
   70    continue
      endif
c
c     --- report findings.
      write(iout,1005)
      if(singpt) then
         write(iout,1010) stpsize
      else
         write(iout,1015) stpsize
      endif
      write(iout,1050)
      call frqprt(zsq,nat3,nat3,nat3,nat3,0,0,' ',' ',
     $            0,frcnst,.false.)
c
c     --- update the second derivatives on disk.
      nnvc=nat3*(nat3+1)/2
      call iosys('write real "cartesian second derivatives" to rwf',
     $            nnvc,frcnst,0,' ')
c
c
      return
      end

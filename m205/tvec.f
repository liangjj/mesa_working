*deck @(#)tvec.f	5.1 11/6/94
      subroutine tvec(t,c,ixyz,lp,lm,natoms3,natoms,t1,t2,nrt)
c***begin prologue     tvec.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)tvec.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       tvec.f
      implicit none
c     --- input variables -----
      integer ixyz,natoms,natoms3,nrt
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 t(natoms3,6),lp(natoms3,6),lm(natoms3,6),c(3,natoms)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 t1(natoms3,natoms3),t2(natoms3)
c     --- local variables ---
      integer inp,iout
      integer d2ecycl
      integer i,j,ii,jj,atom,katom,kxyz,ix,iy,iz
      logical prnt,chkpt,singpt,cartfx
      logical debug
      real*8 energy,rmax,rmin,rlim,stpsize
      real*8 zero,one,two
c
      parameter (debug=.false.)
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $              prnt,chkpt,singpt,stpsize,cartfx
      common /io/ inp,iout
c
 1000 format('  tvec: ixyz,kxyz,katom= ',3i5)
 1010 format(' in tvec bet rot and trans ')
 1020 format(' in tvec: calling schmidt nrt=',i5)
 1030 format(' in tvec: after schmidt ')
c
c
      ii=0
      do 40 atom=1,natoms
         do 35 i=1,3
            ii=ii+1
            if(ii.eq.ixyz) then
               kxyz=i
               katom=atom
            endif
   35    continue
   40 continue
      if(debug) then
         write(iout,1000) ixyz,kxyz,katom
      endif
      c(kxyz,katom)=c(kxyz,katom)+stpsize
c
c
      do 20 jj=1,6
         do 15 ii=1,natoms3
            t(ii,jj)=zero
            lp(ii,jj)=zero
            lm(ii,jj)=zero
   15    continue
   20 continue
c
      do 4 atom=1,natoms
            ix=(atom-1)*3+1
            iy=(atom-1)*3+2
            iz=(atom-1)*3+3
            lp(ix,1)=one
            lp(iy,2)=one
            lp(iz,3)=one
    4 continue
      if(debug) then
         write(iout,1010)
      endif
c
      do 5 atom=1,natoms
         lp(3*atom-1,4)=c(3,atom)
         lp(3*atom  ,4)=-c(2,atom)
    5 continue
c
      do 6 atom=1,natoms
         lp(3*atom-2,5)=c(3,atom)
         lp(3*atom  ,5)=-c(1,atom)
    6 continue
c
      do 7 atom=1,natoms
         lp(3*atom-2,6) =  c(2,atom)
         lp(3*atom-1,6) = -c(1,atom)
    7 continue
      if(debug) then
         write(iout,1020) nrt
      endif
      call schmidt(lp,t1,t2,nrt,natoms3)
c
      if(debug) then
         write(iout,1030)
      endif
c
      c(kxyz,katom)=c(kxyz,katom)-two*stpsize
c
      do 8 atom=1,natoms
            ix=(atom-1)*3+1
            iy=(atom-1)*3+2
            iz=(atom-1)*3+3
            lm(ix,1)=one
            lm(iy,2)=one
            lm(iz,3)=one
    8 continue
c
      do 9 atom=1,natoms
         lm(3*atom-1,4)=c(3,atom)
         lm(3*atom  ,4)=-c(2,atom)
    9 continue
c
      do 10 atom=1,natoms
         lm(3*atom-2,5)=c(3,atom)
         lm(3*atom  ,5)=-c(1,atom)
   10 continue
c
      do 11 atom=1,natoms
         lm(3*atom-2,6) =  c(2,atom)
         lm(3*atom-1,6) = -c(1,atom)
   11 continue
      call schmidt(lm,t1,t2,nrt,natoms3)
c
      do 12 j=1,nrt
         do 13 i=1,natoms3
            t(i,j)=(lp(i,j)-lm(i,j))/(two*stpsize)
   13    continue
   12 continue
c
      c(kxyz,katom)=c(kxyz,katom)+stpsize
c
c
      return
      end

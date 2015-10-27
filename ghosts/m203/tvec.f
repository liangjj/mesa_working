*deck %W%  %G%
      subroutine tvec(t,c,ixyz,lp,lm,natoms3,natoms,t1,t2,nrt)
      implicit integer(a-z)
      logical prnt,chkpt,singpt,cartfx,debug
      integer d2ecycl
      real*8 t(natoms3,6),lp(natoms3,6),lm(natoms3,6),c(3,natoms)
      real*8 energy,rmax,rmin,rlim,stpsize
      real*8 t1(natoms3,natoms3),t2(natoms3)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common /io/ inp,iout
      parameter (debug=.false.)
c
c
      write(iout,*) 'm203:tvec; stpsize:',stpsize
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
c      write(iout,1212) ixyz,kxyz,katom
 1212 format('  tvec: ixyz,kxyz,katom= ',3i5)
c
      c(kxyz,katom)=c(kxyz,katom)+stpsize
c
      do 20 jj=1,6
         do 15 ii=1,natoms3
            t(ii,jj)=0.0d+00
            lp(ii,jj)=0.0d+00
            lm(ii,jj)=0.0d+00
   15    continue
   20 continue
c
      do 4 atom=1,natoms
            ix=(atom-1)*3+1
            iy=(atom-1)*3+2
            iz=(atom-1)*3+3
            lp(ix,1)=1.0d+00
            lp(iy,2)=1.0d+00
            lp(iz,3)=1.0d+00
    4 continue
c      write(iout,1114)
 1114 format(' in tvec bet rot and trans ')
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
c      write(iout,1112) nrt
 1112 format(' in tvec: calling schmidt nrt=',i5)
      call schmidt(lp,t1,t2,nrt,natoms3)
c
c      write(iout,1113)
 1113 format(' in tvec: after schmidt ')
      c(kxyz,katom)=c(kxyz,katom)-2.0d+00*stpsize
c
      do 8 atom=1,natoms
            ix=(atom-1)*3+1
            iy=(atom-1)*3+2
            iz=(atom-1)*3+3
            lm(ix,1)=1.0d+00
            lm(iy,2)=1.0d+00
            lm(iz,3)=1.0d+00
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
            t(i,j)=(lp(i,j)-lm(i,j))/(2.0d+00*stpsize)
   13    continue
   12 continue
c
      c(kxyz,katom)=c(kxyz,katom)+stpsize
c
      return
      end

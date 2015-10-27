      subroutine mchead ( ndab, nklv, mrsf, mrsklv, mklvrs, ipflag,
     $     nrsf, nklvrs, maxxjk, ldafab )
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
cc
cmp   extended dummy mrsf,nklv,mrsklv,mklvrs
cc
c
      dimension mrsf(2), nklv(2), mrsklv(2), mklvrs(20,2)
c
      common /io/ inp,iout
c
c
c---------------------------------------------------
c     this program reads the header record from the
c     the j and k integral tape..the header record
c     contains directory iformation that labels each
c     block of integrals which follows
c---------------------------------------------------
c
      idisk = 1
c
cibm
c     read ( ndab'idisk ) nsym, nrsf, nklvrs, ldard, idad1,
c    1                   ida2, ( ndmy, i=1,14 ),
c    2     ( mrsf(i), i = 1, nrsf ),
c    3     ( nklv(i), i = 1, nrsf ),
c    4     ( mrsklv(i), i=1, nklvrs ),
c    5     ( ( mklvrs(i,j), i=1,20 ), j=1,nklvrs ), ldahab, ldafab
cibm
c
cccccc
c      call exrdhd(ndab,idisk,nsym,nrsf,nklvrs,ldard,idad1,ida2,
c     1 mrsf,nklv,mrsklv,mklvrs,ldahab,ldafab)
c
c     write(iout,202) ldard,idad1,ida2,ldahab,ldafab
c 202 format(' ldard idad1 ida2 ldahab ldafab ',/,5i5)
c
cccccc
c
      maxxjk=0
      do 10 i=1,nklvrs
         nrs=mklvrs(15,i)
         mrs=mklvrs(16,i)
         nkl=mklvrs(14,i)
         mkl=mklvrs(17,i)
         lxjk=nrs*mkl
cc    write(iout,11) mrs,nrs,mkl,nkl,lxjk
cc 11 format('  mchead  mrs nrs mkl nkl lxjk ',5i7)
         maxxjk=max0(maxxjk,lxjk)
 10   continue
c
c
      return
      end

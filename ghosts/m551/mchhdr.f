*deck %W%  %G%
      subroutine mchhdr ( ndab, nklv, mrsf, mrsklv, mklvrs, ipflag,
     $     nrsf, nklvrs, ldahab,ldafab)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy mrsf,nklv,mrsklv,mklvrs
cc
      dimension mrsf(2), nklv(2), mrsklv(2), mklvrs(20,2)
      dimension ihd(20)
c
c---------------------------------------------------
c     this program reads the header record from the
c     the j and k integral tape..the header record
c     contains directory iformation that labels each
c     block of integrals which follows
c---------------------------------------------------
c
      idisk = 1
cibm
c     read ( ndab'idisk ) nsym, nrsf, nklvrs, ldard, idad1,
c    1                   ida2, ( ndmy, i=1,14 ),
c    2     ( mrsf(i), i = 1, nrsf ),
c    3     ( nklv(i), i = 1, nrsf ),
c    4     ( mrsklv(i), i=1, nklvrs ),
c    5     ( ( mklvrs(i,j), i=1,20 ), j=1,nklvrs ), ldahab,ldafab
cibm
c
      call lnkerr('exrdhd')
cps      call exrdhd(ndab,idisk,nsym,nrsf,nklvrs,ldard,idad1,ida2,
cps     $     mrsf,nklv,mrsklv,mklvrs,ldahab,ldafab)
c
      if(ipflag.ne.0) go to 6000
c
      return
c
 6000 continue
c
c
c     print half transformed symmetry block infromation.
c
 7000 write ( 6, 7100 ) ( j, (mklvrs(i,j),i=1,15), j=1,nklvrs )
      write ( 6, 7110 ) ( j, (mklvrs(i,j),i=16,20), j=1,nklvrs )
 7100 format ('1half transformed symmetry block orbital infromation'//
     $     5x,'iklvrs  mk  ml  mr  ms mkl mrs  mi indx type  nok  nol',
     $     '  nbr  nbs  nkl  nrs'/(5x,i6,7i4,8i5))
 7110 format ('0half transformed symmetry block data set blocking  '//
     $     5x,'iklvrs     mrs     mkl   idad0   idad1   idad2'/
     $     (5x,i6,5i8))
c
 7130 ni = 0
      ib = 0
      write ( 6,7132 )
 7132 format ('0symetry blocks sorted according to ao pairs'//
     $     5x,'    i mr ms nklv mrsf sym. bloks with the same mr ',
     $     'ms pair')
      do 7150 mr = 1, nsym
         do 7150 ms = 1, mr
            ni = ni + 1
            mkl = nklv(ni)
            if ( mkl .eq. 0 ) go to 7150
            ia = ib + 1
            ib = ib + mkl
            write ( 6, 7140 ) ni, mr, ms, mkl, mrsf(ni),
     $           ( mrsklv(i), i=ia,ib )
 7140       format (5x,i5,2i3,2i5,10i5/(26x,10i5))
 7150    continue
c
c
         return
         end

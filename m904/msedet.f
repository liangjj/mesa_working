*deck @(#)msedet.f	5.1  11/6/94
      subroutine msedet(iopn,jdbl,ms,ifiga,ifigb)
c
c  generate all of the determinants with ms even needed to make
c  double-group-adapted functions from a spatial configuration.
c
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c2/ msbas,jsym,nopen,kdbl,ndeti,nsefi,ncpi
c
      dimension iopn(*),jdbl(*),ms(*),ifiga(mxocc,*),ifigb(mxocc,*)
c
      nalpha = nel/2
      kbrch = nopen/2 + 1
      go to (100,120,150,210,280), kbrch
c
c  case 1:  no open-shell electrons
c
c          number of double-group functions = 1
c          number of determinants           = 1
c
c  the list of doubly occupied orbitals goes in both the spin-up
c  and spin-down lists
c
  100 ms(1) = msbas
      do 110 jj = 1,kdbl
      ifiga(jj,1) = jdbl(jj)
  110 ifigb(jj,1) = jdbl(jj)
c
      go to 500
c
c  case 2:  two open-shell electrons
c
c          number of double-group functions = 1
c          number of determinants           = 2
c
c  arrangement of the singly occupied orbitals
c
c         singly occupied orbitals
c     det  alpha list  beta list   (sign) spin function
c
c      1        1          2        (+) a b
c      2        2          1        (-) b a
c
c  the spin function and its sign are as obtained when the singly
c  occupied orbitals are arranged in ascending order
c
c  the doubly occupied orbitals are given first in the lists
c
  120 ms(1) = msbas
      ms(2) = msbas
c
      jdbl(nalpha) = iopn(1)
      do 130 jj = 1,nalpha
      ifiga(jj,1) = jdbl(jj)
  130 ifigb(jj,2) = jdbl(jj)
c
      jdbl(nalpha) = iopn(2)
      do 140 jj = 1,nalpha
      ifigb(jj,1) = jdbl(jj)
  140 ifiga(jj,2) = jdbl(jj)
c
      go to 500
c
c  case 3:  four open-shell electrons
c
c          number of double-group functions = 4
c          number of determinants           = 8
c
c  arrangement of the singly occupied orbitals
c
c         singly occupied orbitals
c     det      alpha list    beta list    (sign) spin function
c
c      1       1  2          3  4          (+) a a b b
c      2       3  4          1  2          (+) b b a a
c      3       1  3          2  4          (-) a b a b
c      4       2  4          1  3          (-) b a b a
c      5       2  3          1  4          (+) b a a b
c      6       1  4          2  3          (+) a b b a
c      7       1  2  3  4                  (+) a a a a
c      8                     1  2  3  4    (+) b b b b
c
c  determinants with ms = 0
c
  150 ndeta = 1
      do 170 i2=2,4
      jdbl(nalpha) = iopn(i2)
      do 170 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      ndetb = min(ndeta,13-ndeta)
      ndetc = min(ndeta+1,12-ndeta)
      ms(ndetb) = msbas
      do 160 jj=1,nalpha
      ifiga(jj,ndetb) = jdbl(jj)
  160 ifigb(jj,ndetc) = jdbl(jj)
  170 ndeta = ndeta + 2
c
c  determinants with ms = 2,-2
c
      ms(7) = msbas + 2
      ms(8) = msbas - 2
      do 180 jj = 1,kdbl
      ifigb(jj,7) = jdbl(jj)
  180 ifiga(jj,8) = jdbl(jj)
      do 190 i1=1,4
  190 jdbl(kdbl+i1) = iopn(i1)
      nalpha = nalpha + 2
      do 200 jj=1,nalpha
      ifiga(jj,7) = jdbl(jj)
  200 ifigb(jj,8) = jdbl(jj)
c
      go to 500
c
c  case 4:  six open-shell electrons
c
c          number of double-group functions = 16
c          number of determinants           = 32
c
c  arrangement of the singly occupied orbitals
c
c         singly occupied orbitals
c     det       alpha list        beta list     (sign) spin function
c
c      1       1  2  3          4  5  6          (+) a a a b b b
c      2       4  5  6          1  2  3          (-) b b b a a a
c      3       1  2  4          3  5  6          (-) a a b a b b
c      4       3  5  6          1  2  4          (+) b b a b a a
c      5       1  3  4          2  5  6          (+) a b a a b b
c      6       2  5  6          1  3  4          (-) b a b b a a
c      7       2  3  4          1  5  6          (-) b a a a b b
c      8       1  5  6          2  3  4          (+) a b b b a a
c      9       1  2  5          3  4  6          (+) a a b b a b
c     10       3  4  6          1  2  5          (-) b b a a b a
c     11       1  3  5          2  4  6          (-) a b a b a b
c     12       2  4  6          1  3  5          (+) b a b a b a
c     13       2  3  5          1  4  6          (+) b a a b a b
c     14       1  4  6          2  3  5          (-) a b b a b a
c     15       1  4  5          2  3  6          (+) a b b a a b
c     16       2  3  6          1  4  5          (-) b a a b b a
c     17       2  4  5          1  3  6          (-) b a b a a b
c     18       1  3  6          2  4  5          (+) a b a b b a
c     19       3  4  5          1  2  6          (+) b b a a a b
c     20       1  2  6          3  4  5          (-) a a b b b a
c     21       1  2  3  4  5    6                (+) a a a a a b
c     22       6                1  2  3  4  5    (-) b b b b b a
c     23       1  2  3  4  6    5                (-) a a a a b a
c     24       5                1  2  3  4  6    (+) b b b b a b
c     25       1  2  3  5  6    4                (+) a a a b a a
c     26       4                1  2  3  5  6    (-) b b b a b b
c     27       1  2  4  5  6    3                (-) a a b a a a
c     28       3                1  2  3  5  6    (+) b b a b b b
c     29       1  3  4  5  6    2                (+) a b a a a a
c     30       2                1  3  4  5  6    (-) b a b b b b
c     31       2  3  4  5  6    1                (-) b a a a a a
c     32       1                2  3  4  5  6    (+) a b b b b b
c
c  determinants with ms = 0
c
  210 ndeta = 1
      do 230 i3=3,6
      jdbl(nalpha) = iopn(i3)
      do 230 i2=2,i3-1
      jdbl(nalpha-1) = iopn(i2)
      do 230 i1=1,i2-1
      jdbl(nalpha-2) = iopn(i1)
      ndetb = min(ndeta,41-ndeta)
      ndetc = min(ndeta+1,40-ndeta)
      ms(ndetb) = msbas
      do 220 jj=1,nalpha
      ifiga(jj,ndetb) = jdbl(jj)
  220 ifigb(jj,ndetc) = jdbl(jj)
  230 ndeta = ndeta + 2
c
c  determinants with ms = 2,-2
c
      nalpha = nalpha + 2
      ndeta = 21
      do 250 i5=5,6
      jdbl(nalpha) = iopn(i5)
      do 250 i4=4,i5-1
      jdbl(nalpha-1) = iopn(i4)
      do 250 i3=3,i4-1
      jdbl(nalpha-2) = iopn(i3)
      do 250 i2=2,i3-1
      jdbl(nalpha-3) = iopn(i2)
      do 250 i1=1,i2-1
      jdbl(nalpha-4) = iopn(i1)
      ms(ndeta) = (msbas + 2)
      ms(ndeta+1) = (msbas - 2)
      do 240 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  240 ifigb(jj,ndeta+1) = jdbl(jj)
  250 ndeta = ndeta + 2
c
      nalpha = nalpha - 4
      ndeta = 32
      do 270 i1=1,6
      jdbl(nalpha) = iopn(i1)
      do 260 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  260 ifiga(jj,ndeta) = jdbl(jj)
  270 ndeta = ndeta - 2
c
      go to 500
c
c  case 5:  eight open-shell electrons
c
c          number of double-group functions =  64
c          number of determinants           = 128
c
c  arrangement of the singly occupied orbitals is analogous to
c  that for case 4.
c
c  determinants with ms = 0
c
  280 ndeta = 1
      do 300 i4=4,8
      jdbl(nalpha) = iopn(i4)
      do 300 i3=3,i4-1
      jdbl(nalpha-1) = iopn(i3)
      do 300 i2=2,i3-1
      jdbl(nalpha-2) = iopn(i2)
      do 300 i1=1,i2-1
      jdbl(nalpha-3) = iopn(i1)
      ndetb = min(ndeta,141-ndeta)
      ndetc = min(ndeta+1,140-ndeta)
      ms(ndetb) = msbas
      do 290 jj=1,nalpha
      ifiga(jj,ndetb) = jdbl(jj)
  290 ifigb(jj,ndetc) = jdbl(jj)
  300 ndeta = ndeta + 2
c
c  determinants with ms = 2,-2
c
      nalpha = nalpha + 2
      ndeta = 71
      do 320 i6=6,8
      jdbl(nalpha) = iopn(i6)
      do 320 i5=5,i6-1
      jdbl(nalpha-1) = iopn(i5)
      do 320 i4=4,i5-1
      jdbl(nalpha-2) = iopn(i4)
      do 320 i3=3,i4-1
      jdbl(nalpha-3) = iopn(i3)
      do 320 i2=2,i3-1
      jdbl(nalpha-4) = iopn(i2)
      do 320 i1=1,i2-1
      jdbl(nalpha-5) = iopn(i1)
      ms(ndeta) = (msbas + 2)
      ms(ndeta+1) = (msbas - 2)
      do 310 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  310 ifigb(jj,ndeta+1) = jdbl(jj)
  320 ndeta = ndeta + 2
c
      nalpha = nalpha - 4
      ndeta = 126
      do 340 i2=2,8
      jdbl(nalpha) = iopn(i2)
      do 340 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      do 330 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  330 ifiga(jj,ndeta) = jdbl(jj)
  340 ndeta = ndeta - 2
c
c  determinants with ms = 4,-4
c
      ms(127) = msbas + 4
      ms(128) = msbas - 4
      do 350 i1=1,8
  350 jdbl(kdbl+i1) = iopn(i1)
      do 360 jj=1,kdbl
      ifigb(jj,127) = jdbl(jj)
  360 ifiga(jj,128) = jdbl(jj)
      nalpha = nalpha + 6
      do 370 jj=1,nalpha
      ifiga(jj,127) = jdbl(jj)
  370 ifigb(jj,128) = jdbl(jj)
c
  500 return
c
      end

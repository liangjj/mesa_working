*deck @(#)msodet.f	5.1  11/6/94
      subroutine msodet(iopn,jdbl,ms,ifiga,ifigb)
c
c  generate all of the determinants with ms odd needed to make
c  double-group-adapted functions from a spatial configuration.
c
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c2/ msbas,jsym,nopen,kdbl,ndeti,nsefi,ncpi
c
      dimension iopn(*),jdbl(*),ms(*),ifiga(mxocc,*),ifigb(mxocc,*)
c
      nalpha = nel/2 + 1
      kbrch = nopen/2
      go to (100,130,180,260), kbrch
c
c  case 1:  two open-shell electrons
c
c          number of double-group functions = 1
c          number of determinants           = 2
c
c  arrangement of the singly occupied orbitals
c
c          singly occupied orbitals
c     det   alpha list   beta list    (sign) spin function
c
c      1       1  2                    (+) a a
c      2                    1  2       (+) b b
c
c  the list of doubly occupied orbitals goes in both the spin-up
c  and spin-down lists.
c
  100 ms(1) = msbas + 1
      ms(2) = msbas - 1
      jdbl(kdbl+1) = iopn(1)
      jdbl(kdbl+2) = iopn(2)
      do 110 jj = 1,kdbl
      ifigb(jj,1) = jdbl(jj)
  110 ifiga(jj,2) = jdbl(jj)
      do 120 jj = 1,nalpha
      ifiga(jj,1) = jdbl(jj)
  120 ifigb(jj,2) = jdbl(jj)
c
      go to 500
c
c  case 2:  four open-shell electrons
c
c          number of double-group functions = 4
c          number of determinants           = 8
c
c  arrangement of the singly occupied orbitals
c
c          singly occupied orbitals
c     det   alpha list   beta list    (sign) spin function
c
c      1      1  2  3     4            (+) a a a b
c      2      4           1  2  3      (-) b b b a
c      3      1  2  4     3            (-) a a b a
c      4      3           1  2  4      (+) b b a b
c      5      1  3  4     2            (+) a b a a
c      6      2           1  3  4      (-) b a b b
c      7      2  3  4     1            (-) b a a a
c      8      1           2  3  4      (+) a b b b
c
  130 ndeta = 1
      do 150 i3=3,4
      jdbl(nalpha) = iopn(i3)
      do 150 i2=2,i3-1
      jdbl(nalpha-1) = iopn(i2)
      do 150 i1=1,i2-1
      jdbl(nalpha-2) = iopn(i1)
      ms(ndeta) = (msbas + 1)
      ms(ndeta+1) = (msbas - 1)
      do 140 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  140 ifigb(jj,ndeta+1) = jdbl(jj)
  150 ndeta = ndeta + 2
c
      nalpha = nalpha - 2
      ndeta = 8
      do 170 i1=1,4
      jdbl(nalpha) = iopn(i1)
      do 160 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  160 ifiga(jj,ndeta) = jdbl(jj)
  170 ndeta = ndeta - 2
c
      go to 500
c
c  case 3:  six open-shell electrons
c
c          number of double-group functions = 16
c          number of determinants           = 32
c
c  arrangement of the singly occupied orbitals
c
c             singly occupied orbitals
c     det       alpha list        beta list     (sign) spin function
c
c      1     1  2  3  4       5  6               (+) a a a a b b
c      2     5  6             1  2  3  4         (+) b b a a a a
c      3     1  2  3  5       4  6               (-) a a a b a b
c      4     4  6             1  2  3  5         (-) b b b a b a
c      5     1  2  4  5       3  6               (+) a a b a a b
c      6     3  6             1  2  4  5         (+) b b a b b a
c      7     1  3  4  5       2  6               (-) a b a a a b
c      8     2  6             1  3  4  5         (-) b a b b b a
c      9     2  3  4  5       1  6               (+) b a a a a b
c     10     1  6             2  3  4  5         (+) a b b b b a
c     11     1  2  3  6       4  5               (+) a a a b b a
c     12     4  5             1  2  3  6         (+) b b b a a b
c     13     1  2  4  6       3  5               (-) a a b a b a
c     14     3  5             1  2  4  6         (-) b b a b a b
c     15     1  3  4  6       2  5               (+) a b a a b a
c     16     2  5             1  3  4  6         (+) b a b b a b
c     17     2  3  4  6       1  5               (-) b a a a b a
c     18     1  5             2  3  4  6         (-) a b b b a b
c     19     1  2  5  6       3  4               (+) a a b b a a
c     20     3  4             1  2  5  6         (+) b b a a b b
c     21     1  3  5  6       2  4               (-) a b a b a a
c     22     2  4             1  3  5  6         (-) b a b a b b
c     23     2  3  5  6       1  4               (+) b a a b a a
c     24     1  4             2  3  5  6         (+) a b b a b b
c     25     1  4  5  6       2  3               (+) a b b a a a
c     26     2  3             1  4  5  6         (+) b a a b b b
c     27     2  4  5  6       1  3               (-) b a b a a a
c     28     1  3             2  4  5  6         (-) a b a b b b
c     29     3  4  5  6       1  2               (+) b b a a a a
c     30     1  2             3  4  5  6         (+) a a b b b b
c     31     1  2  3  4  5  6                    (+) a a a a a a
c     32                      1  2  3  4  5  6   (+) b b b b b b
c
c  determinants with ms = 1,-1
c
  180 ndeta = 1
      do 200 i4 = 4,6
      jdbl(nalpha) = iopn(i4)
      do 200 i3=3,i4-1
      jdbl(nalpha-1) = iopn(i3)
      do 200 i2=2,i3-1
      jdbl(nalpha-2) = iopn(i2)
      do 200 i1=1,i2-1
      jdbl(nalpha-3) = iopn(i1)
      ms(ndeta) = (msbas + 1)
      ms(ndeta+1) = (msbas - 1)
      do 190 jj = 1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  190 ifigb(jj,ndeta+1) = jdbl(jj)
  200 ndeta = ndeta + 2
c
      nalpha = nalpha - 2
      ndeta = 30
      do 220 i2 = 2,6
      jdbl(nalpha) = iopn(i2)
      do 220 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      do 210 jj = 1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  210 ifiga(jj,ndeta) = jdbl(jj)
  220 ndeta = ndeta - 2
c
c  determinants with ms = 3,-3
c
      ms(31) = msbas + 3
      ms(32) = msbas - 3
      do 230 jj = 1,kdbl
      ifigb(jj,31) = jdbl(jj)
  230 ifiga(jj,32) = jdbl(jj)
      do 240 i1=1,6
  240 jdbl(kdbl+i1) = iopn(i1)
      nalpha = nalpha + 4
      do 250 jj=1,nalpha
      ifiga(jj,31) = jdbl(jj)
  250 ifigb(jj,32) = jdbl(jj)
c
      go to 500
c
c  case 4:  eight open-shell electrons
c
c          number of double-group functions =  64
c          number of determinants           = 128
c
c  arrangement of the singly occupied orbitals is analogous to
c  that for case 3.
c
c  determinants with ms = 1,-1
c
  260 ndeta = 1
      do 280 i5=5,8
      jdbl(nalpha) = iopn(i5)
      do 280 i4=4,i5-1
      jdbl(nalpha-1) = iopn(i4)
      do 280 i3=3,i4-1
      jdbl(nalpha-2) = iopn(i3)
      do 280 i2=2,i3-1
      jdbl(nalpha-3) = iopn(i2)
      do 280 i1=1,i2-1
      jdbl(nalpha-4) = iopn(i1)
      ms(ndeta) = (msbas + 1)
      ms(ndeta+1) = (msbas - 1)
      do 270 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  270 ifigb(jj,ndeta+1) = jdbl(jj)
  280 ndeta = ndeta + 2
c
      nalpha = nalpha - 2
      ndeta = 112
      do 300 i3=3,8
      jdbl(nalpha) = iopn(i3)
      do 300 i2=2,i3-1
      jdbl(nalpha-1) = iopn(i2)
      do 300 i1=1,i2-1
      jdbl(nalpha-2) = iopn(i1)
      do 290 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  290 ifiga(jj,ndeta) = jdbl(jj)
  300 ndeta = ndeta - 2
c
c  determinants with ms = 3,-3
c
      nalpha = nalpha + 4
      ndeta = 113
      do 320 i7=7,8
      jdbl(nalpha) = iopn(i7)
      do 320 i6=6,i7-1
      jdbl(nalpha-1) = iopn(i6)
      do 320 i5=5,i6-1
      jdbl(nalpha-2) = iopn(i5)
      do 320 i4=4,i5-1
      jdbl(nalpha-3) = iopn(i4)
      do 320 i3=3,i4-1
      jdbl(nalpha-4) = iopn(i3)
      do 320 i2=2,i3-1
      jdbl(nalpha-5) = iopn(i2)
      do 320 i1=1,i2-1
      jdbl(nalpha-6) = iopn(i1)
      ms(ndeta) = (msbas + 3)
      ms(ndeta+1) = (msbas - 3)
      do 310 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  310 ifigb(jj,ndeta+1) = jdbl(jj)
  320 ndeta = ndeta + 2
c
      nalpha = nalpha - 6
      ndeta = 128
      do 340 i1=1,8
      jdbl(nalpha) = iopn(i1)
      do 330 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  330 ifiga(jj,ndeta) = jdbl(jj)
  340 ndeta = ndeta - 2
c
  500 return
c
      end

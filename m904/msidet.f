*deck @(#)msidet.f	5.1  11/6/94
      subroutine msidet(iopn,jdbl,ms,ifiga,ifigb)
c
c  generate all of the determinants with ms integral needed to make
c  double-group-adapted functions from a spatial configuration.
c
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c2/ msbas,jsym,nopen,kdbl,ndeti,nsefi,ncpi
c
      dimension iopn(*),jdbl(*),ms(*),ifiga(mxocc,*),ifigb(mxocc,*)
c
      nalpha = nel/2
      kbrch = nopen/2 + 1
      go to (100,120,170,270,410), kbrch
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
      go to 600
c
c  case 2:  two open-shell electrons
c
c          number of double-group functions = 4 
c          number of determinants           = 4 
c
c  arrangement of the singly occupied orbitals
c
c         singly occupied orbitals
c     det  alpha list  beta list   (sign) spin function
c
c      1        1          2        (+) a b
c      2        2          1        (-) b a
c      3        1  2                (+) a a
c      4                   1  2     (+) b b
c
c  the spin function and its sign are as obtained when the singly
c  occupied orbitals are arranged in ascending order
c
c  the doubly occupied orbitals are given first in the lists
c
c  determinants with ms = 0
c
  120 ms(1) = msbas
      ms(2) = msbas
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
c     determinants with ms = 1,-1
c
      ms(3) = msbas + 1
      ms(4) = msbas - 1
      do 150 jj = 1,kdbl
      ifigb(jj,3) = jdbl(jj)
  150 ifiga(jj,4) = jdbl(jj)
      jdbl(kdbl+1) = iopn(1)
      jdbl(kdbl+2) = iopn(2)
      nalpha = nalpha + 1
      do 160 jj=1,nalpha
      ifiga(jj,3) = jdbl(jj)
  160 ifigb(jj,4) = jdbl(jj)
      go to 600
c
c  case 3:  four open-shell electrons
c
c          number of double-group functions = 16 
c          number of determinants           = 16 
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
c      7       1  2  3       4             (+) a a a b
c      8       4             1  2  3       (-) b b b a
c      9       1  2  4       3             (-) a a b a
c     10       3             1  2  4       (+) b b a b
c     11       1  3  4       2             (+) a b a a
c     12       2             1  3  4       (-) b a b b
c     13       2  3  4       1             (-) b a a a
c     14       1             2  3  4       (+) a a a b
c     15       1  2  3  4                  (+) a a a a
c     16                     1  2  3  4    (+) b b b b
c
c  determinants with ms = 0
c
  170 ndeta = 1
      do 190 i2=2,4
      jdbl(nalpha) = iopn(i2)
      do 190 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      ndetb = min(ndeta,13-ndeta)
      ndetc = min(ndeta+1,12-ndeta)
      ms(ndetb) = msbas
      do 180 jj=1,nalpha
      ifiga(jj,ndetb) = jdbl(jj)
  180 ifigb(jj,ndetc) = jdbl(jj)
  190 ndeta = ndeta + 2
c
c  determinants with ms = 1,-1
c
      nalpha = nalpha + 1
      ndeta = 7
      do 210 i3=3,4
      jdbl(nalpha) = iopn(i3)
      do 210 i2=2,i3-1
      jdbl(nalpha-1) = iopn(i2)
      do 210 i1=1,i2-1
      jdbl(nalpha-2) = iopn(i1)
      ms(ndeta) = (msbas+1)
      ms(ndeta+1) = (msbas-1)
      do 200 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  200 ifigb(jj,ndeta+1) = jdbl(jj)
  210 ndeta = ndeta + 2
c
      nalpha = nalpha - 2
      ndeta = 14
      do 230 i1=1,4
      jdbl(nalpha) = iopn(i1)
      do 220 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  220 ifiga(jj,ndeta) = jdbl(jj)
  230 ndeta = ndeta - 2
c
c  determinants with ms = 2,-2
c
      ms(15) = msbas + 2
      ms(16) = msbas - 2
      do 240 jj = 1,kdbl
      ifigb(jj,15) = jdbl(jj)
  240 ifiga(jj,16) = jdbl(jj)
      do 250 i1=1,4
  250 jdbl(kdbl+i1) = iopn(i1)
      nalpha = nalpha +3 
      do 260 jj=1,nalpha
      ifiga(jj,15) = jdbl(jj)
  260 ifigb(jj,16) = jdbl(jj)
c
      go to 600
c
c  case 4:  six open-shell electrons
c
c          number of double-group functions = 64
c          number of determinants           = 64
c
c  arrangement of the singly occupied orbitals
c
c             singly occupied orbitals
c     det       alpha list        beta list     (sign) spin function
c
c      1     1  2  3          4  5  6            (+) a a a b b b
c      2     4  5  6          1  2  3            (-) b b b a a a
c      3     1  2  4          3  5  6            (-) a a b a b b
c      4     3  5  6          1  2  4            (+) b b a b a a
c      5     1  3  4          2  5  6            (+) a b a a b b
c      6     2  5  6          1  3  4            (-) b a b b a a
c      7     2  3  4          1  5  6            (-) b a a a b b
c      8     1  5  6          2  3  4            (+) a b b b a a
c      9     1  2  5          3  4  6            (+) a a b b a b
c     10     3  4  6          1  2  5            (-) b b a a b a
c     11     1  3  5          2  4  6            (-) a b a b a b
c     12     2  4  6          1  3  5            (+) b a b a b a
c     13     2  3  5          1  4  6            (+) b a a b a b
c     14     1  4  6          2  3  5            (-) a b b a b a
c     15     1  4  5          2  3  6            (+) a b b a a b
c     16     2  3  6          1  4  5            (-) b a a b b a
c     17     2  4  5          1  3  6            (-) b a b a a b
c     18     1  3  6          2  4  5            (+) a b a b b a
c     19     3  4  5          1  2  6            (+) b b a a a b
c     20     1  2  6          3  4  5            (-) a a b b b a
c     21     1  2  3  4       5  6               (+) a a a a b b
c     22     5  6             1  2  3  4         (+) b b a a a a
c     23     1  2  3  5       4  6               (-) a a a b a b
c     24     4  6             1  2  3  5         (-) b b b a b a
c     25     1  2  4  5       3  6               (+) a a b a a b
c     26     3  6             1  2  4  5         (+) b b a b b a
c     27     1  3  4  5       2  6               (-) a b a a a b
c     28     2  6             1  3  4  5         (-) b a b b b a
c     29     2  3  4  5       1  6               (+) b a a a a b
c     30     1  6             2  3  4  5         (+) a b b b b a
c     31     1  2  3  6       4  5               (+) a a a b b a
c     32     4  5             2  3  4  5         (+) b b b a a b
c     33     1  2  4  6       3  5               (-) a a b a b a
c     34     3  5             1  2  4  6         (-) b b a b a b
c     35     1  3  4  6       2  5               (+) a b a a b a
c     36     2  5             1  3  4  6         (+) b a b b a b
c     37     2  3  4  6       1  5               (-) b a a a b a
c     38     1  5             2  3  4  6         (-) a b b b a b
c     39     1  2  5  6       3  4               (+) a a b b a a
c     40     3  4             1  2  5  6         (+) b b a a b b
c     41     1  3  5  6       2  4               (-) a b a b a a
c     42     2  4             1  3  5  6         (-) b a b a b b
c     43     2  3  5  6       1  4               (+) b a a b a a
c     44     1  4             2  3  5  6         (+) a b b a b b
c     45     1  4  5  6       2  3               (+) a b b a a a
c     46     2  3             1  4  5  6         (+) b a a b b b
c     47     2  4  5  6       1  3               (-) b a b a a a
c     48     1  3             2  4  5  6         (-) a b a b b b
c     49     3  4  5  6       1  2               (+) b b a a a a
c     50     1  2             3  4  5  6         (+) a a b b b b
c     51     1  2  3  4  5    6                  (+) a a a a a b
c     52     6                1  2  3  4  5      (-) b b b b b a
c     53     1  2  3  4  6    5                  (-) a a a a b a
c     54     5                1  2  3  4  6      (+) b b b b a b
c     55     1  2  3  5  6    4                  (+) a a a b a a
c     56     4                1  2  3  5  6      (-) b b b a b b
c     57     1  2  4  5  6    3                  (-) a a b a a a
c     58     3                1  2  3  5  6      (+) b b a b b b
c     59     1  3  4  5  6    2                  (+) a b a a a a
c     60     2                1  3  4  5  6      (-) b a b b b b
c     61     2  3  4  5  6    1                  (-) b a a a a a
c     62     1                2  3  4  5  6      (+) a b b b b b
c     63     1  2  3  4  5  6                    (+) a a a a a a
c     64                      1  2  3  4  5  6   (+) b b b b b b
c
c  determinants with ms = 0
c
  270 ndeta = 1
      do 290 i3=3,6
      jdbl(nalpha) = iopn(i3)
      do 290 i2=2,i3-1
      jdbl(nalpha-1) = iopn(i2)
      do 290 i1=1,i2-1
      jdbl(nalpha-2) = iopn(i1)
      ndetb = min(ndeta,41-ndeta)
      ndetc = min(ndeta+1,40-ndeta)
      ms(ndetb) = msbas
      do 280 jj=1,nalpha
      ifiga(jj,ndetb) = jdbl(jj)
  280 ifigb(jj,ndetc) = jdbl(jj)
  290 ndeta = ndeta + 2
c
c  determinants with ms = 1,-1
c
      nalpha = nalpha + 1
      ndeta = 1
      do 310 i4=4,6
      jdbl(nalpha) = iopn(i4)
      do 310 i3=3,i4-1
      jdbl(nalpha-1) = iopn(i3)
      do 310 i2=2,i3-1
      jdbl(nalpha-2) = iopn(i2)
      do 310 i1=1,i2-1
      jdbl(nalpha-3) = iopn(i1)
      ms(ndeta) = (msbas + 1)
      ms(ndeta+1) = (msbas - 1)
      do 300 jj = 1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  300 ifigb(jj,ndeta+1) = jdbl(jj)
  310 ndeta = ndeta + 2
c
      nalpha = nalpha - 2
      ndeta = 50
      do 330 i2=2,6
      jdbl(nalpha) = iopn(i2)
      do 330 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      do 320 jj = 1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  320 ifiga(jj,ndeta) = jdbl(jj)
  330 ndeta = ndeta - 2
c
c  determinants with ms = 2,-2
c
      nalpha = nalpha + 3 
      ndeta = 51
      do 350 i5=5,6
      jdbl(nalpha) = iopn(i5)
      do 350 i4=4,i5-1
      jdbl(nalpha-1) = iopn(i4)
      do 350 i3=3,i4-1
      jdbl(nalpha-2) = iopn(i3)
      do 350 i2=2,i3-1
      jdbl(nalpha-3) = iopn(i2)
      do 350 i1=1,i2-1
      jdbl(nalpha-4) = iopn(i1)
      ms(ndeta) = (msbas + 2)
      ms(ndeta+1) = (msbas - 2)
      do 340 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  340 ifigb(jj,ndeta+1) = jdbl(jj)
  350 ndeta = ndeta + 2
c
      nalpha = nalpha - 4
      ndeta = 62
      do 370 i1=1,6
      jdbl(nalpha) = iopn(i1)
      do 360 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  360 ifiga(jj,ndeta) = jdbl(jj)
  370 ndeta = ndeta - 2
c
c  determinants with ms = 3,-3
c
      ms(63) = msbas + 3
      ms(64) = msbas - 3
      do 380 jj = 1,kdbl
      ifigb(jj,63) = jdbl(jj)
  380 ifiga(jj,64) = jdbl(jj)
      do 390 i1=1,6
  390 jdbl(kdbl+i1) = iopn(i1)
      nalpha = nalpha + 5 
      do 400 jj=1,nalpha
      ifiga(jj,63) = jdbl(jj)
  400 ifigb(jj,64) = jdbl(jj)
c
      go to 600
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
  410 ndeta = 1
      do 430 i4=4,8
      jdbl(nalpha) = iopn(i4)
      do 430 i3=3,i4-1
      jdbl(nalpha-1) = iopn(i3)
      do 430 i2=2,i3-1
      jdbl(nalpha-2) = iopn(i2)
      do 430 i1=1,i2-1
      jdbl(nalpha-3) = iopn(i1)
      ndetb = min(ndeta,141-ndeta)
      ndetc = min(ndeta+1,140-ndeta)
      ms(ndetb) = msbas
      do 420 jj=1,nalpha
      ifiga(jj,ndetb) = jdbl(jj)
  420 ifigb(jj,ndetc) = jdbl(jj)
  430 ndeta = ndeta + 2
c
c  determinants with ms = 2,-2
c
      nalpha = nalpha + 2
      ndeta = 71
      do 490 i6=6,8
      jdbl(nalpha) = iopn(i6)
      do 490 i5=5,i6-1
      jdbl(nalpha-1) = iopn(i5)
      do 490 i4=4,i5-1
      jdbl(nalpha-2) = iopn(i4)
      do 490 i3=3,i4-1
      jdbl(nalpha-3) = iopn(i3)
      do 490 i2=2,i3-1
      jdbl(nalpha-4) = iopn(i2)
      do 490 i1=1,i2-1
      jdbl(nalpha-5) = iopn(i1)
      ms(ndeta) = (msbas + 2)
      ms(ndeta+1) = (msbas - 2)
      do 480 jj=1,nalpha
      ifiga(jj,ndeta) = jdbl(jj)
  480 ifigb(jj,ndeta+1) = jdbl(jj)
  490 ndeta = ndeta + 2
c
      nalpha = nalpha - 4
      ndeta = 126
      do 510 i2=2,8
      jdbl(nalpha) = iopn(i2)
      do 510 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      do 500 jj=1,nalpha
      ifigb(jj,ndeta-1) = jdbl(jj)
  500 ifiga(jj,ndeta) = jdbl(jj)
  510 ndeta = ndeta - 2
c
c  determinants with ms = 4,-4
c
      ms(127) = msbas + 4
      ms(128) = msbas - 4
      do 560 i1=1,8
  560 jdbl(kdbl+i1) = iopn(i1)
      do 570 jj=1,kdbl
      ifigb(jj,127) = jdbl(jj)
  570 ifiga(jj,128) = jdbl(jj)
      nalpha = nalpha + 6
      do 580 jj=1,nalpha
      ifiga(jj,127) = jdbl(jj)
  580 ifigb(jj,128) = jdbl(jj)
c
  600 return
c
      end

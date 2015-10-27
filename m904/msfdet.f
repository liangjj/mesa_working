*deck @(#)msfdet.f	5.1  11/6/94
      subroutine msfdet(iopn,jdbl,ms,ifiga,ifigb,isgn)
c
c  generate all of the determinants with ms fractional needed to make
c  double-group-adapted functions from a spatial configuration.
c
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c2/ msbas,jsym,nopen,kdbl,ndeti,nsefi,ncpi
c
      dimension iopn(*),jdbl(*),ms(*),ifiga(mxocc,*),ifigb(mxocc,*)
c
      nalpha = (nel + 1)/2
      kbrch = (nopen + 1)/2
      go to (100,130,210,330), kbrch
c
c  case 1:  one open-shell electron
c
c          number of double-group functions = 1
c          number of determinants           = 1
c
c  the list of doubly occupied orbitals goes in both the spin-up
c  and spin-down lists.  the singly occupied orbital goes in the
c  spin-up list.
c
  100 ms(1) = msbas + (1 + isgn)/2
      do 110 jj = 1,kdbl
  110 ifigb(jj,1) = jdbl(jj)
c
      jdbl(nalpha) = iopn(1)
      do 120 jj = 1,nalpha
  120 ifiga(jj,1) = jdbl(jj)
c
      go to 500
c
c  case 2:  three open-shell electrons
c
c          number of double-group functions = 4
c          number of determinants           = 4
c
c  arrangement of the singly occupied orbitals
c
c          singly occupied orbitals
c     det   alpha list   beta list   (sign) spin function
c
c      1       1  2       3           (+) a a b
c      2       1  3       2           (-) a b a
c      3       2  3       1           (+) b a a
c      4                  1  2  3     (+) b b b
c
c  (alpha), then (beta), lists for determinants with ms = (1/2)
c
  130 ndeta = 1
      do 150 i2=2,3
      jdbl(nalpha) = iopn(i2)
      do 150 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      ms(ndeta) = (msbas + (1 + isgn)/2)
      do 140 jj=1,nalpha
  140 ifiga(jj,ndeta) = jdbl(jj)
  150 ndeta = ndeta + 1
c
      nalpha = nalpha - 1
      ndeta = 3
      do 170 i1=1,3
      jdbl(nalpha) = iopn(i1)
      do 160 jj=1,nalpha
  160 ifigb(jj,ndeta) = jdbl(jj)
  170 ndeta = ndeta - 1
c
c  determinant with ms = (-3/2)
c
      ms(4) = msbas + (1 - 3*isgn)/2
      do 180 jj = 1,kdbl
  180 ifiga(jj,4) = jdbl(jj)
      do 190 i1=1,3
  190 jdbl(kdbl+i1) = iopn(i1)
      nalpha = nalpha + 2
      do 200 jj = 1,nalpha
  200 ifigb(jj,4) = jdbl(jj)
c
      go to 500
c
c  case 3:  five open-shell electrons
c
c          number of double-group functions = 16
c          number of determinants           = 16
c
c  arrangement of the singly occupied orbitals
c
c          singly occupied orbitals
c     det       alpha list        beta list     (sign) spin function
c
c      1       1  2  3            4  5           (+) a a a b b
c      2       1  2  4            3  5           (-) a a b a b
c      3       1  3  4            2  5           (+) a b a a b
c      4       2  3  4            1  5           (-) b a a a b
c      5       1  2  5            3  4           (+) a a b b a
c      6       1  3  5            2  4           (-) a b a b a
c      7       2  3  5            1  4           (+) b a a b a
c      8       1  4  5            2  3           (+) a b b a a
c      9       2  4  5            1  3           (-) b a b a a
c     10       3  4  5            1  2           (+) b b a a a
c     11       5                  1  2  3  4     (+) b b b b a
c     12       4                  1  2  3  5     (-) b b b a b
c     13       3                  1  2  4  5     (+) b b a b b
c     14       2                  1  3  4  5     (-) b a b b b
c     15       1                  2  3  4  5     (+) a b b b b
c     16       1  2  3  4  5                     (+) a a a a a
c
c  (alpha), then (beta), lists for determinants with ms = (1/2)
c
  210 ndeta = 1
      do 230 i3=3,5
      jdbl(nalpha) = iopn(i3)
      do 230 i2=2,i3-1
      jdbl(nalpha-1) = iopn(i2)
      do 230 i1=1,i2-1
      jdbl(nalpha-2) = iopn(i1)
      ms(ndeta) = (msbas + (1 + isgn)/2)
      do 220 jj=1,nalpha
  220 ifiga(jj,ndeta) = jdbl(jj)
  230 ndeta = ndeta + 1
c
      nalpha = nalpha - 1
      ndeta = 10
      do 250 i2=2,5
      jdbl(nalpha) = iopn(i2)
      do 250 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      do 240 jj=1,nalpha
  240 ifigb(jj,ndeta) = jdbl(jj)
  250 ndeta = ndeta - 1
c
c  (beta), then (alpha), lists for determinants with ms = (-3/2)
c
      nalpha = nalpha + 2
      ndeta = 11
      do 270 i4=4,5
      jdbl(nalpha) = iopn(i4)
      do 270 i3=3,i4-1
      jdbl(nalpha-1) = iopn(i3)
      do 270 i2=2,i3-1
      jdbl(nalpha-2) = iopn(i2)
      do 270 i1=1,i2-1
      jdbl(nalpha-3) = iopn(i1)
      ms(ndeta) = (msbas + (1 - 3*isgn)/2)
      do 260 jj=1,nalpha
  260 ifigb(jj,ndeta) = jdbl(jj)
  270 ndeta = ndeta + 1
c
      nalpha = nalpha - 3
      ndeta = 15
      do 290 i1=1,5
      jdbl(nalpha) = iopn(i1)
      do 280 jj=1,nalpha
  280 ifiga(jj,ndeta) = jdbl(jj)
  290 ndeta = ndeta - 1
c
c  determinant with ms = (5/2)
c
      ms(16) = msbas + (1 + 5*isgn)/2
      do 300 jj=1,kdbl
  300 ifigb(jj,16) = jdbl(jj)
      do 310 i1=1,5
  310 jdbl(kdbl+i1) = iopn(i1)
      nalpha = nalpha + 4
      do 320 jj=1,nalpha
  320 ifiga(jj,16) = jdbl(jj)
c
      go to 500
c
c  case 4:  seven open-shell electrons
c
c          number of double-group functions = 64
c          number of determinants           = 64
c
c  (alpha), then (beta), lists for determinants with ms = (1/2)
c
  330 ndeta = 1
      do 350 i4=4,7
      jdbl(nalpha) = iopn(i4)
      do 350 i3=3,i4-1
      jdbl(nalpha-1) = iopn(i3)
      do 350 i2=2,i3-1
      jdbl(nalpha-2) = iopn(i2)
      do 350 i1=1,i2-1
      jdbl(nalpha-3) = iopn(i1)
      ms(ndeta) = (msbas + (1 + isgn)/2)
      do 340 jj = 1,nalpha
  340 ifiga(jj,ndeta) = jdbl(jj)
  350 ndeta = ndeta + 1
c
      nalpha = nalpha - 1
      ndeta = 35
      do 370 i3=3,7
      jdbl(nalpha) = iopn(i3)
      do 370 i2=2,i3-1
      jdbl(nalpha-1) = iopn(i2)
      do 370 i1=1,i2-1
      jdbl(nalpha-2) = iopn(i1)
      do 360 jj = 1,nalpha
  360 ifigb(jj,ndeta) = jdbl(jj)
  370 ndeta = ndeta - 1
c
c  (beta), then (alpha), lists for determinants with ms = (-3/2)
c
      nalpha = nalpha + 2
      ndeta = 36
      do 390 i5=5,7
      jdbl(nalpha) = iopn(i5)
      do 390 i4=4,i5-1
      jdbl(nalpha-1) = iopn(i4)
      do 390 i3=3,i4-1
      jdbl(nalpha-2) = iopn(i3)
      do 390 i2=2,i3-1
      jdbl(nalpha-3) = iopn(i2)
      do 390 i1=1,i2-1
      jdbl(nalpha-4) = iopn(i1)
      ms(ndeta) = (msbas + (1 - 3*isgn)/2)
      do 380 jj=1,nalpha
  380 ifigb(jj,ndeta) = jdbl(jj)
  390 ndeta = ndeta + 1
c
      nalpha = nalpha - 3
      ndeta = 56
      do 410 i2=2,7
      jdbl(nalpha) = iopn(i2)
      do 410 i1=1,i2-1
      jdbl(nalpha-1) = iopn(i1)
      do 400 jj=1,nalpha
  400 ifiga(jj,ndeta) = jdbl(jj)
  410 ndeta = ndeta - 1
c
c  (alpha), then (beta), lists for determinants with ms = (5/2)
c
      nalpha = nalpha + 4
      ndeta = 57
      do 430 i6=6,7
      jdbl(nalpha) = iopn(i6)
      do 430 i5=5,i6-1
      jdbl(nalpha-1) = iopn(i5)
      do 430 i4=4,i5-1
      jdbl(nalpha-2) = iopn(i4)
      do 430 i3=3,i4-1
      jdbl(nalpha-3) = iopn(i3)
      do 430 i2=2,i3-1
      jdbl(nalpha-4) = iopn(i2)
      do 430 i1=1,i2-1
      jdbl(nalpha-5) = iopn(i1)
      ms(ndeta) = (msbas + (1 + 5*isgn)/2)
      do 420 jj=1,nalpha
  420 ifiga(jj,ndeta) = jdbl(jj)
  430 ndeta = ndeta + 1
c
      nalpha = nalpha - 5
      ndeta = 63
      do 450 i1=1,7
      jdbl(nalpha) = iopn(i1)
      do 440 jj=1,nalpha
  440 ifigb(jj,ndeta) = jdbl(jj)
  450 ndeta = ndeta - 1
c
c  determinant with ms = (-7/2)
c
      ms(64) = msbas + (1 - 7*isgn)/2
      do 460 jj=1,kdbl
  460 ifiga (jj,64) = jdbl(jj)
      do 470 i1=1,7
  470 jdbl(kdbl+i1) = iopn(i1)
      nalpha = nalpha + 6
      do 480 jj=1,nalpha
  480 ifigb(jj,64) = jdbl(jj)
c
  500 return
c
      end

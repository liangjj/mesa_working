*deck dmout
      subroutine dmout (m, n, lda, a, ifmt, idigit)
c***begin prologue  dmout
c***subsidiary
c***purpose  subsidiary to dbocls and dfc
c***library   slatec
c***type      double precision (smout-s, dmout-d)
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c     double precision matrix output routine.
c
c  input..
c
c  m,n,lda,a(*,*) print the double precision array a(i,j),i = 1,...,m,
c                 j=1,...,n, on output unit lout=6. lda is the declared
c                 first dimension of a(*,*) as specified in the calling
c                 program. the heading in the fortran format statement
c                 ifmt(*), described below, is printed as a first step.
c                 the components a(i,j) are indexed, on output, in a
c                 pleasant format.
c  ifmt(*)        a fortran format statement. this is printed on
c                 output unit lout=6 with the variable format fortran
c                 statement
c                       write(lout,ifmt).
c  idigit         print at least abs(idigit) decimal digits per number.
c                 the subprogram will choose that integer 4,6,14,20 or
c                 28 which will print at least abs(idigit) number of
c                 places.  if idigit.lt.0, 72 printing columns are
c                 utilized to write each line of output of the array
c                 a(*,*). (this can be used on most time-sharing
c                 terminals).  if idigit.ge.0, 133 printing columns are
c                 utilized. (this can be used on most line printers).
c
c  example..
c
c  print an array called (simplex tableau   ) of size 10 by 20 showing
c  6 decimal digits per number. the user is running on a time-sharing
c  system with a 72 column output device.
c
c     double precision tableu(20,20)
c     m = 10
c     n = 20
c     ldtabl = 20
c     idigit = -6
c     call dmout(m,n,ldtabl,tableu,21h(16h1simplex tableau),idigit)
c
c***see also  dbocls, dfc
c***routines called  i1mach
c***revision history  (yymmdd)
c   821220  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891107  added comma after 1p edit descriptor in format
c           statements.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910403  updated author section.  (wrb)
c***end prologue  dmout
      double precision a(lda,*)
      character ifmt*(*),icol*3
      save icol
      data icol /'col'/
c***first executable statement  dmout
      lout=i1mach(2)
      write(lout,ifmt)
      if(m.le.0.or.n.le.0.or.lda.le.0) return
      ndigit = idigit
      if(idigit.eq.0) ndigit = 4
      if(idigit.ge.0) go to 80
c
      ndigit = -idigit
      if(ndigit.gt.4) go to 9
c
      do 5 k1=1,n,5
      k2 = min(n,k1+4)
      write(lout,1010) (icol,i,i = k1, k2)
      do 5 i = 1, m
      write(lout,1009) i,(a(i,j),j = k1, k2)
   5  continue
      return
c
   9  continue
      if(ndigit.gt.6) go to 20
c
      do 10 k1=1,n,4
      k2 = min(n,k1+3)
      write(lout,1000) (icol,i,i = k1, k2)
      do 10 i = 1, m
      write(lout,1004) i,(a(i,j),j = k1, k2)
   10 continue
      return
c
   20 continue
      if(ndigit.gt.14) go to 40
c
      do 30 k1=1,n,2
      k2 = min(n,k1+1)
      write(lout,1001) (icol,i,i = k1, k2)
      do 30 i = 1, m
      write(lout,1005) i,(a(i,j),j = k1, k2)
   30 continue
      return
c
   40 continue
      if(ndigit.gt.20) go to 60
c
      do 50 k1=1,n,2
      k2=min(n,k1+1)
      write(lout,1002) (icol,i,i = k1, k2)
      do 50 i = 1, m
      write(lout,1006) i,(a(i,j),j = k1, k2)
   50 continue
      return
c
   60 continue
      do 70 k1=1,n
      k2 = k1
      write(lout,1003) (icol,i,i = k1, k2)
      do 70 i = 1, m
      write(lout,1007) i,(a(i,j),j = k1, k2)
   70 continue
      return
c
   80 continue
      if(ndigit.gt.4) go to 86
c
      do 85 k1=1,n,10
      k2 = min(n,k1+9)
      write(lout,1000) (icol,i,i = k1, k2)
      do 85 i = 1, m
      write(lout,1009) i,(a(i,j),j = k1, k2)
   85 continue
c
86    if (ndigit.gt.6) go to 100
c
      do 90 k1=1,n,8
      k2 = min(n,k1+7)
      write(lout,1000) (icol,i,i = k1, k2)
      do 90 i = 1, m
      write(lout,1004) i,(a(i,j),j = k1, k2)
   90 continue
      return
c
  100 continue
      if(ndigit.gt.14) go to 120
c
      do 110 k1=1,n,5
      k2 = min(n,k1+4)
      write(lout,1001) (icol,i,i = k1, k2)
      do 110 i = 1, m
      write(lout,1005) i,(a(i,j),j = k1, k2)
  110 continue
      return
c
  120 continue
      if(ndigit.gt.20) go to 140
c
      do 130 k1=1,n,4
      k2 = min(n,k1+3)
      write(lout,1002) (icol,i,i = k1, k2)
      do 130 i = 1, m
      write(lout,1006) i,(a(i,j),j = k1, k2)
  130 continue
      return
c
  140 continue
      do 150 k1=1,n,3
      k2 = min(n,k1+2)
      write(lout,1003) (icol,i,i = k1, k2)
      do 150 i = 1, m
      write(lout,1007) i,(a(i,j),j = k1, k2)
  150 continue
      return
 1000 format(10x,8(5x,a,i4,2x))
 1001 format(10x,5(9x,a,i4,6x))
 1002 format(10x,4(12x,a,i4,9x))
 1003 format(10x,3(16x,a,i4,13x))
 1004 format(1x,3hrow,i4,2x,1p,8d14.5)
 1005 format(1x,3hrow,i4,2x,1p,5d22.13)
 1006 format(1x,3hrow,i4,2x,1p,4d28.19)
 1007 format(1x,3hrow,i4,2x,1p,3d36.27)
 1009 format(1x,3hrow,i4,2x,1p,10d12.3)
 1010 format(10x,10(4x,a,i4,1x))
      end

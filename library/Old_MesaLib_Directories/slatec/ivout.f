*deck ivout
      subroutine ivout (n, ix, ifmt, idigit)
c***begin prologue  ivout
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      integer (ivout-i)
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c     integer vector output routine.
c
c  input..
c
c  n,ix(*) print the integer array ix(i),i=1,...,n, on output
c          unit lout. the heading in the fortran format
c          statement ifmt(*), described below, is printed as a first
c          step. the components ix(i) are indexed, on output,
c          in a pleasant format.
c  ifmt(*) a fortran format statement. this is printed on output
c          unit lout with the variable format fortran statement
c                write(lout,ifmt)
c  idigit  print up to abs(idigit) decimal digits per number.
c          the subprogram will choose that integer 4,6,10 or 14
c          which will print at least abs(idigit) number of
c          places.  if idigit.lt.0, 72 printing columns are utilized
c          to write each line of output of the array ix(*). (this
c          can be used on most time-sharing terminals). if
c          idigit.ge.0, 133 printing columns are utilized. (this can
c          be used on most line printers).
c
c  example..
c
c  print an array called (costs of purchases) of length 100 showing
c  6 decimal digits per number. the user is running on a time-sharing
c  system with a 72 column output device.
c
c     dimension icosts(100)
c     n = 100
c     idigit = -6
c     call ivout(n,icosts,'(''1costs of purchases'')',idigit)
c
c***see also  splp
c***routines called  i1mach
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   910403  updated author section.  (wrb)
c***end prologue  ivout
      dimension ix(*)
      character ifmt*(*)
c
c     get the unit number where output will be written.
c***first executable statement  ivout
      j=2
      lout=i1mach(j)
      write(lout,ifmt)
      if(n.le.0) return
      ndigit = idigit
      if(idigit.eq.0) ndigit = 4
      if(idigit.ge.0) go to 80
c
      ndigit = -idigit
      if(ndigit.gt.4) go to 20
c
      do 10 k1=1,n,10
      k2 = min(n,k1+9)
      write(lout,1000) k1,k2,(ix(i),i=k1,k2)
   10 continue
      return
c
   20 continue
      if(ndigit.gt.6) go to 40
c
      do 30 k1=1,n,7
      k2 = min(n,k1+6)
      write(lout,1001) k1,k2,(ix(i),i=k1,k2)
   30 continue
      return
c
   40 continue
      if(ndigit.gt.10) go to 60
c
      do 50 k1=1,n,5
      k2=min(n,k1+4)
      write(lout,1002) k1,k2,(ix(i),i=k1,k2)
   50 continue
      return
c
   60 continue
      do 70 k1=1,n,3
      k2 = min(n,k1+2)
      write(lout,1003) k1,k2,(ix(i),i=k1,k2)
   70 continue
      return
c
   80 continue
      if(ndigit.gt.4) go to 100
c
      do 90 k1=1,n,20
      k2 = min(n,k1+19)
      write(lout,1000) k1,k2,(ix(i),i=k1,k2)
   90 continue
      return
c
  100 continue
      if(ndigit.gt.6) go to 120
c
      do 110 k1=1,n,15
      k2 = min(n,k1+14)
      write(lout,1001) k1,k2,(ix(i),i=k1,k2)
  110 continue
      return
c
  120 continue
      if(ndigit.gt.10) go to 140
c
      do 130 k1=1,n,10
      k2 = min(n,k1+9)
      write(lout,1002) k1,k2,(ix(i),i=k1,k2)
  130 continue
      return
c
  140 continue
      do 150 k1=1,n,7
      k2 = min(n,k1+6)
      write(lout,1003) k1,k2,(ix(i),i=k1,k2)
  150 continue
      return
 1000 format(1x,i4,' - ',i4,20(1x,i5))
 1001 format(1x,i4,' - ',i4,15(1x,i7))
 1002 format(1x,i4,' - ',i4,10(1x,i11))
 1003 format(1x,i4,' - ',i4,7(1x,i15))
      end

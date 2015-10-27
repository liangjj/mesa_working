*deck dbnslv
      subroutine dbnslv (w, nroww, nrow, nbandl, nbandu, b)
c***begin prologue  dbnslv
c***subsidiary
c***purpose  subsidiary to dbint4 and dbintk
c***library   slatec
c***type      double precision (bnslv-s, dbnslv-d)
c***author  (unknown)
c***description
c
c  dbnslv is the banslv routine from
c        * a practical guide to splines *  by c. de boor
c
c  dbnslv is a double precision routine
c
c  companion routine to  dbnfac . it returns the solution  x  of the
c  linear system  a*x = b  in place of  b , given the lu-factorization
c  for  a  in the work array  w from dbnfac.
c
c *****  i n p u t  ****** w,b are double precision
c  w, nroww,nrow,nbandl,nbandu.....describe the lu-factorization of a
c        banded matrix  a  of order  nrow  as constructed in  dbnfac .
c        for details, see  dbnfac .
c  b.....right side of the system to be solved .
c
c *****  o u t p u t  ****** b is double precision
c  b.....contains the solution  x , of order  nrow .
c
c *****  m e t h o d  ******
c     (with  a = l*u, as stored in  w,) the unit lower triangular system
c  l(u*x) = b  is solved for  y = u*x, and  y  stored in  b . then the
c  upper triangular system  u*x = y  is solved for  x  . the calcul-
c  ations are so arranged that the innermost loops stay within columns.
c
c***see also  dbint4, dbintk
c***routines called  (none)
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dbnslv
c
      integer nbandl, nbandu, nrow, nroww, i, j, jmax, middle, nrowm1
      double precision w(nroww,*), b(*)
c***first executable statement  dbnslv
      middle = nbandu + 1
      if (nrow.eq.1) go to 80
      nrowm1 = nrow - 1
      if (nbandl.eq.0) go to 30
c                                 forward pass
c            for i=1,2,...,nrow-1, subtract  right side(i)*(i-th column
c            of  l )  from right side  (below i-th row) .
      do 20 i=1,nrowm1
        jmax = min(nbandl,nrow-i)
        do 10 j=1,jmax
          b(i+j) = b(i+j) - b(i)*w(middle+j,i)
   10   continue
   20 continue
c                                 backward pass
c            for i=nrow,nrow-1,...,1, divide right side(i) by i-th diag-
c            onal entry of  u, then subtract  right side(i)*(i-th column
c            of  u)  from right side  (above i-th row).
   30 if (nbandu.gt.0) go to 50
c                                a  is lower triangular .
      do 40 i=1,nrow
        b(i) = b(i)/w(1,i)
   40 continue
      return
   50 i = nrow
   60 b(i) = b(i)/w(middle,i)
      jmax = min(nbandu,i-1)
      do 70 j=1,jmax
        b(i-j) = b(i-j) - b(i)*w(middle-j,i)
   70 continue
      i = i - 1
      if (i.gt.1) go to 60
   80 b(1) = b(1)/w(middle,1)
      return
      end

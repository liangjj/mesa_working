*deck @(#)initds.f	5.1  4/18/95
      function initds (dos, nos, eta)
      integer initds
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c initialize the double precision orthogonal series dos so that initds
c is the number of terms needed to insure the error is no larger than
c eta.  ordinarily eta will be chosen to be one-tenth machine precision.
c
c             input arguments --
c dos    dble prec array of nos coefficients in an orthogonal series.
c nos    number of coefficients in dos.
c eta    requested accuracy of series.
c
      real*8 dos(nos)
c
      if (nos.lt.1) call lnkerr (
     1  'initds: number of coefficients lt 1')
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call lnkerr ('initds: eta may be too small')
      initds = i
c
      return
      end

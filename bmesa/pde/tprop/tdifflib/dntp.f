*deck dntp
      subroutine dntp (x, y, xout, yout, ypout, neqn, kold, phi, ivc,
     +   iv, kgi, gi, alpha, og, ow, ox, oy)
c***begin prologue  dintp
c***purpose  approximate the solution at xout by evaluating the
c            polynomial computed in dsteps at xout.  must be used in
c            conjunction with dsteps.
c***library   slatec (depac)
c***category  i1a1b
c***type      double precision (sintrp-s, dintp-d)
c***keywords  adams method, depac, initial value problems, ode,
c             ordinary differential equations, predictor-corrector,
c             smooth interpolant
c***author  watts, h. a., (snla)
c***description
c
c   the methods in subroutine  dsteps  approximate the solution near  x
c   by a polynomial.  subroutine  dintp  approximates the solution at
c   xout  by evaluating the polynomial there.  information defining this
c   polynomial is passed from  dsteps  so  dintp  cannot be used alone.
c
c   subroutine dsteps is completely explained and documented in the text
c   "computer solution of ordinary differential equations, the initial
c   value problem"  by l. f. shampine and m. k. gordon.
c
c   input to dintp --
c
c   the user provides storage in the calling program for the arrays in
c   the call list
c      dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),oy(neqn)
c                and alpha(12),og(13),ow(12),gi(11),iv(10)
c   and defines
c      xout -- point at which solution is desired.
c   the remaining parameters are defined in  dsteps  and passed to
c   dintp  from that subroutine
c
c   output from  dintp --
c
c      yout(*) -- solution at  xout
c      ypout(*) -- derivative of solution at  xout
c   the remaining parameters are returned unaltered from their input
c   values.  integration with  dsteps  may be continued.
c
c***references  h. a. watts, a smoother interpolant for de/step, intrp
c                 ii, report sand84-0293, sandia laboratories, 1984.
c***routines called  (none)
c***revision history  (yymmdd)
c   840201  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dintp
c
      implicit integer(a-z)
      real*8 alp, alpha, c, g, gdi, gdif, gi, gamma, h, hi,
     1       hmu, og, ow, ox, oy, phi, rmu, sigma, temp1, temp2, temp3,
     2       w, x, xi, xim1, xiq, xout, y, yout, ypout
c
      dimension y(*),yout(*),ypout(*),phi(neqn,16),oy(*)
      dimension g(13),c(13),w(13),og(13),ow(12),alpha(12),gi(11),iv(10)
c
c***first executable statement  dintp
      kp1 = kold + 1
      kp2 = kold + 2
c
      hi = xout - ox
      h = x - ox
      xi = hi/h
      xim1 = xi - 1.d0
c
c   initialize w(*) for computing g(*)
c
      xiq = xi
      do 10 iq = 1,kp1
        xiq = xi*xiq
        temp1 = iq*(iq+1)
 10     w(iq) = xiq/temp1
c
c   compute the double integral term gdi
c
      if (kold .le. kgi) go to 50
      if (ivc .gt. 0) go to 20
      gdi = 1.0d0/temp1
      m = 2
      go to 30
 20   iw = iv(ivc)
      gdi = ow(iw)
      m = kold - iw + 3
 30   if (m .gt. kold) go to 60
      do 40 i = m,kold
 40     gdi = ow(kp2-i) - alpha(i)*gdi
      go to 60
 50   gdi = gi(kold)
c
c   compute g(*) and c(*)
c
 60   g(1) = xi
      g(2) = 0.5d0*xi*xi
      c(1) = 1.0d0
      c(2) = xi
      if (kold .lt. 2) go to 90
      do 80 i = 2,kold
        alp = alpha(i)
        gamma = 1.0d0 + xim1*alp
        l = kp2 - i
        do 70 jq = 1,l
 70       w(jq) = gamma*w(jq) - alp*w(jq+1)
        g(i+1) = w(1)
 80     c(i+1) = gamma*c(i)
c
c   define interpolation parameters
c
 90   sigma = (w(2) - xim1*w(1))/gdi
      rmu = xim1*c(kp1)/gdi
      hmu = rmu/h
c
c   interpolate for the solution -- yout
c   and for the derivative of the solution -- ypout
c
      do 100 l = 1,neqn
        yout(l) = 0.0d0
 100    ypout(l) = 0.0d0
      do 120 j = 1,kold
        i = kp2 - j
        gdif = og(i) - og(i-1)
        temp2 = (g(i) - g(i-1)) - sigma*gdif
        temp3 = (c(i) - c(i-1)) + rmu*gdif
        do 110 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
 110      ypout(l) = ypout(l) + temp3*phi(l,i)
 120    continue
      do 130 l = 1,neqn
        yout(l) = ((1.0d0 - sigma)*oy(l) + sigma*y(l)) +
     1             h*(yout(l) + (g(1) - sigma*og(1))*phi(l,1))
 130    ypout(l) = hmu*(oy(l) - y(l)) +
     1                (ypout(l) + (c(1) + rmu*og(1))*phi(l,1))
c
      return
      end

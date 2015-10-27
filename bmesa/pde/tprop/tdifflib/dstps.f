*deck dstps
      subroutine dstps ( neqn, y, x, h, eps, wt, start, hold, k,
     +   kold, crash, phi, p, yp, psi, alpha, beta, sig, v, w, g,
     +   phase1, ns, nornd, ksteps, twou, fouru, xold, kprev, ivc, iv,
     +   kgi, gi, rpar, ipar,
     +   h01, h02, h03, eig1, eig2, eig3, v0, vt, tpert, efield, omega, 
     +   width, shift, hbar, dim, nmax, n)
c***begin prologue  dsteps
c***purpose  integrate a system of first order ordinary differential
c            equations one step.
c***library   slatec (depac)
c***category  i1a1b
c***type      double precision (steps-s, dsteps-d)
c***keywords  adams method, depac, initial value problems, ode,
c             ordinary differential equations, predictor-corrector
c***author  shampine, l. f., (snla)
c           gordon, m. k., (snla)
c             modified by h.a. watts
c***description
c
c   written by l. f. shampine and m. k. gordon
c
c   abstract
c
c   subroutine  dsteps  is normally used indirectly through subroutine
c   ddeabm .  because  ddeabm  suffices for most problems and is much
c   easier to use, using it should be considered before using  dsteps
c   alone.
c
c   subroutine dsteps integrates a system of  neqn  first order ordinary
c   differential equations one step, normally from x to x+h, using a
c   modified divided difference form of the adams pece formulas.  local
c   extrapolation is used to improve absolute stability and accuracy.
c   the code adjusts its order and step size to control the local error
c   per unit step in a generalized sense.  special devices are included
c   to control roundoff error and to detect when the user is requesting
c   too much accuracy.
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations, the initial
c   value problem  by l. f. shampine and m. k. gordon.
c   further details on use of this code are available in "solving
c   ordinary differential equations with ode, step, and intrp",
c   by l. f. shampine and m. k. gordon, sla-73-1060.
c
c
c   the parameters represent --
c      df -- subroutine to evaluate derivatives
c      neqn -- number of equations to be integrated
c      y(*) -- solution vector at x
c      x -- independent variable
c      h -- appropriate step size for next step.  normally determined by
c           code
c      eps -- local error tolerance
c      wt(*) -- vector of weights for error criterion
c      start -- logical variable set .true. for first step,  .false.
c           otherwise
c      hold -- step size used for last successful step
c      k -- appropriate order for next step (determined by code)
c      kold -- order used for last successful step
c      crash -- logical variable set .true. when no step can be taken,
c           .false. otherwise.
c      yp(*) -- derivative of solution vector at  x  after successful
c           step
c      ksteps -- counter on attempted steps
c      twou -- 2.*u where u is machine unit roundoff quantity
c      fouru -- 4.*u where u is machine unit roundoff quantity
c      rpar,ipar -- parameter arrays which you may choose to use
c            for communication between your program and subroutine f.
c            they are not altered or used by dsteps.
c   the variables x,xold,kold,kgi and ivc and the arrays y,phi,alpha,g,
c   w,p,iv and gi are required for the interpolation subroutine sintrp.
c   the remaining variables and arrays are included in the call list
c   only to eliminate local retention of variables between calls.
c
c   input to dsteps
c
c      first call --
c
c   the user must provide storage in his calling program for all arrays
c   in the call list, namely
c
c     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12),
c    1  alpha(12),beta(12),sig(13),v(12),w(12),g(13),gi(11),iv(10),
c    2  rpar(*),ipar(*)
c
c    **note**
c
c   the user must also declare  start ,  crash ,  phase1  and  nornd
c   logical variables and  df  an external subroutine, supply the
c   subroutine  df(x,y,yp)  to evaluate
c      dy(i)/dx = yp(i) = df(x,y(1),y(2),...,y(neqn))
c   and initialize only the following parameters.
c      neqn -- number of equations to be integrated
c      y(*) -- vector of initial values of dependent variables
c      x -- initial value of the independent variable
c      h -- nominal step size indicating direction of integration
c           and maximum size of step.  must be variable
c      eps -- local error tolerance per step.  must be variable
c      wt(*) -- vector of non-zero weights for error criterion
c      start -- .true.
c      yp(*) -- vector of initial derivative values
c      ksteps -- set ksteps to zero
c      twou -- 2.*u where u is machine unit roundoff quantity
c      fouru -- 4.*u where u is machine unit roundoff quantity
c   define u to be the machine unit roundoff quantity by calling
c   the function routine  r1mach,  u = r1mach(4), or by
c   computing u so that u is the smallest positive number such
c   that 1.0+u .gt. 1.0.
c
c   dsteps  requires that the l2 norm of the vector with components
c   local error(l)/wt(l)  be less than  eps  for a successful step.  the
c   array  wt  allows the user to specify an error test appropriate
c   for his problem.  for example,
c      wt(l) = 1.0  specifies absolute error,
c            = abs(y(l))  error relative to the most recent value of the
c                 l-th component of the solution,
c            = abs(yp(l))  error relative to the most recent value of
c                 the l-th component of the derivative,
c            = max(wt(l),abs(y(l)))  error relative to the largest
c                 magnitude of l-th component obtained so far,
c            = abs(y(l))*relerr/eps + abserr/eps  specifies a mixed
c                 relative-absolute test where  relerr  is relative
c                 error,  abserr  is absolute error and  eps =
c                 max(relerr,abserr) .
c
c      subsequent calls --
c
c   subroutine  dsteps  is designed so that all information needed to
c   continue the integration, including the step size  h  and the order
c   k , is returned with each step.  with the exception of the step
c   size, the error tolerance, and the weights, none of the parameters
c   should be altered.  the array  wt  must be updated after each step
c   to maintain relative error tests like those above.  normally the
c   integration is continued just beyond the desired endpoint and the
c   solution interpolated there with subroutine  sintrp .  if it is
c   impossible to integrate beyond the endpoint, the step size may be
c   reduced to hit the endpoint since the code will not take a step
c   larger than the  h  input.  changing the direction of integration,
c   i.e., the sign of  h , requires the user set  start = .true. before
c   calling  dsteps  again.  this is the only situation in which  start
c   should be altered.
c
c   output from dsteps
c
c      successful step --
c
c   the subroutine returns after each successful step with  start  and
c   crash  set .false. .  x  represents the independent variable
c   advanced one step of length  hold  from its value on input and  y
c   the solution vector at the new value of  x .  all other parameters
c   represent information corresponding to the new  x  needed to
c   continue the integration.
c
c      unsuccessful step --
c
c   when the error tolerance is too small for the machine precision,
c   the subroutine returns without taking a step and  crash = .true. .
c   an appropriate step size and error tolerance for continuing are
c   estimated and all other information is restored as upon input
c   before returning.  to continue with the larger tolerance, the user
c   just calls the code again.  a restart is neither required nor
c   desirable.
c
c***references  l. f. shampine and m. k. gordon, solving ordinary
c                 differential equations with ode, step, and intrp,
c                 report sla-73-1060, sandia laboratories, 1973.
c***routines called  r1mach, dhstrt
c***revision history  (yymmdd)
c   740101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dsteps
c
      implicit integer(a-z)
      real*8 absh, alpha, beta, big, r1mach,
     1      eps, erk, erkm1, erkm2, erkp1, err,
     2      fouru, g, gi, gstr, h, hnew, hold, p, p5eps, phi, psi, r,
     3      reali, realns, rho, round, rpar, sig, tau, temp1,
     4      temp2, temp3, temp4, temp5, temp6, two, twou, u, v, w, wt,
     5      x, xold, y, yp
      logical start,crash,phase1,nornd
      dimension y(*),wt(*),phi(neqn,16),p(*),yp(*),psi(12),
     1  alpha(12),beta(12),sig(13),v(12),w(12),g(13),gi(11),iv(10),
     2  rpar(*),ipar(*)
      dimension two(13),gstr(13)
      real*8 h01, h02, h03, eig1, eig2, eig3
      real*8 v0, vt, omega, efield, width, shift, hbar
      character*(*) tpert
      dimension nmax(3)
      dimension h01(nmax(1),*), h02(nmax(2),*), h03(nmax(3),*)
      dimension eig1(nmax(1)), eig2(nmax(2)), eig3(nmax(3))
      dimension v0(n), vt(n)
      common/io/inp, iout
c      external df
      save two, gstr
c
      data two(1),two(2),two(3),two(4),two(5),two(6),two(7),two(8),
     1     two(9),two(10),two(11),two(12),two(13)
     2     /2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,
     3      512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
      data gstr(1),gstr(2),gstr(3),gstr(4),gstr(5),gstr(6),gstr(7),
     1     gstr(8),gstr(9),gstr(10),gstr(11),gstr(12),gstr(13)
     2     /0.5d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0,0.0114d0,
     3      0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,0.00468d0/
c
c       ***     begin block 0     ***
c   check if step size or error tolerance is too small for machine
c   precision.  if first step, initialize phi array and estimate a
c   starting step size.
c                   ***
c
c   if step size is too small, determine an acceptable one
c
c***first executable statement  dsteps
      crash = .true.
      if(abs(h) .ge. fouru*abs(x)) go to 5
      h = sign(fouru*abs(x),h)
      return
 5    p5eps = 0.5d0*eps
c
c   if error tolerance is too small, increase it to an acceptable value
c
      round = 0.0d0
      do 10 l = 1,neqn
 10     round = round + (y(l)/wt(l))**2
      round = twou*sqrt(round)
      if(p5eps .ge. round) go to 15
      eps = 2.0d0*round*(1.0d0 + fouru)
      return
 15   crash = .false.
      g(1) = 1.0d0
      g(2) = 0.5d0
      sig(1) = 1.0d0
      if(.not.start) go to 99
c
c   initialize.  compute appropriate step size for first step
c
c     sum = 0.0
      do 20 l = 1,neqn
        phi(l,1) = yp(l)
   20   phi(l,2) = 0.0d0
c20     sum = sum + (yp(l)/wt(l))**2
c     sum = sqrt(sum)
c     absh = abs(h)
c     if(eps .lt. 16.0*sum*h*h) absh = 0.25*sqrt(eps/sum)
c     h = sign(max(absh,fouru*abs(x)),h)
c
      u = r1mach(4)
      big = sqrt(r1mach(2))
      call dstrt(neqn,x,x+h,y,yp,wt,1,u,big,
     1           phi(1,3),phi(1,4),phi(1,5),phi(1,6),h,rpar,ipar,
     2           h01,h02,h03,eig1,eig2,eig3,v0,vt,tpert,
     3           efield,omega,width,shift,hbar,dim,nmax,n)
      hold = 0.0d0
      k = 1
      kold = 0
      kprev = 0
      start = .false.
      phase1 = .true.
      nornd = .true.
      if(p5eps .gt. 100.0d0*round) go to 99
      nornd = .false.
      do 25 l = 1,neqn
 25     phi(l,15) = 0.0d0
 99   ifail = 0
c       ***     end block 0     ***
c
c       ***     begin block 1     ***
c   compute coefficients of formulas for this step.  avoid computing
c   those quantities not changed when step size is not changed.
c                   ***
c
 100  kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
c
c   ns is the number of dsteps taken with size h, including the current
c   one.  when k.lt.ns, no coefficients change
c
      if(h .ne. hold) ns = 0
      if (ns.le.kold) ns = ns+1
      nsp1 = ns+1
      if (k .lt. ns) go to 199
c
c   compute those components of alpha(*),beta(*),psi(*),sig(*) which
c   are changed
c
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      if(k .lt. nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
 105    sig(i+1) = reali*alpha(i)*sig(i)
 110  psi(k) = temp1
c
c   compute coefficients g(*)
c
c   initialize v(*) and set w(*).
c
      if(ns .gt. 1) go to 120
      do 115 iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
 115    w(iq) = v(iq)
      ivc = 0
      kgi = 0
      if (k .eq. 1) go to 140
      kgi = 1
      gi(1) = w(2)
      go to 140
c
c   if order was raised, update diagonal part of v(*)
c
 120  if(k .le. kprev) go to 130
      if (ivc .eq. 0) go to 122
      jv = kp1 - iv(ivc)
      ivc = ivc - 1
      go to 123
 122  jv = 1
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      w(k) = v(k)
      if (k .ne. 2) go to 123
      kgi = 1
      gi(1) = w(2)
 123  nsm2 = ns-2
      if(nsm2 .lt. jv) go to 130
      do 125 j = jv,nsm2
        i = k-j
        v(i) = v(i) - alpha(j+1)*v(i+1)
 125    w(i) = v(i)
      if (i .ne. 2) go to 130
      kgi = ns - 1
      gi(kgi) = w(2)
c
c   update v(*) and set w(*)
c
 130  limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
 135    w(iq) = v(iq)
      g(nsp1) = w(1)
      if (limit1 .eq. 1) go to 137
      kgi = ns
      gi(kgi) = w(2)
 137  w(limit1+1) = v(limit1+1)
      if (k .ge. kold) go to 140
      ivc = ivc + 1
      iv(ivc) = limit1 + 2
c
c   compute the g(*) in the work vector w(*)
c
 140  nsp2 = ns + 2
      kprev = k
      if(kp1 .lt. nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
 145      w(iq) = w(iq) - temp6*w(iq+1)
 150    g(i) = w(1)
 199    continue
c       ***     end block 1     ***
c
c       ***     begin block 2     ***
c   predict a solution p(*), evaluate derivatives using predicted
c   solution, estimate local error at order k and errors at orders k,
c   k-1, k-2 as if constant step size were used.
c                   ***
c
c   increment counter on attempted dsteps
c
      ksteps = ksteps + 1
c
c   change phi to phi star
c
      if(k .lt. nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
 205      phi(l,i) = temp1*phi(l,i)
 210    continue
c
c   predict solution and differences
c
 215  do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
 220    p(l) = 0.0d0
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
 225      phi(l,i) = phi(l,i) + phi(l,ip1)
 230    continue
      if(nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
 235    phi(l,16) = (p(l) - y(l)) - tau
      go to 250
 240  do 245 l = 1,neqn
 245    p(l) = y(l) + h*p(l)
 250  xold = x
      x = x + h
      absh = abs(h)
      call df(x,p,yp,rpar,ipar,neqn,
     +        h01,h02,h03,eig1,eig2,eig3,v0,vt,tpert,efield,
     +        omega,width,shift,hbar,dim,nmax,n)
c
c   estimate errors at orders k,k-1,k-2
c
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      do 265 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
        if(km2)265,260,255
 255    erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
 260    erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
 265    erk = erk + (temp4*temp3)**2
      if(km2)280,275,270
 270  erkm2 = absh*sig(km1)*gstr(km2)*sqrt(erkm2)
 275  erkm1 = absh*sig(k)*gstr(km1)*sqrt(erkm1)
 280  temp5 = absh*sqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
c
c   test if order should be lowered
c
      if(km2)299,290,285
 285  if(max(erkm1,erkm2) .le. erk) knew = km1
      go to 299
 290  if(erkm1 .le. 0.5d0*erk) knew = km1
c
c   test if step successful
c
 299  if(err .le. eps) go to 400
c       ***     end block 2     ***
c
c       ***     begin block 3     ***
c   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
c   if third consecutive failure, set order to one.  if step fails more
c   than three times, consider an optimal step size.  double error
c   tolerance and return if estimated step size is too small for machine
c   precision.
c                   ***
c
c   restore x, phi(*,*) and psi(*)
c
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
 305      phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
 310    continue
      if(k .lt. 2) go to 320
      do 315 i = 2,k
 315    psi(i-1) = psi(i) - h
c
c   on third failure, set order to one.  thereafter, use optimal step
c   size
c
 320  ifail = ifail + 1
      temp2 = 0.5d0
      if(ifail - 3) 335,330,325
 325  if(p5eps .lt. 0.25d0*erk) temp2 = sqrt(p5eps/erk)
 330  knew = 1
 335  h = temp2*h
      k = knew
      ns = 0
      if(abs(h) .ge. fouru*abs(x)) go to 340
      crash = .true.
      h = sign(fouru*abs(x),h)
      eps = eps + eps
      return
 340  go to 100
c       ***     end block 3     ***
c
c       ***     begin block 4     ***
c   the step is successful.  correct the predicted solution, evaluate
c   the derivatives using the corrected solution and update the
c   differences.  determine best order and step size for next step.
c                   ***
 400  kold = k
      hold = h
c
c   correct and evaluate
c
      temp1 = h*g(kp1)
      if(nornd) go to 410
      do 405 l = 1,neqn
        temp3 = y(l)
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
        phi(l,15) = (y(l) - p(l)) - rho
 405    p(l) = temp3
      go to 420
 410  do 415 l = 1,neqn
        temp3 = y(l)
        y(l) = p(l) + temp1*(yp(l) - phi(l,1))
 415    p(l) = temp3
 420  call df(x,y,yp,rpar,ipar,neqn,
     +        h01,h02,h03,eig1,eig2,eig3,v0,vt,tpert,efield,omega,
     +        width,shift,hbar,dim,nmax,n)
c
c   update differences for next step
c
      do 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
 425    phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      do 435 i = 1,k
        do 430 l = 1,neqn
 430      phi(l,i) = phi(l,i) + phi(l,kp1)
 435    continue
c
c   estimate error at order k+1 unless:
c     in first phase when always raise order,
c     already decided to lower order,
c     step size not constant so estimate unreliable
c
      erkp1 = 0.0d0
      if(knew .eq. km1  .or.  k .eq. 12) phase1 = .false.
      if(phase1) go to 450
      if(knew .eq. km1) go to 455
      if(kp1 .gt. ns) go to 460
      do 440 l = 1,neqn
 440    erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*sqrt(erkp1)
c
c   using estimated error at order k+1, determine appropriate order
c   for next step
c
      if(k .gt. 1) go to 445
      if(erkp1 .ge. 0.5d0*erk) go to 460
      go to 450
 445  if(erkm1 .le. min(erk,erkp1)) go to 455
      if(erkp1 .ge. erk  .or.  k .eq. 12) go to 460
c
c   here erkp1 .lt. erk .lt. max(erkm1,erkm2) else order would have
c   been lowered in block 2.  thus order is to be raised
c
c   raise order
c
 450  k = kp1
      erk = erkp1
      go to 460
c
c   lower order
c
 455  k = km1
      erk = erkm1
c
c   with new order determine appropriate step size for next step
c
 460  hnew = h + h
      if(phase1) go to 465
      if(p5eps .ge. erk*two(k+1)) go to 465
      hnew = h
      if(p5eps .ge. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = absh*max(0.5d0,min(0.9d0,r))
      hnew = sign(max(hnew,fouru*abs(x)),h)
 465  h = hnew
      return
c       ***     end block 4     ***
      end

*deck dqnc79
      subroutine dqnc79 (fun, a, b, err, ans, ierr, k)
c***begin prologue  dqnc79
c***purpose  integrate a function using a 7-point adaptive newton-cotes
c            quadrature rule.
c***library   slatec
c***category  h2a1a1
c***type      double precision (qnc79-s, dqnc79-d)
c***keywords  adaptive quadrature, integration, newton-cotes
c***author  kahaner, d. k., (nbs)
c           jones, r. e., (snla)
c***description
c
c     abstract  *** a double precision routine ***
c       dqnc79 is a general purpose program for evaluation of
c       one dimensional integrals of user defined functions.
c       dqnc79 will pick its own points for evaluation of the
c       integrand and these will vary from problem to problem.
c       thus, dqnc79 is not designed to integrate over data sets.
c       moderately smooth integrands will be integrated efficiently
c       and reliably.  for problems with strong singularities,
c       oscillations etc., the user may wish to use more sophis-
c       ticated routines such as those in quadpack.  one measure
c       of the reliability of dqnc79 is the output parameter k,
c       giving the number of integrand evaluations that were needed.
c
c     description of arguments
c
c     --input--* fun, a, b, err are double precision *
c       fun  - name of external function to be integrated.  this name
c              must be in an external statement in your calling
c              program.  you must write a fortran function to evaluate
c              fun.  this should be of the form
c                    double precision function fun (x)
c              c
c              c     x can vary from a to b
c              c     fun(x) should be finite for all x on interval.
c              c
c                    fun = ...
c                    return
c                    end
c       a    - lower limit of integration
c       b    - upper limit of integration (may be less than a)
c       err  - is a requested error tolerance.  normally, pick a value
c              0 .lt. err .lt. 1.0d-8.
c
c     --output--
c       ans  - computed value of the integral.  hopefully, ans is
c              accurate to within err * integral of abs(fun(x)).
c       ierr - a status code
c            - normal codes
c               1  ans most likely meets requested error tolerance.
c              -1  a equals b, or a and b are too nearly equal to
c                  allow normal integration.  ans is set to zero.
c            - abnormal code
c               2  ans probably does not meet requested error tolerance.
c       k    - the number of function evaluations actually used to do
c              the integration.  a value of k .gt. 1000 indicates a
c              difficult problem; other programs may be more efficient.
c              dqnc79 will gracefully give up if k exceeds 2000.
c
c***references  (none)
c***routines called  d1mach, i1mach, xermsg
c***revision history  (yymmdd)
c   790601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920218  code redone to parallel qnc79.  (wrb)
c   930120  increase array size 80->99, and kmx 2000->5000 for sun -r8
c           wordlength.  (rwc)
c***end prologue  dqnc79
c     .. scalar arguments ..
      double precision a, ans, b, err
      integer ierr, k
c     .. function arguments ..
      double precision fun
      external fun
c     .. local scalars ..
      double precision ae, area, bank, blocal, c, ce, ee, ef, eps, q13,
     +                 q7, q7l, sq2, test, tol, vr, w1, w2, w3, w4
      integer i, kml, kmx, l, lmn, lmx, nbits, nib, nlmn, nlmx
      logical first
c     .. local arrays ..
      double precision aa(99), f(13), f1(99), f2(99), f3(99), f4(99),
     +                 f5(99), f6(99), f7(99), hh(99), q7r(99), vl(99)
      integer lr(99)
c     .. external functions ..
      double precision d1mach
      integer i1mach
      external d1mach, i1mach
c     .. external subroutines ..
      external xermsg
c     .. intrinsic functions ..
      intrinsic abs, log, max, min, sign, sqrt
c     .. save statement ..
      save nbits, nlmx, first, sq2, w1, w2, w3, w4
c     .. data statements ..
      data kml /7/, kmx /5000/, nlmn /2/
      data first /.true./
c***first executable statement  dqnc79
      if (first) then
        w1 = 41.0d0/140.0d0
        w2 = 216.0d0/140.0d0
        w3 = 27.0d0/140.0d0
        w4 = 272.0d0/140.0d0
        nbits = d1mach(5)*i1mach(14)/0.30102000d0
        nlmx = min(99,(nbits*4)/5)
        sq2 = sqrt(2.0d0)
      endif
      first = .false.
      ans = 0.0d0
      ierr = 1
      ce = 0.0d0
      if (a .eq. b) go to 260
      lmx = nlmx
      lmn = nlmn
      if (b .eq. 0.0d0) go to 100
      if (sign(1.0d0,b)*a .le. 0.0d0) go to 100
      c = abs(1.0d0-a/b)
      if (c .gt. 0.1d0) go to 100
      if (c .le. 0.0d0) go to 260
      nib = 0.5d0 - log(c)/log(2.0d0)
      lmx = min(nlmx,nbits-nib-4)
      if (lmx .lt. 2) go to 260
      lmn = min(lmn,lmx)
  100 tol = max(abs(err),2.0d0**(5-nbits))
      if (err .eq. 0.0d0) tol = sqrt(d1mach(4))
      eps = tol
      hh(1) = (b-a)/12.0d0
      aa(1) = a
      lr(1) = 1
      do 110 i = 1,11,2
        f(i) = fun(a+(i-1)*hh(1))
  110 continue
      blocal = b
      f(13) = fun(blocal)
      k = 7
      l = 1
      area = 0.0d0
      q7 = 0.0d0
      ef = 256.0d0/255.0d0
      bank = 0.0d0
c
c     compute refined estimates, estimate the error, etc.
c
  120 do 130 i = 2,12,2
        f(i) = fun(aa(l)+(i-1)*hh(l))
  130 continue
      k = k + 6
c
c     compute left and right half estimates
c
      q7l = hh(l)*((w1*(f(1)+f(7))+w2*(f(2)+f(6)))+
     +      (w3*(f(3)+f(5))+w4*f(4)))
      q7r(l) = hh(l)*((w1*(f(7)+f(13))+w2*(f(8)+f(12)))+
     +         (w3*(f(9)+f(11))+w4*f(10)))
c
c     update estimate of integral of absolute value
c
      area = area + (abs(q7l)+abs(q7r(l))-abs(q7))
c
c     do not bother to test convergence before minimum refinement level
c
      if (l .lt. lmn) go to 180
c
c     estimate the error in new value for whole interval, q13
c
      q13 = q7l + q7r(l)
      ee = abs(q7-q13)*ef
c
c     compute nominal allowed error
c
      ae = eps*area
c
c     borrow from bank account, but not too much
c
      test = min(ae+0.8d0*bank,10.0d0*ae)
c
c     don't ask for excessive accuracy
c
      test = max(test,tol*abs(q13),0.00003d0*tol*area)
c
c     now, did this interval pass or not?
c
      if (ee-test) 150,150,170
c
c     have hit maximum refinement level -- penalize the cumulative error
c
  140 ce = ce + (q7-q13)
      go to 160
c
c     on good intervals accumulate the theoretical estimate
c
  150 ce = ce + (q7-q13)/255.0d0
c
c     update the bank account.  don't go into debt.
c
  160 bank = bank + (ae-ee)
      if (bank .lt. 0.0d0) bank = 0.0d0
c
c     did we just finish a left half or a right half?
c
      if (lr(l)) 190,190,210
c
c     consider the left half of next deeper level
c
  170 if (k .gt. kmx) lmx = min(kml,lmx)
      if (l .ge. lmx) go to 140
  180 l = l + 1
      eps = eps*0.5d0
      if (l .le. 17) ef = ef/sq2
      hh(l) = hh(l-1)*0.5d0
      lr(l) = -1
      aa(l) = aa(l-1)
      q7 = q7l
      f1(l) = f(7)
      f2(l) = f(8)
      f3(l) = f(9)
      f4(l) = f(10)
      f5(l) = f(11)
      f6(l) = f(12)
      f7(l) = f(13)
      f(13) = f(7)
      f(11) = f(6)
      f(9) = f(5)
      f(7) = f(4)
      f(5) = f(3)
      f(3) = f(2)
      go to 120
c
c     proceed to right half at this level
c
  190 vl(l) = q13
  200 q7 = q7r(l-1)
      lr(l) = 1
      aa(l) = aa(l) + 12.0d0*hh(l)
      f(1) = f1(l)
      f(3) = f2(l)
      f(5) = f3(l)
      f(7) = f4(l)
      f(9) = f5(l)
      f(11) = f6(l)
      f(13) = f7(l)
      go to 120
c
c     left and right halves are done, so go back up a level
c
  210 vr = q13
  220 if (l .le. 1) go to 250
      if (l .le. 17) ef = ef*sq2
      eps = eps*2.0d0
      l = l - 1
      if (lr(l)) 230,230,240
  230 vl(l) = vl(l+1) + vr
      go to 200
  240 vr = vl(l+1) + vr
      go to 220
c
c     exit
c
  250 ans = vr
      if (abs(ce) .le. 2.0d0*tol*area) go to 270
      ierr = 2
      call xermsg ('slatec', 'dqnc79',
     +   'ans is probably insufficiently accurate.', 2, 1)
      go to 270
  260 ierr = -1
      call xermsg ('slatec', 'dqnc79',
     +   'a and b are too nearly equal to allow normal integration. $$'
     +   // 'ans is set to zero and ierr to -1.', -1, -1)
  270 return
      end

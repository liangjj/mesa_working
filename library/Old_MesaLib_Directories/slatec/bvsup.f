*deck bvsup
      subroutine bvsup (y, nrowy, ncomp, xpts, nxpts, a, nrowa, alpha,
     +   nic, b, nrowb, beta, nfc, igofx, re, ae, iflag, work, ndw,
     +   iwork, ndiw, neqivp)
c***begin prologue  bvsup
c***purpose  solve a linear two-point boundary value problem using
c            superposition coupled with an orthonormalization procedure
c            and a variable-step integration scheme.
c***library   slatec
c***category  i1b1
c***type      single precision (bvsup-s, dbvsup-d)
c***keywords  orthonormalization, shooting,
c             two-point boundary value problem
c***author  scott, m. r., (snla)
c           watts, h. a., (snla)
c***description
c
c **********************************************************************
c     subroutine bvsup solves a linear two-point boundary-value problem
c     of the form
c                        dy/dx = matrix(x,u)*y(x) + g(x,u)
c                a*y(xinitial) = alpha ,  b*y(xfinal) = beta
c
c     coupled with the solution of the initial value problem
c
c                        du/dx = f(x,u)
c                      u(xinitial) = eta
c
c **********************************************************************
c     abstract
c        the method of solution uses superposition coupled with an
c     orthonormalization procedure and a variable-step integration
c     scheme.  each time the superposition solutions start to
c     lose their numerical linear independence, the vectors are
c     reorthonormalized before integration proceeds.  the underlying
c     principle of the algorithm is then to piece together the
c     intermediate (orthogonalized) solutions, defined on the various
c     subintervals, to obtain the desired solutions.
c
c **********************************************************************
c     input to bvsup
c **********************************************************************
c
c     nrowy = actual row dimension of y in calling program.
c             nrowy must be .ge. ncomp
c
c     ncomp = number of components per solution vector.
c             ncomp is equal to number of original differential
c             equations.  ncomp = nic + nfc.
c
c     xpts = desired output points for solution. they must be monotonic.
c            xinitial = xpts(1)
c            xfinal = xpts(nxpts)
c
c     nxpts = number of output points
c
c     a(nrowa,ncomp) = boundary condition matrix at xinitial,
c                      must be contained in (nic,ncomp) sub-matrix.
c
c     nrowa = actual row dimension of a in calling program,
c             nrowa must be .ge. nic.
c
c     alpha(nic+neqivp) = boundary conditions at xinitial.
c                         if neqivp .gt. 0 (see below), the boundary
c                         conditions at xinitial for the initial value
c                         equations must be stored starting in
c                         position (nic + 1) of alpha.
c                         thus,  alpha(nic+k) = eta(k).
c
c     nic = number of boundary conditions at xinitial.
c
c     b(nrowb,ncomp) = boundary condition matrix at xfinal,
c                      must be contained in (nfc,ncomp) sub-matrix.
c
c     nrowb = actual row dimension of b in calling program,
c             nrowb must be .ge. nfc.
c
c     beta(nfc) = boundary conditions at xfinal.
c
c     nfc = number of boundary conditions at xfinal
c
c     igofx =0 -- the inhomogeneous term g(x) is identically zero.
c           =1 -- the inhomogeneous term g(x) is not identically zero.
c                 (if igofx=1, then subroutine gvec (or uvec) must be
c                  supplied).
c
c     re = relative error tolerance used by the integrator
c          (see one of the integrators)
c
c     ae = absolute error tolerance used by the integrator
c          (see one of the integrators)
c **note-  re and ae should not both be zero.
c
c     iflag = a status parameter used principally for output.
c             however, for efficient solution of problems which
c             are originally defined as complex valued (but
c             converted to real systems to use this code), the
c             user must set iflag=13 on input. see the comment below
c             for more information on solving such problems.
c
c     work(ndw) = floating point array used for internal storage.
c
c     ndw = actual dimension of work array allocated by user.
c           an estimate for ndw can be computed from the following
c            ndw = 130 + ncomp**2 * (6 + nxpts/2 + expected number of
c                                                orthonormalizations/8)
c             for the disk or tape storage mode,
c            ndw = 6 * ncomp**2 + 10 * ncomp + 130
c  however, when the adams integrator is to be used, the estimates are
c            ndw = 130 + ncomp**2 * (13 + nxpts/2 + expected number of
c                                                orthonormalizations/8)
c    and     ndw = 13 * ncomp**2 + 22 * ncomp + 130   , respectively.
c
c     iwork(ndiw) = integer array used for internal storage.
c
c     ndiw = actual dimension of iwork array allocated by user.
c            an estimate for ndiw can be computed from the following
c            ndiw = 68 + ncomp * (1 + expected number of
c                                        orthonormalizations)
c **note --  the amount of storage required is problem dependent and may
c            be difficult to predict in advance. experience has shown
c            that for most problems 20 or fewer orthonormalizations
c            should suffice. if the problem cannot be completed with the
c            allotted storage, then a message will be printed which
c            estimates the amount of storage necessary. in any case, the
c            user can examine the iwork array for the actual storage
c            requirements, as described in the output information below.
c
c     neqivp = number of auxiliary initial value equations being added
c              to the boundary value problem.
c **note -- occasionally the coefficients  matrix  and/or  g  may be
c           functions which depend on the independent variable  x  and
c           on  u, the solution of an auxiliary initial value problem.
c           in order to avoid the difficulties associated with
c           interpolation, the auxiliary equations may be solved
c           simultaneously with the given boundary value problem.
c           this initial value problem may be linear or nonlinear.
c                 see sand77-1328 for an example.
c
c
c     the user must supply subroutines fmat, gvec, uivp and uvec, when
c     needed (they must be so named), to evaluate the derivatives
c     as follows
c
c        a. fmat must be supplied.
c
c              subroutine fmat(x,y,yp)
c              x = independent variable (input to fmat)
c              y = dependent variable vector (input to fmat)
c              yp = dy/dx = derivative vector (output from fmat)
c
c            compute the derivatives for the homogeneous problem
c              yp(i) = dy(i)/dx = matrix(x) * y(i)  , i = 1,...,ncomp
c
c            when (neqivp .gt. 0) and  matrix  is dependent on  u  as
c            well as on  x, the following common statement must be
c            included in fmat
c                    common /mlivp/ nofst
c            for convenience, the  u  vector is stored at the bottom
c            of the  y  array.  thus, during any call to fmat,
c            u(i) is referenced by  y(nofst + i).
c
c
c            subroutine bvder calls fmat nfc times to evaluate the
c            homogeneous equations and, if necessary, it calls fmat once
c            in evaluating the particular solution. since x remains
c            unchanged in this sequence of calls it is possible to
c            realize considerable computational savings for complicated
c            and expensive evaluations of the matrix entries. to do this
c            the user merely passes a variable, say xs, via common where
c            xs is defined in the main program to be any value except
c            the initial x. then the non-constant elements of matrix(x)
c            appearing in the differential equations need only be
c            computed if x is unequal to xs, whereupon xs is reset to x.
c
c
c        b. if  neqivp .gt. 0 ,  uivp must also be supplied.
c
c              subroutine uivp(x,u,up)
c              x = independent variable (input to uivp)
c              u = dependent variable vector (input to uivp)
c              up = du/dx = derivative vector (output from uivp)
c
c            compute the derivatives for the auxiliary initial value eqs
c              up(i) = du(i)/dx, i = 1,...,neqivp.
c
c            subroutine bvder calls uivp once to evaluate the
c            derivatives for the auxiliary initial value equations.
c
c
c        c. if  neqivp = 0  and  igofx = 1 ,  gvec must be supplied.
c
c              subroutine gvec(x,g)
c              x = independent variable (input to gvec)
c              g = vector of inhomogeneous terms g(x) (output from gvec)
c
c            compute the inhomogeneous terms g(x)
c                g(i) = g(x) values for i = 1,...,ncomp.
c
c            subroutine bvder calls gvec in evaluating the particular
c            solution provided g(x) is not identically zero. thus, when
c            igofx=0, the user need not write a gvec subroutine. also,
c            the user does not have to bother with the computational
c            savings scheme for gvec as this is automatically achieved
c            via the bvder subroutine.
c
c
c        d. if  neqivp .gt. 0  and  igofx = 1 ,  uvec must be supplied.
c
c              subroutine uvec(x,u,g)
c              x = independent variable (input to uvec)
c              u = dependent variable vector from the auxiliary initial
c                  value problem    (input to uvec)
c              g = array of inhomogeneous terms g(x,u)(output from uvec)
c
c            compute the inhomogeneous terms g(x,u)
c                g(i) = g(x,u) values for i = 1,...,ncomp.
c
c            subroutine bvder calls uvec in evaluating the particular
c            solution provided g(x,u) is not identically zero.  thus,
c            when igofx=0, the user need not write a uvec subroutine.
c
c
c
c     the following is optional input to bvsup to give the user more
c     flexibility in use of the code.  see sand75-0198 , sand77-1328 ,
c     sand77-1690,sand78-0522, and sand78-1501 for more information.
c
c ****caution -- the user must zero out iwork(1),...,iwork(15)
c                prior to calling bvsup. these locations define optional
c                input and must be zero unless set to special values by
c                the user as described below.
c
c     iwork(1) -- number of orthonormalization points.
c                 a value need be set only if iwork(11) = 1
c
c     iwork(9) -- integrator and orthonormalization parameter
c                 (default value is 1)
c                 1 = runge-kutta-fehlberg code using gram-schmidt test.
c                 2 = adams code using gram-schmidt test.
c
c     iwork(11) -- orthonormalization points parameter
c                  (default value is 0)
c                  0 - orthonormalization points not pre-assigned.
c                  1 - orthonormalization points pre-assigned in
c                      the first iwork(1) positions of work.
c
c     iwork(12) -- storage parameter
c                  (default value is 0)
c                  0 - all storage in core
c                lun - homogeneous and inhomogeneous solutions at
c                     output points and orthonormalization information
c                     are stored on disk.  the logical unit number to be
c                     used for disk i/o (ntape) is set to iwork(12).
c
c     work(1),... -- pre-assigned orthonormalization points, stored
c                    monotonically, corresponding to the direction
c                    of integration.
c
c
c
c                 ******************************
c                 *** complex valued problem ***
c                 ******************************
c **note***
c       suppose the original boundary value problem is nc equations
c     of the form
c                   dw/dx = mat(x,u)*w(x) + h(x,u)
c                 r*w(xinitial)=gamma , s*w(xfinal)=delta
c
c     where all variables are complex valued. the bvsup code can be
c     used by converting to a real system of size 2*nc. to solve the
c     larger dimensioned problem efficiently,  the user must initialize
c     iflag=13 on input and order the vector components according to
c     y(1)=real(w(1)),...,y(nc)=real(w(nc)),y(nc+1)=imag(w(1)),....,
c     y(2*nc)=imag(w(nc)). then define
c                        ...........................
c                        . real(mat)    -imag(mat) .
c            matrix  =   .                         .
c                        . imag(mat)     real(mat) .
c                        ...........................
c
c     the matrices a,b and vectors g,alpha,beta must be defined
c     similarly. further details can be found in sand78-1501.
c
c
c **********************************************************************
c     output from bvsup
c **********************************************************************
c
c     y(nrowy,nxpts) = solution at specified output points.
c
c     iflag output values
c            =-5 algorithm ,for obtaining starting vectors for the
c                special complex problem structure, was unable to obtain
c                the initial vectors satisfying the necessary
c                independence criteria.
c            =-4 rank of boundary condition matrix a is less than nic,
c                as determined by lssuds.
c            =-2 invalid input parameters.
c            =-1 insufficient number of storage locations allocated for
c                work or iwork.
c
c            =0 indicates successful solution
c
c            =1 a computed solution is returned but uniqueness of the
c               solution of the boundary-value problem is questionable.
c               for an eigenvalue problem, this should be treated as a
c               successful execution since this is the expected mode
c               of return.
c            =2 a computed solution is returned but the existence of the
c               solution to the boundary-value problem is questionable.
c            =3 a nontrivial solution approximation is returned although
c               the boundary condition matrix b*y(xfinal) is found to be
c               nonsingular (to the desired accuracy level) while the
c               right hand side vector is zero. to eliminate this type
c               of return, the accuracy of the eigenvalue parameter
c               must be improved.
c           ***note- we attempt to diagnose the correct problem behavior
c               and report possible difficulties by the appropriate
c               error flag.  however, the user should probably resolve
c               the problem using smaller error tolerances and/or
c               perturbations in the boundary conditions or other
c               parameters. this will often reveal the correct
c               interpretation for the problem posed.
c
c            =13 maximum number of orthonormalizations attained before
c                reaching xfinal.
c            =20-flag from integrator (derkf or deabm) values can range
c                from 21 to 25.
c            =30 solution vectors form a dependent set.
c
c     work(1),...,work(iwork(1)) = orthonormalization points
c                                  determined by bvpor.
c
c     iwork(1) = number of orthonormalizations performed by bvpor.
c
c     iwork(2) = maximum number of orthonormalizations allowed as
c                calculated from storage allocated by user.
c
c     iwork(3),iwork(4),iwork(5),iwork(6)   give information about
c                actual storage requirements for work and iwork
c                arrays.  in particular,
c                       required storage for  work array is
c        iwork(3) + iwork(4)*(expected number of orthonormalizations)
c
c                       required storage for iwork array is
c        iwork(5) + iwork(6)*(expected number of orthonormalizations)
c
c     iwork(8) = final value of exponent parameter used in tolerance
c                test for orthonormalization.
c
c     iwork(16) = number of independent vectors returned from mgsbv.
c                 it is only of interest when iflag=30 is obtained.
c
c     iwork(17) = numerically estimated rank of the boundary
c                 condition matrix defined from b*y(xfinal)
c
c **********************************************************************
c
c     necessary machine constants are defined in the function
c     routine r1mach. the user must make sure that the values
c     set in r1mach are relevant to the computer being used.
c
c **********************************************************************
c
c***references  m. r. scott and h. a. watts, suport - a computer code
c                 for two-point boundary-value problems via
c                 orthonormalization, siam journal of numerical
c                 analysis 14, (1977), pp. 40-70.
c               b. l. darlow, m. r. scott and h. a. watts, modifications
c                 of suport, a linear boundary value problem solver
c                 part i - pre-assigning orthonormalization points,
c                 auxiliary initial value problem, disk or tape storage,
c                 report sand77-1328, sandia laboratories, albuquerque,
c                 new mexico, 1977.
c               b. l. darlow, m. r. scott and h. a. watts, modifications
c                 of suport, a linear boundary value problem solver
c                 part ii - inclusion of an adams integrator, report
c                 sand77-1690, sandia laboratories, albuquerque,
c                 new mexico, 1977.
c               m. e. lord and h. a. watts, modifications of suport,
c                 a linear boundary value problem solver part iii -
c                 orthonormalization improvements, report sand78-0522,
c                 sandia laboratories, albuquerque, new mexico, 1978.
c               h. a. watts, m. r. scott and m. e. lord, computational
c                 solution of complex*16 valued boundary problems,
c                 report sand78-1501, sandia laboratories,
c                 albuquerque, new mexico, 1978.
c***routines called  exbvp, macon, xermsg
c***common blocks    ml15to, ml17bw, ml18jr, ml5mco, ml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   890921  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bvsup
c **********************************************************************
c
c
      dimension y(nrowy,*),a(nrowa,*),alpha(*),b(nrowb,*),
     1          beta(*),work(*),iwork(*),xpts(*)
      character*8 xern1, xern2, xern3, xern4
c
c **********************************************************************
c     the common block below is used to communicate with subroutine
c     bvder.  the user should not alter or use this common block in the
c     calling program.
c
      common /ml8sz/ c,xsav,igofxd,inhomo,ivp,ncompd,nfcd
c
c **********************************************************************
c     these common blocks aid in reducing the number of subroutine
c     arguments prevalent in this modular structure
c
      common /ml18jr/ aed,red,tol,nxptsd,nicd,nopg,mxnon,ndisk,ntape,
     1                neq,indpvt,integ,nps,ntp,neqivd,numort,nfcc,
     2                icoco
      common /ml17bw/ kkkzpw,needw,neediw,k1,k2,k3,k4,k5,k6,k7,k8,k9,
     1                k10,k11,l1,l2,kkkint,lllint
c
c **********************************************************************
c     this common block is used in subroutines bvsup,bvpor,rkfab,
c     reort, and stway. it contains information necessary
c     for the orthonormalization testing procedure and a backup
c     restarting capability.
c
      common /ml15to/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
c
c **********************************************************************
c     this common block contains the machine dependent parameters
c     used by the code
c
      common /ml5mco/ uro,sru,eps,sqovfl,twou,fouru,lpar
c
c **********************************************************************
c     set up machine dependent constants.
c
c***first executable statement  bvsup
      call macon
c
c **********************************************************************
c     test for invalid input
c
      if (nrowy .lt. ncomp)  go to 20
      if (ncomp .ne. nic+nfc)  go to 20
      if (nxpts .lt. 2)  go to 20
      if (nic .le. 0)  go to 20
      if (nrowa .lt. nic)  go to 20
      if (nfc .le. 0)  go to 20
      if (nrowb .lt. nfc)  go to 20
      if (igofx .lt. 0  .or.  igofx .gt. 1) go to 20
      if (re .lt. 0.0)  go to 20
      if (ae .lt. 0.0)  go to 20
      if (re .eq. 0.0  .and.  ae .eq. 0.0)  go to 20
      is = 1
      if (xpts(nxpts) .lt. xpts(1))  is = 2
      nxptsm = nxpts - 1
      do 13 k = 1,nxptsm
      if (is .eq. 2) go to 12
      if (xpts(k+1) .le. xpts(k))  go to 20
      go to 13
   12 if (xpts(k) .le. xpts(k+1))  go to 20
   13 continue
      go to 30
   20 iflag = -2
      return
   30 continue
c
c **********************************************************************
c     check for disk storage
c
      kpts = nxpts
      ndisk = 0
      if (iwork(12) .eq. 0)  go to 35
      ntape = iwork(12)
      kpts = 1
      ndisk = 1
   35 continue
c
c **********************************************************************
c     set integ parameter according to choice of integrator.
c
      integ = 1
      if (iwork(9) .eq. 2)  integ = 2
c
c **********************************************************************
c     compute inhomo
c
      if (igofx .eq. 1)  go to 43
      do 40 j = 1,nic
      if (alpha(j) .ne. 0.0)  go to 43
   40 continue
      do 41 j = 1,nfc
      if (beta(j) .ne. 0.0)  go to 42
   41 continue
      inhomo = 3
      go to 45
   42 inhomo = 2
      go to 45
   43 inhomo = 1
   45 continue
c
c **********************************************************************
c     to take advantage of the special structure when solving a
c     complex valued problem,we introduce nfcc=nfc while changing
c     the internal value of nfc
c
      nfcc=nfc
      if (iflag .eq. 13) nfc=nfc/2
c
c **********************************************************************
c     determine necessary storage requirements
c
c for basic arrays in bvpor
      kkkyhp = ncomp*(nfc+1) + neqivp
      kkku   = ncomp*nfc*kpts
      kkkv   = ncomp*kpts
      kkkcoe = nfcc
      kkks   = nfc+1
      kkksto = ncomp*(nfc+1) + neqivp + 1
      kkkg   = ncomp
c
c for orthonormalization related matters
      ntp = (nfcc*(nfcc+1))/2
      kkkzpw = 1 + ntp + nfcc
      lllip  = nfcc
c
c for additional required work space
c   (lssuds)
      kkksud = 4*nic + (nrowa+1)*ncomp
      lllsud = nic
c   (svecs)
      kkksvc = 1 + 4*nfcc + 2*nfcc**2
      lllsvc = 2*nfcc
c
      ndeq=ncomp*nfc+neqivp
      if (inhomo .eq. 1) ndeq=ndeq+ncomp
      go to (51,52),integ
c   (derkf)
   51 kkkint = 33 + 7*ndeq
      lllint = 34
      go to 55
c   (deabm)
   52 kkkint = 130 + 21*ndeq
      lllint = 51
c
c   (coef)
   55 kkkcof = 5*nfcc + nfcc**2
      lllcof = 3 + nfcc
c
      kkkws  = max(kkksud,kkksvc,kkkint,kkkcof)
      llliws = max(lllsud,lllsvc,lllint,lllcof)
c
      needw  = kkkyhp + kkku + kkkv + kkkcoe + kkks + kkksto + kkkg +
     1         kkkzpw + kkkws
      neediw = 17 + lllip + llliws
c **********************************************************************
c     compute the number of possible orthonormalizations with the
c     allotted storage
c
      iwork(3) = needw
      iwork(4) = kkkzpw
      iwork(5) = neediw
      iwork(6) = lllip
      nrtemp = ndw - needw
      nitemp = ndiw - neediw
      if (nrtemp .lt. 0)  go to 70
      if (nitemp .ge. 0)  go to 75
c
   70 iflag = -1
      if (ndisk .ne. 1) then
         write (xern1, '(i8)') needw
         write (xern2, '(i8)') kkkzpw
         write (xern3, '(i8)') neediw
         write (xern4, '(i8)') lllip
         call xermsg ('slatec', 'bvsup',
     *      'required storage for work array is '  // xern1 // ' + ' //
     *      xern2 // '*(expected number of orthonormalizations) $$'  //
     *      'required storage for iwork array is ' // xern3 // ' + ' //
     *      xern4 // '*(expected number of orthonormalizations)', 1, 0)
      else
         write (xern1, '(i8)') needw
         write (xern2, '(i8)') neediw
         call xermsg ('slatec', 'bvsup',
     *      'required storage for work array is '  // xern1 //
     *      ' + number of orthonomalizations. $$'  //
     *      'required storage for iwork array is ' // xern2, 1, 0)
      endif
      return
c
   75 if (ndisk .eq. 0)  go to 77
      non = 0
      mxnon = nrtemp
      go to 78
c
   77 mxnonr = nrtemp / kkkzpw
      mxnoni = nitemp / lllip
      mxnon = min(mxnonr,mxnoni)
      non = mxnon
c
   78 iwork(2) = mxnon
c
c **********************************************************************
c     check for pre-assigned orthonormalization points
c
      nopg = 0
      if (iwork(11) .ne. 1)  go to 85
      if (mxnon .lt. iwork(1))  go to 70
      nopg = 1
      mxnon = iwork(1)
      work(mxnon+1) = 2. * xpts(nxpts)  -  xpts(1)
   85 continue
c
c **********************************************************************
c     allocate storage from work and iwork arrays
c
c  (z)
      k1 = 1 + (mxnon+1)
c  (p)
      k2 = k1 + ntp*(non+1)
c  (w)
      k3 = k2 + nfcc*(non+1)
c  (yhp)
      k4 = k3 + kkkyhp
c  (u)
      k5 = k4 + kkku
c  (v)
      k6 = k5 + kkkv
c  (coef)
      k7 = k6 + kkkcoe
c  (s)
      k8 = k7 + kkks
c  (stowa)
      k9 = k8 + kkksto
c  (g)
      k10 = k9 + kkkg
      k11 = k10 + kkkws
c            required additional real work space starts at work(k10)
c            and extends to work(k11-1)
c
c     first 17 locations of iwork are used for optional
c     input and output items
c  (ip)
      l1 = 18 + nfcc*(non+1)
      l2 = l1 + llliws
c            required integer work space starts at iwork(l1)
c            and extends to iwork(l2-1)
c
c **********************************************************************
c     set indicator for normalization of particular solution
c
      nps = 0
      if (iwork(10) .eq. 1)  nps = 1
c
c **********************************************************************
c     set pivoting parameter
c
      indpvt=0
      if (iwork(15) .eq. 1) indpvt=1
c
c **********************************************************************
c     set other common block parameters
c
      nfcd = nfc
      ncompd = ncomp
      igofxd = igofx
      nxptsd = nxpts
      nicd = nic
      red = re
      aed = ae
      neqivd = neqivp
      mnswot = 20
      if (iwork(13) .eq. -1) mnswot=max(1,iwork(14))
      xbeg=xpts(1)
      xend=xpts(nxpts)
      xsav=xend
      icoco=1
      if (inhomo .eq. 3  .and.  nopg .eq. 1) work(mxnon+1)=xend
c
c **********************************************************************
c
      call exbvp(y,nrowy,xpts,a,nrowa,alpha,b,nrowb,beta,iflag,work,
     1           iwork)
      nfc=nfcc
      iwork(17)=iwork(l1)
      return
      end

*deck dbvsup
      subroutine dbvsup (y, nrowy, ncomp, xpts, nxpts, a, nrowa, alpha,
     +   nic, b, nrowb, beta, nfc, igofx, re, ae, iflag, work, ndw,
     +   iwork, ndiw, neqivp)
c***begin prologue  dbvsup
c***purpose  solve a linear two-point boundary value problem using
c            superposition coupled with an orthonormalization procedure
c            and a variable-step integration scheme.
c***library   slatec
c***category  i1b1
c***type      double precision (bvsup-s, dbvsup-d)
c***keywords  orthonormalization, shooting,
c             two-point boundary value problem
c***author  scott, m. r., (snla)
c           watts, h. a., (snla)
c***description
c
c **********************************************************************
c
c     subroutine dbvsup solves a linear two-point boundary-value problem
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
c     input to dbvsup
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
c     nxpts = number of output points.
c
c     a(nrowa,ncomp) = boundary condition matrix at xinitial
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
c     b(nrowb,ncomp) = boundary condition matrix at xfinal.
c                      must be contained in (nfc,ncomp) sub-matrix.
c
c     nrowb = actual row dimension of b in calling program,
c             nrowb must be .ge. nfc.
c
c     beta(nfc) = boundary conditions at xfinal.
c
c     nfc = number of boundary conditions at xfinal.
c
c     igofx =0 -- the inhomogeneous term g(x) is identically zero.
c           =1 -- the inhomogeneous term g(x) is not identically zero.
c                 (if igofx=1, then subroutine dgvec (or duvec) must be
c                  supplied).
c
c     re = relative error tolerance used by the integrator.
c          (see one of the integrators)
c
c     ae = absolute error tolerance used by the integrator.
c          (see one of the integrators)
c **note-  re and ae should not both be zero.
c
c     iflag = a status parameter used principally for output.
c             however, for efficient solution of problems which
c             are originally defined as complex*16 valued (but
c             converted to double precision systems to use this code),
c             the user must set iflag=13 on input. see the comment
c             below for more information on solving such problems.
c
c     work(ndw) = floating point array used for internal storage.
c
c     ndw = actual dimension of work array allocated by user.
c           an estimate for ndw can be computed from the following
c            ndw = 130 + ncomp**2 * (6 + nxpts/2 + expected number of
c                                           orthonormalizations/8)
c           for the disk or tape storage mode,
c            ndw = 6 * ncomp**2 + 10 * ncomp + 130
c  however, when the adams integrator is to be used, the estimates are
c            ndw = 130 + ncomp**2 * (13 + nxpts/2 + expected number of
c                                           orthonormalizations/8)
c    and     ndw = 13 * ncomp**2 + 22 * ncomp + 130   , respectively.
c
c     iwork(ndiw) = integer array used for internal storage.
c
c     ndiw = actual dimension of iwork array allocated by user.
c            an estimate for ndiw can be computed from the following
c            ndiw = 68 + ncomp * (1 + expected number of
c                                            orthonormalizations)
c **note --  the amount of storage required is problem dependent and may
c            be difficult to predict in advance.  experience has shown
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
c     the user must supply subroutines dfmat, dgvec, duivp and duvec,
c     when needed (they must be so named), to evaluate the derivatives
c     as follows
c
c        a. dfmat must be supplied.
c
c              subroutine dfmat(x,y,yp)
c              x = independent variable (input to dfmat)
c              y = dependent variable vector (input to dfmat)
c              yp = dy/dx = derivative vector (output from dfmat)
c
c            compute the derivatives for the homogeneous problem
c              yp(i) = dy(i)/dx = matrix(x) * y(i)  , i = 1,...,ncomp
c
c            when (neqivp .gt. 0) and  matrix  is dependent on  u  as
c            well as on  x, the following common statement must be
c            included in dfmat
c                    common /dmlivp/ nofst
c            for convenience, the  u  vector is stored at the bottom
c            of the  y  array.  thus, during any call to dfmat,
c            u(i) is referenced by  y(nofst + i).
c
c
c            subroutine dbvder calls dfmat nfc times to evaluate the
c            homogeneous equations and, if necessary, it calls dfmat
c            once in evaluating the particular solution. since x remains
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
c        b. if  neqivp .gt. 0 ,  duivp must also be supplied.
c
c              subroutine duivp(x,u,up)
c              x = independent variable (input to duivp)
c              u = dependent variable vector (input to duivp)
c              up = du/dx = derivative vector (output from duivp)
c
c            compute the derivatives for the auxiliary initial value eqs
c              up(i) = du(i)/dx, i = 1,...,neqivp.
c
c            subroutine dbvder calls duivp once to evaluate the
c            derivatives for the auxiliary initial value equations.
c
c
c        c. if  neqivp = 0  and  igofx = 1 ,  dgvec must be supplied.
c
c              subroutine dgvec(x,g)
c              x = independent variable (input to dgvec)
c              g = vector of inhomogeneous terms g(x) (output from
c              dgvec)
c
c            compute the inhomogeneous terms g(x)
c                g(i) = g(x) values for i = 1,...,ncomp.
c
c            subroutine dbvder calls dgvec in evaluating the particular
c            solution provided g(x) is not identically zero. thus, when
c            igofx=0, the user need not write a dgvec subroutine. also,
c            the user does not have to bother with the computational
c            savings scheme for dgvec as this is automatically achieved
c            via the dbvder subroutine.
c
c
c        d. if  neqivp .gt. 0  and  igofx = 1 ,  duvec must be supplied.
c
c             subroutine duvec(x,u,g)
c             x = independent variable (input to duvec)
c             u = dependent variable vector from the auxiliary initial
c                 value problem    (input to duvec)
c             g = array of inhomogeneous terms g(x,u)(output from duvec)
c
c            compute the inhomogeneous terms g(x,u)
c                g(i) = g(x,u) values for i = 1,...,ncomp.
c
c            subroutine dbvder calls duvec in evaluating the particular
c            solution provided g(x,u) is not identically zero.  thus,
c            when igofx=0, the user need not write a duvec subroutine.
c
c
c
c     the following is optional input to dbvsup to give user more
c     flexibility in use of code.  see sand75-0198, sand77-1328,
c     sand77-1690, sand78-0522, and sand78-1501 for more information.
c
c ****caution -- the user must zero out iwork(1),...,iwork(15)
c                prior to calling dbvsup. these locations define
c                optional input and must be zero unless set to special
c                values by the user as described below.
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
c                  0 - all storage in core.
c                  lun - homogeneous and inhomogeneous solutions at
c                      output points and orthonormalization information
c                      are stored on disk.  the logical unit number to
c                      be used for disk i/o (ntape) is set to iwork(12).
c
c     work(1),... -- pre-assigned orthonormalization points, stored
c                    monotonically, corresponding to the direction
c                    of integration.
c
c
c
c                 ******************************************************
c                 *** complex*16 valued problem ***
c                 ******************************************************
c **note***
c       suppose the original boundary value problem is nc equations
c     of the form
c                   dw/dx = mat(x,u)*w(x) + h(x,u)
c                 r*w(xinitial)=gamma , s*w(xfinal)=delta
c     where all variables are complex*16 valued. the dbvsup code can be
c     used by converting to a double precision system of size 2*nc. to
c     solve the larger dimensioned problem efficiently, the user must
c     initialize iflag=13 on input and order the vector components
c     according to y(1)=double precision(w(1)),...,y(nc)=double
c     precision(w(nc)),y(nc+1)=imag(w(1)),...., y(2*nc)=imag(w(nc)).
c     then define
c                        ...............................................
c                        . double precision(mat)    -imag(mat) .
c            matrix  =   .                         .
c                        . imag(mat)     double precision(mat) .
c                        ...............................................
c
c     the matrices a,b and vectors g,alpha,beta must be defined
c     similarly. further details can be found in sand78-1501.
c
c
c **********************************************************************
c     output from dbvsup
c **********************************************************************
c
c     y(nrowy,nxpts) = solution at specified output points.
c
c     iflag output values
c            =-5 algorithm ,for obtaining starting vectors for the
c                special complex*16 problem structure, was unable to
c                obtain the initial vectors satisfying the necessary
c                independence criteria.
c            =-4 rank of boundary condition matrix a is less than nic,
c                as determined by dlssud.
c            =-2 invalid input parameters.
c            =-1 insufficient number of storage locations allocated for
c                work or iwork.
c
c            =0 indicates successful solution.
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
c            ***note-we attempt to diagnose the correct problem behavior
c               and report possible difficulties by the appropriate
c               error flag.  however, the user should probably resolve
c               the problem using smaller error tolerances and/or
c               perturbations in the boundary conditions or other
c               parameters. this will often reveal the correct
c               interpretation for the problem posed.
c
c            =13 maximum number of orthonormalizations attained before
c                reaching xfinal.
c            =20-flag from integrator (dderkf or ddeabm) values can
c                range from 21 to 25.
c            =30 solution vectors form a dependent set.
c
c     work(1),...,work(iwork(1)) = orthonormalization points
c                                  determined by dbvpor.
c
c     iwork(1) = number of orthonormalizations performed by dbvpor.
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
c     iwork(16) = number of independent vectors returned from dmgsbv.
c                it is only of interest when iflag=30 is obtained.
c
c     iwork(17) = numerically estimated rank of the boundary
c                 condition matrix defined from b*y(xfinal)
c
c **********************************************************************
c
c     necessary machine constants are defined in the function
c     routine d1mach. the user must make sure that the values
c     set in d1mach are relevant to the computer being used.
c
c **********************************************************************
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
c***routines called  dexbvp, dmacon, xermsg
c***common blocks    dml15t, dml17b, dml18j, dml5mc, dml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   890921  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  convert xerrwv calls to xermsg calls, remove some extraneous
c           comments.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbvsup
c **********************************************************************
c
      integer icoco, iflag, igofx, igofxd, indpvt, info, inhomo, integ,
     1     is, istkop, ivp, iwork(*), j, k, k1, k10, k11, k2,
     2     k3, k4, k5, k6, k7, k8, k9, kkkcoe, kkkcof, kkkg, kkkint,
     3     kkks, kkksto, kkksud, kkksvc, kkku, kkkv, kkkws, kkkyhp,
     4     kkkzpw, knswot, kop, kpts, l1, l2, lllcof, lllint, lllip,
     5     llliws, lllsud, lllsvc, lotjp, lpar, mnswot,
     6     mxnon, mxnoni, mxnonr, ncomp, ncompd, ndeq, ndisk, ndiw,
     7     ndw, neediw, needw, neq, neqivd, neqivp, nfc, nfcc,
     8     nfcd, nic, nicd, nitemp, non, nopg, nps, nrowa, nrowb,
     9     nrowy, nrtemp, nswot, ntape, ntp, numort, nxpts, nxptsd,
     1     nxptsm
      double precision a(nrowa,*), ae, aed, alpha(*),
     1     b(nrowb,*), beta(*), c, eps, fouru, pwcnd, px, re,
     2     red, sqovfl, sru, tnd, tol, twou, uro, work(ndw), x, xbeg,
     3     xend, xop, xot, xpts(*), xsav, y(nrowy,*)
      character*8 xern1, xern2, xern3, xern4
c
c     ******************************************************************
c         the common block below is used to communicate with subroutine
c         dbvder.  the user should not alter or use this common block in
c         the calling program.
c
      common /dml8sz/ c,xsav,igofxd,inhomo,ivp,ncompd,nfcd
c
c     ******************************************************************
c         these common blocks aid in reducing the number of subroutine
c         arguments prevalent in this modular structure
c
      common /dml18j/ aed,red,tol,nxptsd,nicd,nopg,mxnon,ndisk,ntape,
     1                neq,indpvt,integ,nps,ntp,neqivd,numort,nfcc,
     2                icoco
      common /dml17b/ kkkzpw,needw,neediw,k1,k2,k3,k4,k5,k6,k7,k8,k9,
     1                k10,k11,l1,l2,kkkint,lllint
c
c     ******************************************************************
c         this common block is used in subroutines dbvsup,dbvpor,drkfab,
c         dreort, and dstway. it contains information necessary
c         for the orthonormalization testing procedure and a backup
c         restarting capability.
c
      common /dml15t/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
c
c     ******************************************************************
c         this common block contains the machine dependent parameters
c         used by the code
c
      common /dml5mc/ uro,sru,eps,sqovfl,twou,fouru,lpar
c
c      *****************************************************************
c          set up machine dependent constants.
c
c***first executable statement  dbvsup
                        call dmacon
c
c                       ************************************************
c                           test for invalid input
c
                        if (nrowy .lt. ncomp) go to 80
                        if (ncomp .ne. nic + nfc) go to 80
                        if (nxpts .lt. 2) go to 80
                        if (nic .le. 0) go to 80
                        if (nrowa .lt. nic) go to 80
                        if (nfc .le. 0) go to 80
                        if (nrowb .lt. nfc) go to 80
                        if (igofx .lt. 0 .or. igofx .gt. 1) go to 80
                        if (re .lt. 0.0d0) go to 80
                        if (ae .lt. 0.0d0) go to 80
                        if (re .eq. 0.0d0 .and. ae .eq. 0.0d0) go to 80
c                          begin block permitting ...exits to 70
                              is = 1
                              if (xpts(nxpts) .lt. xpts(1)) is = 2
                              nxptsm = nxpts - 1
                              do 30 k = 1, nxptsm
                                 if (is .eq. 2) go to 10
c                          .........exit
                                    if (xpts(k+1) .le. xpts(k)) go to 70
                                 go to 20
   10                            continue
c                          .........exit
                                    if (xpts(k) .le. xpts(k+1)) go to 70
   20                            continue
   30                         continue
c
c                             ******************************************
c                                 check for disk storage
c
                              kpts = nxpts
                              ndisk = 0
                              if (iwork(12) .eq. 0) go to 40
                                 ntape = iwork(12)
                                 kpts = 1
                                 ndisk = 1
   40                         continue
c
c                             ******************************************
c                                 set integ parameter according to
c                                 choice of integrator.
c
                              integ = 1
                              if (iwork(9) .eq. 2) integ = 2
c
c                             ******************************************
c                                 compute inhomo
c
c                 ............exit
                              if (igofx .eq. 1) go to 100
                              do 50 j = 1, nic
c                 ...............exit
                                 if (alpha(j) .ne. 0.0d0) go to 100
   50                         continue
                              do 60 j = 1, nfc
c                    ............exit
                                 if (beta(j) .ne. 0.0d0) go to 90
   60                         continue
                              inhomo = 3
c              ...............exit
                              go to 110
   70                      continue
   80                   continue
                        iflag = -2
c     ..................exit
                        go to 220
   90                continue
                     inhomo = 2
c              ......exit
                     go to 110
  100             continue
                  inhomo = 1
  110          continue
c
c              *********************************************************
c                  to take advantage of the special structure when
c                  solving a complex*16 valued problem,we introduce
c                  nfcc=nfc while changing the internal value of nfc
c
               nfcc = nfc
               if (iflag .eq. 13) nfc = nfc/2
c
c              *********************************************************
c                  determine necessary storage requirements
c
c              for basic arrays in dbvpor
               kkkyhp = ncomp*(nfc + 1) + neqivp
               kkku = ncomp*nfc*kpts
               kkkv = ncomp*kpts
               kkkcoe = nfcc
               kkks = nfc + 1
               kkksto = ncomp*(nfc + 1) + neqivp + 1
               kkkg = ncomp
c
c              for orthonormalization related matters
               ntp = (nfcc*(nfcc + 1))/2
               kkkzpw = 1 + ntp + nfcc
               lllip = nfcc
c
c              for additional required work space
c                (dlssud)
               kkksud = 4*nic + (nrowa + 1)*ncomp
               lllsud = nic
c              (dvecs)
               kkksvc = 1 + 4*nfcc + 2*nfcc**2
               lllsvc = 2*nfcc
c
               ndeq = ncomp*nfc + neqivp
               if (inhomo .eq. 1) ndeq = ndeq + ncomp
               go to (120,130), integ
c              (dderkf)
  120          continue
                  kkkint = 33 + 7*ndeq
                  lllint = 34
               go to 140
c              (ddeabm)
  130          continue
                  kkkint = 130 + 21*ndeq
                  lllint = 51
  140          continue
c
c              (coef)
               kkkcof = 5*nfcc + nfcc**2
               lllcof = 3 + nfcc
c
               kkkws = max(kkksud,kkksvc,kkkint,kkkcof)
               llliws = max(lllsud,lllsvc,lllint,lllcof)
c
               needw = kkkyhp + kkku + kkkv + kkkcoe + kkks + kkksto
     1                 + kkkg + kkkzpw + kkkws
               neediw = 17 + lllip + llliws
c              *********************************************************
c                  compute the number of possible orthonormalizations
c                  with the allotted storage
c
               iwork(3) = needw
               iwork(4) = kkkzpw
               iwork(5) = neediw
               iwork(6) = lllip
               nrtemp = ndw - needw
               nitemp = ndiw - neediw
c           ...exit
               if (nrtemp .lt. 0) go to 180
c           ...exit
               if (nitemp .lt. 0) go to 180
c
               if (ndisk .eq. 0) go to 150
                  non = 0
                  mxnon = nrtemp
               go to 160
  150          continue
c
                  mxnonr = nrtemp/kkkzpw
                  mxnoni = nitemp/lllip
                  mxnon = min(mxnonr,mxnoni)
                  non = mxnon
  160          continue
c
               iwork(2) = mxnon
c
c              *********************************************************
c                  check for pre-assigned orthonormalization points
c
               nopg = 0
c        ......exit
               if (iwork(11) .ne. 1) go to 210
               if (mxnon .lt. iwork(1)) go to 170
                  nopg = 1
                  mxnon = iwork(1)
                  work(mxnon+1) = 2.0d0*xpts(nxpts) - xpts(1)
c        .........exit
                  go to 210
  170          continue
  180       continue
c
            iflag = -1
      if (ndisk .ne. 1) then
         write (xern1, '(i8)') needw
         write (xern2, '(i8)') kkkzpw
         write (xern3, '(i8)') neediw
         write (xern4, '(i8)') lllip
         call xermsg ('slatec', 'dbvsup',
     *      'required storage for work array is '  // xern1 // ' + ' //
     *      xern2 // '*(expected number of orthonormalizations) $$'  //
     *      'required storage for iwork array is ' // xern3 // ' + ' //
     *      xern4 // '*(expected number of orthonormalizations)', 1, 0)
      else
         write (xern1, '(i8)') needw
         write (xern2, '(i8)') neediw
         call xermsg ('slatec', 'dbvsup',
     *      'required storage for work array is '  // xern1 //
     *      ' + number of orthonomalizations. $$'  //
     *      'required storage for iwork array is ' // xern2, 1, 0)
      endif
      return
c
c        ***************************************************************
c            allocate storage from work and iwork arrays
c
c         (z)
  210    k1 = 1 + (mxnon + 1)
c        (p)
         k2 = k1 + ntp*(non + 1)
c        (w)
         k3 = k2 + nfcc*(non + 1)
c        (yhp)
         k4 = k3 + kkkyhp
c        (u)
         k5 = k4 + kkku
c        (v)
         k6 = k5 + kkkv
c        (coef)
         k7 = k6 + kkkcoe
c        (s)
         k8 = k7 + kkks
c        (stowa)
         k9 = k8 + kkksto
c        (g)
         k10 = k9 + kkkg
         k11 = k10 + kkkws
c                  required additional double precision work space
c                  starts at work(k10) and extends to work(k11-1)
c
c           first 17 locations of iwork are used for optional
c           input and output items
c        (ip)
         l1 = 18 + nfcc*(non + 1)
         l2 = l1 + llliws
c                   required integer work space starts at iwork(l1)
c                   and extends to iwork(l2-1)
c
c        ***************************************************************
c            set indicator for normalization of particular solution
c
         nps = 0
         if (iwork(10) .eq. 1) nps = 1
c
c        ***************************************************************
c            set pivoting parameter
c
         indpvt = 0
         if (iwork(15) .eq. 1) indpvt = 1
c
c        ***************************************************************
c            set other common block parameters
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
         if (iwork(13) .eq. -1) mnswot = max(1,iwork(14))
         xbeg = xpts(1)
         xend = xpts(nxpts)
         xsav = xend
         icoco = 1
         if (inhomo .eq. 3 .and. nopg .eq. 1) work(mxnon+1) = xend
c
c        ***************************************************************
c
         call dexbvp(y,nrowy,xpts,a,nrowa,alpha,b,nrowb,beta,iflag,work,
     1               iwork)
         nfc = nfcc
         iwork(17) = iwork(l1)
  220 continue
      return
      end

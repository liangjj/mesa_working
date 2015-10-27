*deck blktri
      subroutine blktri (iflg, np, n, an, bn, cn, mp, m, am, bm, cm,
     +   idimy, y, ierror, w)
c***begin prologue  blktri
c***purpose  solve a block tridiagonal system of linear equations
c            (usually resulting from the discretization of separable
c            two-dimensional elliptic equations).
c***library   slatec (fishpack)
c***category  i2b4b
c***type      single precision (blktri-s, cblktr-c)
c***keywords  elliptic pde, fishpack, tridiagonal linear system
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine blktri solves a system of linear equations of the form
c
c          an(j)*x(i,j-1) + am(i)*x(i-1,j) + (bn(j)+bm(i))*x(i,j)
c
c          + cn(j)*x(i,j+1) + cm(i)*x(i+1,j) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     i+1 and i-1 are evaluated modulo m and j+1 and j-1 modulo n, i.e.,
c
c          x(i,0) = x(i,n),  x(i,n+1) = x(i,1),
c          x(0,j) = x(m,j),  x(m+1,j) = x(1,j).
c
c     these equations usually result from the discretization of
c     separable elliptic equations.  boundary conditions may be
c     dirichlet, neumann, or periodic.
c
c
c     * * * * * * * * * *     on input     * * * * * * * * * *
c
c     iflg
c       = 0  initialization only.  certain quantities that depend on np,
c            n, an, bn, and cn are computed and stored in the work
c            array  w.
c       = 1  the quantities that were computed in the initialization are
c            used to obtain the solution x(i,j).
c
c       note   a call with iflg=0 takes approximately one half the time
c              as a call with iflg = 1  .  however, the
c              initialization does not have to be repeated unless np, n,
c              an, bn, or cn change.
c
c     np
c       = 0  if an(1) and cn(n) are not zero, which corresponds to
c            periodic boundary conditions.
c       = 1  if an(1) and cn(n) are zero.
c
c     n
c       the number of unknowns in the j-direction. n must be greater
c       than 4. the operation count is proportional to mnlog2(n), hence
c       n should be selected less than or equal to m.
c
c     an,bn,cn
c       one-dimensional arrays of length n that specify the coefficients
c       in the linear equations given above.
c
c     mp
c       = 0  if am(1) and cm(m) are not zero, which corresponds to
c            periodic boundary conditions.
c       = 1  if am(1) = cm(m) = 0  .
c
c     m
c       the number of unknowns in the i-direction. m must be greater
c       than 4.
c
c     am,bm,cm
c       one-dimensional arrays of length m that specify the coefficients
c       in the linear equations given above.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling blktri.  this parameter is
c       used to specify the variable dimension of y.  idimy must be at
c       least m.
c
c     y
c       a two-dimensional array that specifies the values of the right
c       side of the linear system of equations given above.  y must be
c       dimensioned at least m*n.
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.
c             if np=1 define k=int(log2(n))+1 and set l=2**(k+1) then
c                     w must have dimension (k-2)*l+k+5+max(2n,6m)
c
c             if np=0 define k=int(log2(n-1))+1 and set l=2**(k+1) then
c                     w must have dimension (k-2)*l+k+5+2n+max(2n,6m)
c
c       **important** for purposes of checking, the required dimension
c                     of w is computed by blktri and stored in w(1)
c                     in floating point format.
c
c     * * * * * * * * * *     on output     * * * * * * * * * *
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m is less than 5.
c       = 2  n is less than 5.
c       = 3  idimy is less than m.
c       = 4  blktri failed while computing results that depend on the
c            coefficient arrays an, bn, cn.  check these arrays.
c       = 5  an(j)*cn(j-1) is less than 0 for some j. possible reasons
c            for this condition are
c            1. the arrays an and cn are not correct.
c            2. too large a grid spacing was used in the discretization
c               of the elliptic equation.
c            3. the linear equations resulted from a partial
c               differential equation which was not elliptic.
c
c     w
c       contains intermediate values that must not be destroyed if
c       blktri will be called again with iflg=1.  w(1) contains the
c       number of locations required by w in floating point format.
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   an(n),bn(n),cn(n),am(m),bm(m),cm(m),y(idimy,n)
c     arguments      w(see argument list)
c
c     latest         june 1979
c     revision
c
c     required       blktri,blktri,prod,prodp,cprod,cprodp,compb,indxa,
c     subprograms    indxb,indxc,ppadd,psgf,ppsgf,ppspf,bsrh,tevls,
c                    r1mach
c
c     special        the algorithm may fail if abs(bm(i)+bn(j)) is less
c     conditions     than abs(am(i))+abs(an(j))+abs(cm(i))+abs(cn(j))
c                    for some i and j. the algorithm will also fail if
c                    an(j)*cn(j-1) is less than zero for some j.
c                    see the description of the output parameter ierror.
c
c     common         cblkt
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     paul swarztrauber
c
c     language       fortran
c
c     history        version 1 september 1973
c                    version 2 april     1976
c                    version 3 june      1979
c
c     algorithm      generalized cyclic reduction (see reference below)
c
c     space
c     required       control data 7600
c
c     portability    american national standards institute fortran.
c                    the machine accuracy is set using function r1mach.
c
c     required       none
c     resident
c     routines
c
c     references     swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c                    swarztrauber p. ,'a direct method for the discrete
c                    solution of separable elliptic equations', s.i.a.m.
c                    j. numer. anal.,11(1974) pp. 1136-1150.
c
c***references  p. n. swarztrauber and r. sweet, efficient fortran
c                 subprograms for the solution of elliptic equations,
c                 ncar tn/ia-109, july 1975, 138 pp.
c               p. n. swarztrauber, a direct method for the discrete
c                 solution of separable elliptic equations, siam journal
c                 on numerical analysis 11, (1974), pp. 1136-1150.
c***routines called  blktr1, compb, cprod, cprodp, prod, prodp
c***common blocks    cblkt
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  blktri
c
      dimension       an(*)      ,bn(*)      ,cn(*)      ,am(*)      ,
     1                bm(*)      ,cm(*)      ,y(idimy,*) ,w(*)
      external        prod       ,prodp      ,cprod      ,cprodp
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  blktri
      nm = n
      ierror = 0
      if (m-5) 101,102,102
  101 ierror = 1
      go to 119
  102 if (nm-3) 103,104,104
  103 ierror = 2
      go to 119
  104 if (idimy-m) 105,106,106
  105 ierror = 3
      go to 119
  106 nh = n
      npp = np
      if (npp) 107,108,107
  107 nh = nh+1
  108 ik = 2
      k = 1
  109 ik = ik+ik
      k = k+1
      if (nh-ik) 110,110,109
  110 nl = ik
      ik = ik+ik
      nl = nl-1
      iwah = (k-2)*ik+k+6
      if (npp) 111,112,111
c
c     divide w into working sub arrays
c
  111 iw1 = iwah
      iwbh = iw1+nm
      w(1) = iw1-1+max(2*nm,6*m)
      go to 113
  112 iwbh = iwah+nm+nm
      iw1 = iwbh
      w(1) = iw1-1+max(2*nm,6*m)
      nm = nm-1
c
c subroutine comp b computes the roots of the b polynomials
c
  113 if (ierror) 119,114,119
  114 iw2 = iw1+m
      iw3 = iw2+m
      iwd = iw3+m
      iww = iwd+m
      iwu = iww+m
      if (iflg) 116,115,116
  115 call compb (nl,ierror,an,bn,cn,w(2),w(iwah),w(iwbh))
      go to 119
  116 if (mp) 117,118,117
c
c subroutine blktr1 solves the linear system
c
  117 call blktr1 (nl,an,bn,cn,m,am,bm,cm,idimy,y,w(2),w(iw1),w(iw2),
     1             w(iw3),w(iwd),w(iww),w(iwu),prod,cprod)
      go to 119
  118 call blktr1 (nl,an,bn,cn,m,am,bm,cm,idimy,y,w(2),w(iw1),w(iw2),
     1             w(iw3),w(iwd),w(iww),w(iwu),prodp,cprodp)
  119 continue
      return
      end

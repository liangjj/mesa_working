*deck @(#)simeig.f	5.1  11/6/94
      subroutine simeig(r,v,thresh,almax,u,p,q,y,row)
c.01/25/77
c  iterative calculation of several eigenvalues and eigenvectors of
c  a symmetric matrix by the method of optimized relaxations (mor),
c  with successive over-relaxation, applied simultaneously to the
c  several trial vectors, with each iteration followed by
c  diagonalization of the interaction matrix in order to improve the
c  separation of the eigenvector components in the trial vectors.
c
c  this method is suitable primarily for sparse, diagonally dominant
c  matrices.  convergence is slow if there are large off-diagonal
c  elements, but closely spaced roots should not cause serious
c  difficulties, provided all roots of a nearly degenerate set are
c  included in the calculation if any one of them is sought.
c
c  for references on mor see -
c     i. shavitt et al., j. comput. phys. 11, 90 (1973),
c     s. falk, z. angew. math. mech. 53, 73 (1973).
c  for references on over-relaxation see -
c     r. m. nisbet, j. comput. phys. 10, 614 (1972),
c     a. ruhe, math. comp. 28, 695 (1974),
c     h. r. schwarz, comp. meth. appl. mech. eng. 3, 11 (1974).
c  for references on the interaction matrix approach see -
c     m. clint and a. jennings, comput. j. 13, 76 (1970),
c     r. c. raffenetti, j. comput. phys. 32, 403 (1979),
c     i. shavitt, unpublished.
c
c  the matrix is accessed through the subroutine 'sums', leaving
c  'simeig' independent of most details of the matrix storage scheme.
c  in all cases the matrix is assumed to be arranged by rows of the
c  lower triangle, with zero elements omitted (computation time is
c  roughly proportional to the number of elements present), the
c  remaining elements being suitably identified by their column
c  index, and with the diagonal element last in each row.
c
c  the subroutine 'sums' can easily be modified according to the
c  details of the storage scheme used.
c
c  the subroutine 'geneig' is used to solve the generalized
c  eigenvalue problem for the interaction matrix, using schmidt
c  orthogonalization and the jacobi method.
c
c  scalar arguments are in common /args/; array arguments are in the
c  subroutine argument list.  the calling subroutine should provide
c  values for variables n, k, m, l, itmax, omega, for array thresh,
c  for initial-guess trial vectors in array v, and should supply
c  working space for arrays u, p, q, y, row.  information is
c  returned in variables it, noconv, and arrays r, v, almax.  the
c  contents and dimensions of the arguments are -
c
c  n      = order of the matrix (matrix may actually be of higher
c           order, but only the upper-left n*n part will be used).
c  k      = number of roots sought.
c  m,l    = dimensions of the eigenvector array v(l,m).  normally
c           m=n and l=k, but in some cases it may be convenient to
c           solve a smaller problem first (n.lt.m and/or k.lt.l)
c           without having to reorganize the data in array v.
c  itmax  = maximum number of iterations allowed.  a typical value
c           is 50 or 100.
c  it     = number of iterations actually used.
c  omega  = over-relaxation factor, by which each computed increment
c           is multiplied before being added to the trial vector.
c           the suggested value is 1.4, and a reasonable range is
c           1.0-1.8, but in cases of slow convergence a change in
c           omega may help.  (it is suggested that omega be
c           decreased if the q values in the convergence process
c           reports vary erratically or oscillate, and increased
c           if they change very smoothly.)
c  noconv = non-convergence flag.  noconv=1 if convergence has not
c           been achieved in itmax iterations, noconv=0 otherwise.
c  mu     = matrix row currently being processed (used for
c           communication with subroutine 'sums').
c  amumu  = diagonal element in row mu (used for communication with
c           subroutine 'sums').
c  r      = eigenvalues array (length=k)
c  v      = eigenvectors array (length=l*m).  on entry it should
c           contain the initial-guess trial vectors, and on exit it
c           contains the eigenvectors.  storage organization is as a
c           two-dimensional array of dimensions (l,m), where l.ge.k
c           and m.ge.n (see description of m,l above), but it is
c           used here as a one-dimensional array (elements 1,2,...,k
c           contain the first component of the vectors, elements
c           l+1,l+2,...,l+k contain the second component etc., and
c           any unused elements are not disturbed).
c  thresh = convergence thresholds array (length=k).  convergence is
c           assumed if the maximum change in a normalized vector
c           vector component during an iteration is less than the
c           corresponding threshold.  the thresholds for the
c           different vectors are allowed to be different, since not
c           all may be needed to the same accuracy (in particular,
c           some vectors may be included only to help separate
c           nearly degenerate roots).
c  almax  = maximum increments array (length=k).  contains the
c           maximum change in the normalized vector components
c           during the last iteration.
c  u      = working-space array (length=k*n); contains column sums.
c  p      = working-space array (length=k*(k+1)/2); contains v*a*v',
c           where a is the matrix and v' is the transpose of v.  (p
c           is the interaction matrix, with q its associated overlap
c           matrix.)
c  q      = working-space array (length=k*(k+1)/2); contains v*v'.
c  y      = working-space array (length=k**2); contains eigenvectors
c           of the interaction matrix.
c  row    = working-space array (length should be adequate to hold the
c           maximum number of matrix elements present in any row of
c           the lower triangle of the matrix); used by 'sums' to save
c           one row of the matrix between 'sumsb' and 'sumsc' calls.
c
c  in the subroutine, the variables i, j index different vectors,
c  while mu, nu index matrix elements and vector components.
c
      implicit real*8 (a-h,o-z)
c
      common /args/ omega,amumu,bi,n,k,m,l,itmax,it,noconv,mu,idgwrt
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      dimension u(*),v(*),r(*),row(*),almax(*),p(*),q(*),y(*),thresh(*)
c
      parameter (pone=0.98d0)
c  if the square of one component of a normalized trial vector is
c  greater than pone (i.e., the vector is almost parallel to a
c  coordinate axis), this component is not incremented.
c
      parameter (plineq=0.1d0)
c  this determines whether linear or quadratic calculation of the
c  increment alpha is to be used.  the equation
c            a*alpha**2 + b*alpha + c = 0
c  is solved if (4*a*c/b**2).gt.plineq, otherwise
c            alpha = -c/b
c  is used.
c
      parameter (tolmin=0.001d0)
c  this is the square root of the maximum convergence threshold
c  allowed for the jacobi diagonalization of the interaction matrix.
c  the normal threshold is the square of the smallest max(almax(i),
c  thresh(i)) of the previous iteration (but min(thresh(i))**2 is used
c  during the initialization), unless this is bigger than tolmin**2
c
      parameter (ratio=1.3d0)
c  this is used in the calculation of the root shifts.  for the
c  calculation of the i-th root, the j-th root (j.lt.i) is shifted
c  by   ratio*(r(i)-r(1)) - (r(j)-r(1)) .
c
c  initialization section
c
      ratiom=ratio-1.0d0
      tol=tolmin
      do 12 i=1,k
 12     tol=min(tol,thresh(i))
      tol=tol**2
c     (initial tolerance for interaction matrix diagonalization)
c
c  print heading for convergence reports
c
      if (idgwrt.gt.0) write (iw,1) n,k,omega,(i,i=1,k)
 1      format (29h1convergence progress reports//
     1  15h matrix order =,i7/18h number of roots =,i4/
     2  25h over-relaxation factor =,f7.3//
     3  37h t = convergence thresholds specified/
     4  48h i = rayleigh quotients of initial guess vectors/
     5  49h d = rayleigh quotients after k*k diagonalization/
     6  43h m = rayleigh quotients after mor iteration/
     7  50h q = squared length of vectors after mor iteration/
     8  56h v = maximum relative vector increments in mor iteration//
     9  10h iteration,i11,4i22/(i21,4i22))
      if (idgwrt.gt.0) write (iw,11) (thresh(i), i=1,k)
 11     format (/2h t,5x,1p,5d22.4/(7x,5d22.4))
c
c  compute column sums u(i,nu)=sum(mu=nu+1,n)a(mu,nu)*v(i,mu)
c  and interaction matrices q(i,j)=sum(mu)v(i,mu)*v(j,mu),
c  p(i,j)=sum(mu,nu)v(i,mu)*a(mu,nu)*v(j,nu)
c
      ij=0
      do 2 i=1,k
        do 21 j=1,i
          ij=ij+1
          p(ij)=0.0d0
 21       q(ij)=0.0d0
 2      continue
      muz=0
      muy=0
c     (during the mu loops, muz=(mu-1)*l and muy=(mu-1)*k.  thus
c     muz+i and muy+i are the indices for accessing v(i,mu) and
c     u(i,mu), respectively)
      rewind iunt2a
      read (iunt2a)
      do 3 mu=1,n
c
        call sumsa(u,v(muz+1),row)
c
c       (this call adds a(mu,nu)*v(i,mu) to u(i,nu) for nu=1,mu-1
c       and i=1,k, and places a(mu,mu) in 'amumu')
c
        ij=0
        do 31 i=1,k
          do 311 j=1,i
            ij=ij+1
            vij=v(muz+i)*v(muz+j)
            q(ij)=q(ij)+vij
 311        p(ij)=p(ij)+vij*amumu
c         (only the diagonal contributions to p have been computed here)
 31       u(muy+i)=0.0d0
c       (the mu-th column of u is initialized in preparation for
c       incrementation during subsequent sweeps over the mu loop)
        muy=muy+k
 3      muz=muz+l
      muz=0
      muy=0
      do 4 mu=1,n
        ij=0
        do 41 i=1,k
          do 411 j=1,i
            ij=ij+1
 411        p(ij)=p(ij)+u(muy+i)*v(muz+j)+v(muz+i)*u(muy+j)
 41       r(i)=p(ij)/q(ij)
c       (r contains the rayleigh quotients)
        muy=muy+k
 4      muz=muz+l
c
c  this completes the calculation of p and q
c
c  print initial guess results
c
      if (idgwrt.gt.0) write (iw,5) (r(i), i=1,k)
 5      format (2h i,5x,1p,5d22.12/(7x,5d22.12))
c
c  begin iterations
c
      do 6 it=1,itmax
        write(iw,*)'simeig : it, idgwrt ', it,idgwrt
        if (k.eq.1) then
c         (if only one root is sought, skip diagonalization)
          fnorm=1.0d0/sqrt(q(1))
          p(1)=r(1)
          q(1)=1.0d0
          muz=1
          do 62 mu=1,n
            v(muz)=v(muz)*fnorm
            u(mu)=u(mu)*fnorm
 62         muz=muz+l
        else
c
c         first solve the generalized eigenvalue problem for the k*k
c         interaction matrix
c                              p*y = q*y*r  ,
c         where y'*q*y=1 and r is a diagonal matrix of eigenvalues.
c         the diagonalized matrix r=y'*p*y replaces p, and q is
c         destroyed.  the diagonalization (by the jacobi method) is
c         carried out to the specified tolerance 'tol'.
c
          call geneig(k,p,q,y,tol)
c
c         transform v and u by the y matrix,
c                     v(i,mu)=sum(j)v(j,mu)*y(j,i)  ,
c                     u(i,mu)=sum(j)u(j,mu)*y(j,i)  ,
c         where y(j,i) is the j-th component of the i-th eigenvector
c         of the interaction matrix
c
          muz=0
          muy=0
          do 64 mu=1,n
            iz=0
c           (during the i loop, iz=(i-1)*k.  thus iz+j is the index
c           for accessing y(j,i))
            do 641 i=1,k
              r(i)=0.0d0
              almax(i)=0.0d0
c             (r and almax are used here for temporary storage)
              do 6411 j=1,k
                r(i)=r(i)+v(muz+j)*y(iz+j)
 6411           almax(i)=almax(i)+u(muy+j)*y(iz+j)
 641          iz=iz+k
            do 642 i=1,k
              v(muz+i)=r(i)
 642          u(muy+i)=almax(i)
            muy=muy+k
 64         muz=muz+l
c
c         replace q by a unit matrix and copy the rayleigh quotients r
c
          ij=0
          do 65 i=1,k
            do 651 j=1,i
              ij=ij+1
 651          q(ij)=0.0d0
            q(ij)=1.0d0
 65         r(i)=p(ij)
c
c         print diagonalization report
c
          if (idgwrt.gt.0) write (iw,652) (r(i),i=1,k)
 652        format (/2h d,5x,1p,5d22.12/(7x,5d22.12))
c
        endif
        do 631 i=1,k
 631      almax(i)=0.0d0
c       (almax(i) will contain the maximum alpha(i) of the iteration)
c
c       begin relaxation steps
c
        muz=0
        muy=0
        rewind iunt2a
        read (iunt2a)
        do 66 mu=1,n
c
          call sumsb(v,u(muy+1),row)
c
c         (this call adds sum(nu=1,mu-1)a(mu,nu)*v(i,nu) to u(i,mu)
c         for i=j,k, and places a(mu,mu) in 'amumu'.  at this point
c         u(i,mu) = sum(nu=1,n)a(mu,nu)*v(i,nu) - a(mu,mu)*v(i,mu))
c
          iz=0
          do 661 i=1,k
            ii=iz+i
c           (ii is the index for accessing q(i,i) and p(i,i))
            alpha=0.0d0
            mui=muz+i
            if ((v(mui)**2).gt.q(ii)*pone) go to 661
c           (if the trial vector is essentially parallel to the
c           mu-th coordinate axis, do not increment)
c
c           compute the coefficients of the quadratic equation
c                   a*alpha**2 + b*alpha + c = 0
c           for the vector component increment alpha
c
            a=-u(muy+i)
            b=amumu*q(ii)-p(ii)
            if (i.eq.1) go to 6612
c
c           unless first root, compute shifts in lower roots,
c              shift(j) = (ratio*(r(i)-r(1))-(r(j)-r(1)))/q(j,j)
c                       = (ratio*r(i)-y(j))/q(j,j)
c           where  y(j) = r(j)+(ratio-1)*r(1) .
c
            y(i-1)=ratiom*r(1)+r(i-1)
            shifti=ratio*r(i)
            svv=0.0d0
            jj=0
            do 6613 j=1,i-1
              jj=jj+j
              shiftj=(shifti-y(j))/q(jj)
c
c             adjust a, b for root shifts
c
              sv=shiftj*v(muz+j)
              svv=svv+sv*v(muz+j)
              a=a-sv*q(iz+j)
 6613         b=b-shiftj*q(iz+j)**2
            a=a+svv*v(mui)
            b=b+svv*q(ii)
 6612       c=b*v(mui)-a*q(ii)
            cc=-(c+c)
            ac=(a+a)*cc
c           (ac=-4*a*c)
            bsq=b**2
c
c           use quadratic equation only if 4*a*c is not very
c           small compared to b**2
c
            if (abs(ac).le.bsq*plineq) then
c
c             section for linear equation
c
              if (b.eq.0.0d0) go to 661
c             (if a=b=c=0, no incrementation)
              alpha=-c/b
            else
c
c             section for quadratic equation
c
              alpha=cc/(b+sign(sqrt(bsq+ac),b))
            endif
c
c           joint section
c
            alpha=alpha*omega
c           (omega is the over-relaxation factor)
            almax(i)=max(almax(i),abs(alpha))
c           (almax contains the maximum increments of the iteration)
c
c           update v, q, p
c
            ij=iz
            do 6616 j=1,i
              ij=ij+1
              q(ij)=q(ij)+v(muz+j)*alpha
 6616         p(ij)=p(ij)+(u(muy+j)+amumu*v(muz+j))*alpha
            v(mui)=v(mui)+alpha
            do 6617 j=i,k
              q(ij)=q(ij)+v(muz+j)*alpha
              p(ij)=p(ij)+(u(muy+j)+amumu*v(muz+j))*alpha
 6617         ij=ij+j
            r(i)=p(ii)/q(ii)
 661        iz=ii
c
c         compute column sums
c
          do 662 i=1,k
 662        u(muy+i)=0.0d0
c
          call sumsc(u,v(muz+1),row)
c
c         (this call adds a(mu,nu)*v(i,mu) to u(i,nu) for
c         nu=j,mu-1 and i=1,k)
c
          muy=muy+k
 66       muz=muz+l
c
c       end of relaxation steps
c
c       print progress report
c
        ii=0
        do 67 i=1,k
          ii=ii+i
          y(i)=q(ii)
c         (y is used here for temporary storage)
 67       almax(i)=almax(i)/sqrt(q(ii))

        write(iw,*)'simeig: it,r(1)',it,r(1)

        if (idgwrt.gt.0) write (iw,68) it,(r(i), i=1,k)
 68       format (2h m,i5,1p,5d22.12/(7x,5d22.12))
        if (idgwrt.gt.0) write (iw,69) (y(i), i=1,k)
 69       format (2h q,5x,5f22.10/(7x,5f22.10))
        if (idgwrt.gt.0) write (iw,610) (almax(i), i=1,k)
 610      format (2h v,5x,1p,5d22.3/(7x,5d22.3))
c
c       test convergence and set tolerance for 'geneig'
c
        noconv=0
        tol=tolmin
        do 620 i=1,k
          if (almax(i).gt.thresh(i)) noconv=1
c         (if not yet converged, set flag)
          tol=min(tol,max(almax(i),thresh(i)))
 620      continue
        tol=tol**2
        if (noconv.eq.0) return
c       (return if convergence has been achieved)
 6      continue
c
c  end of iteration loop.  if came through here, convergence has not
c  been achieved in itmax iterations.
c
      it=itmax
      write (iw,630)
 630    format(/'********** no convergence **********')
c
      return
c
      end

*deck @(#)solver.f	5.1  11/6/94
      subroutine solver(oper,status,v1,v2,n,maxvec,i3,a,eigvec)
c
c***begin prologue     solver
c***date written       851007   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           linear equations, large matrices,
c                      linear algorithm, linear algebra,
c                      davidson algorithm
c***author             saxe, paul,    (lanl)
c***source             @(#)solver.f	5.1   11/6/94
c***purpose            iterative solution of an system of linear
c                      equations.
c***description
c
c    solver is intended for solving one or a few sets of linear
c    equations with a large symmetric and diagonally
c    dominant matrix. only the product of the matrix and a vector
c    are needed, not explicitly the matrix. the work required is
c    roughly proportional to the number of solutions desired times
c    the dimension of the matrix squared.
c
c    the arguments depend on the first argument, which controls
c    essentially different entries into the routine. each entry
c    will be described individually.
c
c    currently the entries or modes of operation are:
c       'initialize'
c       'with rhs'
c       'with vector'
c       'solve'
c       'new trial'
c       'get vector'
c       'finish'
c    and they must be called in essentially this order. 'initialize'
c    is used to pass in thresholds and cutoffs, as well as the
c    diagonals of the matrix. after initializing, at least one right
c    hand side must be passed in with 'with rhs'. after the rhs are
c    passed in, set up an iterative loop, asking for 'trial vector'
c    and passing in the product of the trial vector with the matrix
c    (less diagonal) until status is returned as 'done'. then call
c    'solve' to set up and solve the effective linear equations.
c    after calling 'solve', repeat by calling 'trial vector' to get the
c    next trial vector, then form the product of the trial vector
c    with the matrix (less diagonals). enter the new trial and product
c    vectors into the sequence with a call to 'with vector'. when
c    a call to 'trial vector' returns a status of 'done', there are
c    no more orthonormal trial vectors to the cutoff 'thresh', so call
c    'solve' and start another cycle. finally, a call to 'solve' will
c    return a status of 'converged', signifying that all solutions
c    have converged to the requested tolerance.
c
c    a note about the various cutoffs. new trial vectors are normalized
c    and then schmidt orthogonalized to the previous trial vectors.
c    if the norm of the new trial vector falls below 'thresh' during
c    the orthogonalization, the new trial vector is discarded. the
c    total convergence is controlled by 'convgnc'. when none of the
c    norms of the solutions for all of the equations change less
c    than 'cnvgnc', the iterations are considered complete. if all of
c    the trial vectors are thrown out by 'thresh' before the
c    convergence criterion is met, a single trial vector is constructed
c    from the least converged equation. thus 'thresh' will not affect
c    the final convergence, but just the way it is reached.
c
c
c----------------------------------------------------------------------
c 'initialize'
c
c   call solver('initialize',char,diag,thresh,n,maxvec,iout,
c             efixed,cnvgnc)
c
c  on input:
c
c    diag      real (n)
c              the diagonal of the matrix.
c
c    thresh    real
c              the cutoff on the norm of the schmidt orthogonalized
c              potential trial vector. those smaller than "thresh" will
c
c    n         integer
c              the dimension of the matrix and all vectors.
c
c    maxvec    integer
c              the maximum number of expansion vectors permissable.
c              if this number is exceeded, lnkerr is called, aborting
c              the run.
c
c    iout      integer
c              the fortran unit for output about the iterations.
c
c    efixed    real
c              a constant to subtract from the diagonals of the matrix,
c              so one is solving  (a-efixed)x=y
c
c    cnvgnc    real
c              the convergence criterion for the vectors, essentially
c              the change in the norm of the vectors.
c
c  unused arguments:  char (character)
c
c---------------------------------------------------------------------
c  'with rhs'
c
c     call solver('with rhs',char,rhs,r1,n,i1,i2,r2,r3)
c
c
c   on input:
c
c     rhs       real (n)
c               one of the right hand sides for the sets of equations.
c
c     n         integer
c               dimension of the vectors, which is the number of
c               variables in the linear equations.
c
c   unused arguments:
c
c     char      character
c     r1,r2,r3  real arrays
c     i1,i2     integers
c----------------------------------------------------------------------
c
c  'with vectors'
c
c     call solver('with vectors',char,trial,product,n,i1,i2,r1,r2)
c
c  on input:
c
c    trial      real (n)
c               the trial vector, which must be orthonormal to any
c               other trial vectors.
c
c    product    real (n)
c               the product of the matrix, less the diagonals, times
c               the trial vector.
c
c    n          integer
c               the size of the matrix and vectors.
c
c unused arguments:
c
c    i1,i2      integers
c    r1,r2      real arrays
c----------------------------------------------------------------------
c
c  'solve'
c
c    call solver('solve',status,v1,v2,n,maxvec,i1,a,eigvec)
c
c  on input:
c
c   n        integer
c            the size of a side of the matrix, and the length of the
c            vectors.
c
c   maxvec   integer
c            the maximum number of expansion vectors permitted.
c
c
c  on output:
c
c   status   character
c            either 'continue' or 'converged'.
c
c   eigvec   real (maxvec,neqn)
c            the eigenvectors of the effective eigenvalue problem.
c
c  scratch:
c
c   v1,v2    real (n)
c   a        real (maxvec*(maxvec+1)/2)
c
c  unused arguments:
c
c   i1       integer
c----------------------------------------------------------------------
c
c  'new trial'
c
c     call solver('new trial',status,trial,v1,n,maxvec,i1,r1,eigvec)
c
c
c  on entry:
c
c   n        integer
c            the length of the vectors.
c
c   maxvec   integer
c            the maximum number of expansion vectors permissible.
c
c   eigvec   real (maxvec,number of sets of equations)
c
c  on exit:
c
c   status   character
c            either 'trial' if there is another trial vecotr, or 'done'
c            if there are no more trial vectors for this cycle.
c
c   trial    real (n)
c            the next orthonormal trial vector in the expansion set.
c
c  scratch:
c
c   v1       real (n)
c
c  unused arguments:
c
c   i1       integer
c   r1       real array
c----------------------------------------------------------------------
c
c 'get vector'
c
c     call solver('get vector',status,vector,r1,n,ieqn,i1,r2,r3)
c
c  on input:
c
c   n        integer
c            the length of the vectors.
c
c   ieqn     integer
c            the equation for which the vector is desired.
c
c  on output:
c
c   status   character
c            'ok', 'error', or 'not converged'
c
c   vector   real (n)
c            the vector, or a vector of 1.0d+99 if an error occurred.
c
c  unreferenced arguments:
c
c   r1,r2    real arrays
c   r3
c   i1       integer
c---------------------------------------------------------------------
c
c  'finish'
c
c     call solver('finish',c1,r1,r2,i1,i2,i3,r3,r4)
c
c all arguments are unreferenced. c is for character; i, integer; and
c r, real arrays.
c----------------------------------------------------------------------
c
c***references         b. liu, "the simultaneous expansion method for
c          the iterative solution of several of the lowest eigenvalues
c          and corresponding eigenvectors of large real-symmetric
c          matrices" in numerical algorithms in chemistry:
c          algebraic methods, a report of an nrcc workshop, 1978.
c
c          e. r. davidson, j. computational phys. 17, p87 (1975).
c
c          see also the routine bliu for eigen equations.
c
c***routines called    unqfil (mdutil)
c                      iosys  (io)
c                      rzero  (math)
c                      sdot   (calmath)
c                      saxpy  (calmath)
c                      lnkerr (mdutil)
c                      sscal  (calmath)
c                      itoc   (chr)
c                      sspfa  (clams)
c                      sspsl  (clams)
c
c***end prologue       solver
c
c
      implicit integer (a-z)
c
      real*8 v1(n),v2(n),a(*),eigvec(maxvec,*)
      real*8 sdot,sqcdif,maxdif,nrzero,norm,thresh,efixed,t
      real*8 ovrlap,cnvgnc
      character*(*) oper,status
      character*8 cjunk,itoc*4
c
      data nrzero / 1.0d-06 /
c
      save added,teqn,neqn,nrzero,thresh,cycle,efixed,iout
      save cnvgnc,oldnum,maxdif,numdif
c
c     ----- decide what to do based on "oper" passed in -----
c
c-------------------------------------------------------- 'initialize'
c
      if (oper.eq.'initialize') then
c
c        ----- save local copies of important variables -----
c
         iout=i3
         thresh=v2(1)
         efixed=a(1)
         cnvgnc=eigvec(1,1)
c
c        ----- create a scratch file for storage of intermediates -----
c
         junk=min(262144,(2*maxvec+1)*((n+511)/512*512))
         call iosys('open solver as scratch on ssd',junk,0,0,cjunk)
c
c        ----- store the diagonals on the scratch file -----
c
         call iosys('write real diagonals on solver',n,v1,0,' ')
c
c        ----- form the inverse of the diagonals, and store -----
c
         count=0
         do 1 i=1,n
            t=v1(i)-efixed
            if (t.eq.0.0d+00) then
               count=count+1
               v1(i)=1.0d+00
           else
               v1(i)=1.0d+00/t
            end if
    1    continue
c
         if (count.gt.0) write (iout,2) count
    2    format (' encountered ',i4,' singularities in the inverse '//
     #           'of the diagonal of the hessian')
c
         call iosys('write real inverse on solver',n,v1,0,' ')
c
c        ----- and initialize some local variables -----
c
         oldnum=0
         added=0
         cycle=0
         neqn=0
c
c----------------------------------------------------- 'with rhs' -----
c
      else if (oper.eq.'with rhs') then
         neqn=neqn+1
         call iosys('write real "rhs:'//itoc(neqn)//'" to solver',
     #               n,v1,0,' ')
c
c----------------------------------------------------- 'with vectors'
c
      else if (oper.eq.'with vectors') then
c
c        ----- check for exceeding maximum number of expansion vectors
c
         if (oldnum+added.ge.maxvec) then
            call lnkerr('solver exceeded maximum number, '//
     #                  itoc(maxvec)//' of expansion vectors')
         end if
c
c        ----- save this trial and product vector to disk -----
c
         added=added+1
         call iosys('write real "trial:'//itoc(oldnum+added)//
     #              '" to solver', n,v1,0,' ')
         call iosys('write real "product:'//itoc(oldnum+added)//
     #              '" to solver',n,v2,0,' ')
c
c        ----- add in product of trial times diagonal of hessian -----
c
         call iosys('read real diagonals from solver',-1,v2,0,' ')
         do 3 i=1,n
            v1(i)=v1(i)*(v2(i)-efixed)
    3    continue
c
         call iosys('read real "product:'//itoc(oldnum+added)//
     #              '" from solver',-1,v2,0,' ')
c
         do 4 i=1,n
            v2(i)=v2(i)+v1(i)
    4    continue
c
         call iosys('write real "product:'//itoc(oldnum+added)//
     #              '" to solver',n,v2,0,' ')
c
c------------------------------------------------------- 'solve' ------
c
      else if (oper.eq.'solve') then
         status='continue'
c
c        ----- set up and solve the small set of equations -----
c
         if (oldnum.gt.0) then
            call iosys('read real "small matrix" from solver',-1,a,0,
     $           ' ')
            call iosys('read real rhs from solver',-1,eigvec,0,' ')
         else
            call rzero(a,maxvec*(maxvec+1)/2)
            call rzero(eigvec,maxvec*neqn)
         end if
c
c        ----- first read in new product vector -----
c
         do 10 vector=oldnum+1,oldnum+added
c
            call iosys('read real "trial:'//itoc(vector)//
     #                 '" from solver',-1,v1,0,' ')
c
c           ----- add to the right hand side <x/y> -----
c
            do 8 i=1,neqn
               call iosys('read real "rhs:'//itoc(i)//'" from solver',
     #                     -1,v2,0,' ')
               eigvec(vector,i)=sdot(n,v1,1,v2,1)
    8       continue
c
c           ----- add the next row to the small matrix -----
c
            call iosys('read real "product:'//itoc(vector)//
     #                 '" from solver',-1,v1,0,' ')
            ij=vector*(vector-1)/2
            do 9 trial=1,vector
               call iosys('read real "trial:'//itoc(trial)//
     #                    '" from solver',-1,v2,0,' ')
               ij=ij+1
               a(ij)=sdot(n,v1,1,v2,1)
    9       continue
   10    continue
c
c        ----- save the small matrix -----
c
         oldnum=oldnum+added
         call iosys('write real "small matrix" to solver',
     #              maxvec*(maxvec+1)/2,a,0,' ')
         call iosys('write real rhs to solver',maxvec*neqn,eigvec,0,' ')
c
c        ----- solve the small set of equations -----
c
         call sspfa(a,oldnum,v1,info)
         if (info.ne.0) call lnkerr('error factoring matrix:'//
     #          itoc(info))
         do 7 i=1,neqn
            call sspsl(a,oldnum,v1,eigvec(1,i))
    7    continue
c
c        ----- print out results for this cycle -----
c
         cycle=cycle+1
         write (iout,11) cycle,oldnum
   11    format (/,1x,25('-'),' cycle ',i2,' using ',i4,' vectors ',
     #           10('-'),/,t5,'equation',t30,'convergence')
         maxdif=0.0d+00
         numdif=0
         do 15 eqn=1,neqn
            junk=oldnum-added+1
            sqcdif=sqrt(sdot(added,eigvec(junk,eqn),1,
     #                        eigvec(junk,eqn),1)/
     #             sdot(oldnum,eigvec(1,eqn),1,eigvec(1,eqn),1))
            write (iout,12) eqn,sqcdif
   12       format (t4,i4,t30,f15.9)
            if (sqcdif.gt.maxdif) then
               maxdif=sqcdif
               numdif=eqn
            end if
   15    continue
c
c        ----- form the new "best" vectors -----
c
         do 20 eqn=1,min(neqn,oldnum)
            call rzero(v1,n)
            do 16 i=1,oldnum
               call iosys('read real "trial:'//itoc(i)//
     #                    '" from solver',-1,v2,0,' ')
               call saxpy(n,eigvec(i,eqn),v2,1,v1,1)
   16       continue
            call iosys('write real "vector:'//itoc(eqn)//'" to solver',
     #                 n,v1,0,' ')
   20    continue
c
c        ----- check if all equations are converged -----
c
         if (maxdif.lt.cnvgnc) then
            status='converged'
            return
         end if
c
         added=0
         teqn=0
c
c---------------------------------------------------------- 'new trial'
c
      else if (oper.eq.'new trial') then
c
c        ----- create a new trial vector, start by forming (y-hv)/h(ii)
c
         status=' '
   21    continue
         teqn=teqn+1
c
c        ----- if first cycle, then trial vectors simply formed
c              from rhs
c
         if (oldnum.eq.0) then
            if (teqn.gt.neqn) then
               status='done'
               return
            end if
c
            call iosys('read real "rhs:'//itoc(teqn)//'" from solver',
     #                  -1,v1,0,' ')
         else
c
c           ----- check if we are done with new trial vectors -----
c
            if (teqn.gt.min(neqn,oldnum)) then
               status='done'
               if (added.gt.0) then
               return
               else
                  teqn=numdif
               end if
            end if
c
            call iosys('read real "rhs:'//itoc(teqn)//'" from solver',
     #                  -1,v1,0,' ')
c
            do 25 i=1,oldnum
               call iosys('read real "product:'//itoc(i)//
     #                    '" from solver',-1,v2,0,' ')
               call saxpy(n,-eigvec(i,teqn),v2,1,v1,1)
   25       continue
         end if
c
         call iosys('read real inverse from solver',-1,v2,0,' ')
         do 24 i=1,n
            v1(i)=v1(i)*v2(i)
   24    continue
c
c        ----- and orthonormalize to previous vectors -----
c
         norm=sdot(n,v1,1,v1,1)
         call sscal(n,1.0d+00/sqrt(norm),v1,1)
         norm=1.0d+00
         do 35 i=1,oldnum+added
            call iosys('read real "trial:'//itoc(i)//'" from solver',
     #                 -1,v2,0,' ')
            ovrlap=sdot(n,v1,1,v2,1)
            norm=norm-ovrlap**2
            if (sqrt(norm).lt.thresh.and.status.ne.'done') go to 21
            call saxpy(n,-ovrlap,v2,1,v1,1)
   35    continue
c
         norm=sqrt(sdot(n,v1,1,v1,1))
         if (norm.lt.thresh.and.status.ne.'done') go to 21
         if (status.eq.'done') teqn=neqn+1
         call sscal(n,1.0d+00/norm,v1,1)
c
c        ----- return this vector with an 'all okay' -----
c
         status='trial'
c
c--------------------------------------------------------- 'get vector'
c
      else if (oper.eq.'get vector') then
         if (i3.le.0.or.i3.gt.neqn) then
            status='error'
            do 40 i=1,n
               v1(i)=1.0d+99
   40       continue
            return
         end if
c
         if (maxdif.ge.cnvgnc) then
            status='not converged'
         else
            status='ok'
         end if
         call iosys('read real "vector:'//itoc(i3)//'" from solver',
     #              -1,v1,0,' ')
c
c----------------------------------------------------------- 'finish'
c
      else if (oper.eq.'finish') then
         call iosys('destroy solver',0,0,0,' ')
c
c------------------------------------------------------------- errors
c
      else
         call lnkerr('bad call to solver')
      end if
c
c
      return
      end

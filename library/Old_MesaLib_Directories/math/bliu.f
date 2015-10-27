*deck @(#)bliu.f	5.1  11/6/94
      subroutine bliu(oper,status,v1,v2,n,maxvec,i3,i4,i5,
     $                a,root,eigvec)
c
c***begin prologue     bliu
c***date written       850925   (yymmdd)
c***revision date      910507   (yymmdd)
c
c     7 may    1991 bhl at llnl
c        adding a second orthogonalization for new trial vectors.
c     7 march  1991 rlm at lanl
c        fixing a mixed mode problem with 'initialize'.
c        added the integer argument i5, also adding variables
c        to the list which must be saved from one call to another.
c     9 august 1987 pws at lanl
c        adding the ability to do just a few roots at a time.
c        changing the call to givens to one for rsp.
c
c***keywords           eigenvalue, large matrices, davidson algorithm,
c                      liu algorithm, diagonalization
c***author             saxe, paul,    (lanl)
c***purpose            iterative solution of an eigenvalue problem for
c                      one or a few of the lowest eigenvalues and
c                      vectors.
c***description
c
c    bliu is intended for extracting one or a few of the lowest
c    eigenvalues and vectors of a large symmetric and diagonally
c    dominant matrix. only the product of the matrix and a vector
c    are needed, not explicitly the matrix. the work required is
c    roughly proportional to the number of eigenvalues desired times
c    the dimension of the matrix squared.
c
c    the arguments depend on the first argument, which controls
c    essentially different entries into the routine. each entry
c    will be described individually.
c
c    currently the four entries or modes of operation are:
c       'initialize'
c       'with vector'
c       'solve'
c       'new trial'
c       'get vector'
c       'finish'
c    and they must be called in essentially this order. 'initialize'
c    is used to pass in thresholds and cutoffs, as well as the
c    diagonals of the matrix. after initializing, at least one guess
c    vector and its product with the matrix (less the diagonal) must
c    be passed in using 'with vector'. there is no restriction on the
c    number of guess vectors that may be used--there may be more or
c    less than nroots. however, the user must ensure that the guess
c    vectors are orthonormal. after the guess vectors are read in,
c    then call 'solve' to set up and solve the effective eigenvalue
c    equation. after calling 'solve', call 'trial vector' to get the
c    next trial vector, then form the product of the trial vector
c    with the matrix (less diagonals). enter the new trial and product
c    vectors into the sequence with a call to 'with vector'. when
c    a call to 'trial vector' returns a status of 'done', there are
c    no more orthonormal trial vectors to the cutoff 'thresh', so call
c    'solve' and start another cycle.
c
c    a note about the various cutoffs. new trial vectors are normalized
c    and then schmidt orthogonalized to the previous trial vectors.
c    if the norm of the new trial vector falls below 'thresh' during
c    the orthogonalization, the new trial vector is discarded. the
c    total convergence is controlled by 'convgnc'. when none of the
c    norms of the eigenvectors for all of the roots change less than
c    'cnvgnc', the iterations are considered complete. if all of the
c    trial vectors are thrown out by 'thresh' before the convergence
c    criterion is met, a single trial vector is constructed from the
c    least converged root. thus 'thresh' will not affect the final
c    convergence, but just the way it is reached.
c
c    also, if there are not as many iterates as roots requested, e.g.
c    if five roots are requested but only two guess vectors provided,
c    the routine will only give the eigenvalues and trial vectors for
c    as many roots as possible. as the iterations continue and more
c    vectors are adeed to the set, the other roots will be included.
c
c----------------------------------------------------------------------
c 'initialize'
c
c   call bliu('initialize',print,diag,thresh,n,maxvec,nroots,iout,
c             nattim,efixed,cnvgnc,junk)
c
c  on input:
c
c    print     character
c              'noprint' for no printing.
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
c    nroots    integer
c              the number of eigenvalues to find.
c
c    iout      integer
c              the fortran unit for output about the iterations.
c
c    nattim
c              the number of roots to determine at a time.
c
c    efixed    real
c              a constant to add to the eigenvalues, e.g. nuclear
c              repulsion and frozen-core energies.
c
c    cnvgnc    real
c              the convergence criterion for the vectors, essentially
c              the change in the norm of the vectors.
c
c    junk      real array
c              unused
c----------------------------------------------------------------------
c
c  'with vectors'
c
c     call bliu('with vectors',char,trial,product,n,i1,i2,i3,i4,
c               r1,r2,r3)
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
c    i1,i2,i3,i4   integers
c    r1,r2,r3   real arrays
c----------------------------------------------------------------------
c
c  'solve'
c
c    call bliu('solve',status,v1,v2,n,maxvec,i1,i2,i3,a,root,eigvec)
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
c   root     real (nroots)
c            the unaltered eigenvalue approximations.
c
c   eigvec   real (maxvec,nroots)
c            the eigenvectors of the effective eigenvalue problem.
c
c  scratch:
c
c   v1,v2    real (n)
c   a        real (maxvec*(maxvec+1)/2)
c
c  unused arguments:
c
c   i1,i2,i3 integer
c----------------------------------------------------------------------
c
c  'new trial'
c
c     call bliu('new trial',status,trial,v1,n,maxvec,i1,i2,i3,r1,
c                r1,root,eigvec)
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
c   root     real (nroots)
c            the eigenvalues of the effective eigenvalue problem.
c            this is output when running in 'solve' mode.
c
c   eigvec   real (maxvec,nroots)
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
c   i1,i2,i3 integer
c   r1       real array
c----------------------------------------------------------------------
c
c 'get vector'
c
c     call bliu('get vector',status,vector,r1,n,maxvec,iroot,i1,i2,
c                r2,r3,r4)
c
c  on input:
c
c   n        integer
c            the length of the vectors.
c
c   maxvec   integer
c            the maximum number of expansion vectors permissible.
c
c   iroot    integer
c            the root for which the vector is desired.
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
c   i1,i2    integer
c
c   r1,r2    real arrays
c   r3,r4
c---------------------------------------------------------------------
c
c  'finish'
c
c     call bliu('finish',c1,r1,r2,i1,i2,i3,i4,i5,r3,r4,r5)
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
c***routines called    unqfil (mdutil)
c                      iosys  (io)
c                      rzero  (math)
c                      sdot   (calmath)
c                      saxpy  (calmath)
c                      lnkerr (mdutil)
c                      sscal  (calmath)
c                      itoc   (chr)
c                      givens (math)
c
c***end prologue       bliu
c
c
      implicit integer (a-z)
c
      parameter (mxatim=50)
c
      real*8 v1(n),v2(n),a(*),root(*),eigvec(maxvec,*)
      real*8 sdot,sqcdif,maxdif,eroot,nrzero,norm,thresh,efixed
      real*8 ovrlap,cnvgnc
      character*(*) oper,status
      character*4 itoc
      logical prnt,diag
      logical debug
      integer ptroot(mxatim)
c
      common /io/ inp,itape6
c
      parameter (debug=.false.)
      data nrzero / 1.0d-06 /
c
      save oldnum,added,totnum,nguess
      save diag,cycle,mxroot,maxdif,cnvgnc,efixed,nroots
      save troot,ntodo,newtr,numdif,thresh
      save nrzero,ptroot,nattim,prnt,iout
c
c
c     ----- decide what to do based on "oper" passed in -----
c
c-------------------------------------------------------- 'initialize'
      if (oper(1:10).eq.'initialize') then
c
c        ----- save local copies of important variables -----
c
         nattim=i5
         nattim=min(max(1,nattim),mxatim)
         mnroot=1
         mxroot=mnroot+nattim-1
         nroots=i3
         iout=i4
         thresh=v2(1)
         efixed=a(1)
         cnvgnc=root(1)
         prnt=status.ne.'noprint'
         diag=oper(12:).eq.'with diagonals'
c
c        ----- create a scratch file for storage of intermediates -----
c
         junk=max(262144,(2*maxvec+1)*((n+511)/512*512)/20)
         call iosys('open bliu as scratch on ssd',junk,0,0,' ')
c
c        ----- store the diagonals on the scratch file -----
c
         call iosys('write real diagonals on bliu',n,v1,0,' ')
c
c        ----- and initialize some local variables -----
c
         do 401 i=1,nattim
            ptroot(i)=i
 401     continue
c
         oldnum=0
         newtr=0
         added=0
         cycle=0
         totnum=0
         nguess=0
c----------------------------------------------------- 'with vectors'
      else if (oper.eq.'with vectors') then
c
c        ----- check for exceeding maximum number of expansion vectors
c
         if (oldnum+added.ge.maxvec) then
            call lnkerr('bliu exceeded maximum number, '//itoc(maxvec)
     #                  //' of expansion vectors')
         end if
c
c        ----- save this trial and product vector to disk -----
c
         added=added+1
         totnum=totnum+1
         if (oldnum.eq.0) then
            nguess=nguess+1
            call iosys('write real "trial:'//itoc(oldnum+added)//
     #                 '" to bliu',n,v1,0,' ')
            call iosys('write real "guess:'//itoc(nguess)//'" to bliu',
     $                  n,v1,0,' ')
            call iosys('write real "h*guess:'//itoc(nguess)//
     $                  '" to bliu',n,v2,0,' ')
         end if
         call iosys('write real "product:'//itoc(oldnum+added)//
     #              '" to bliu',n,v2,0,' ')
c------------------------------------------------------- 'solve'
      else if (oper.eq.'solve') then
         status='continue'
c
c        ----- set up and solve the small set of equations -----
c
         if (oldnum.gt.0) then
            call iosys('read real "small matrix" from bliu',-1,a,0,' ')
         else
            call rzero(a,maxvec*(maxvec+1)/2)
         end if
c
c        ----- first, add product of trial and diagonals to product
c
         do 10 vector=oldnum+1,oldnum+added
            if (.not.diag) then
               call iosys('read real diagonals from bliu',-1,v1,0,' ')
               call iosys('read real "trial:'//itoc(vector)//
     #                    '" from bliu',-1,v2,0,' ')
               do 5 i=1,n
                  v2(i)=v2(i)*v1(i)
    5          continue
            end if
            call iosys('read real "product:'//itoc(vector)//
     #                 '" from bliu',-1,v1,0,' ')
            if (.not.diag) then
               do 6 i=1,n
                  v1(i)=v1(i)+v2(i)
    6          continue
c
c              ----- write total product vector back to disc -----
c
               call iosys('write real "product:'//itoc(vector)//
     #                    '" to bliu',-1,v1,0,' ')
            end if
c
c           ----- add the next row to the small matrix -----
c
            ij=vector*(vector-1)/2
            do 9 trial=1,vector
               call iosys('read real "trial:'//itoc(trial)//
     #                    '" from bliu',-1,v2,0,' ')
               ij=ij+1
               a(ij)=sdot(n,v1,1,v2,1)
    9       continue
   10    continue
c
c        ----- save the small matrix -----
c
         oldnum=oldnum+added
         nnpnum=oldnum*(oldnum+1)/2
         if(debug) then
            write(iout,*) ' bliu:small matrix'
            ij=1
            do 12 i=1,oldnum
               write(iout,*) 'i:',i,(a(k),k=ij,ij+i-1)
               ij=ij+i
   12       continue
         endif
         call iosys('write real "small matrix" to bliu',
     #              maxvec*(maxvec+1)/2,a,0,' ')
c
c        ----- solve the small eigenproblem -----
c
         call rsp(maxvec,oldnum,nnpnum,a,root,1,eigvec,v1,v2,error)
         if(error.ne.0) then
            call lnkerr('bliu: error code from rsp')
         end if
c
c        ----- print out results for this cycle -----
c
         cycle=cycle+1
         if (prnt) then
            write (iout,11) cycle,oldnum,totnum
   11       format (t6,' cycle',i4,' using ',i4,' vectors ',
     #           'for a total of',i4,/,
     $           t10,'root',t34,'energy',t50,'convergence')
         end if
         maxdif=0.0d+00
         numdif=0
c
c        ----- form the new "best" vectors -----
c
         pt=0
         do 20 nroot=1,min(mxroot,oldnum)
            call rzero(v1,n)
            do 16 i=1,oldnum
               call iosys('read real "trial:'//itoc(i)//
     #                    '" from bliu',-1,v2,0,' ')
               call saxpy(n,eigvec(i,nroot),v2,1,v1,1)
   16       continue
            call iosys('write real "vector:'//itoc(nroot)//'" to bliu',
     #                 n,v1,0,' ')
            if(debug) then
               write(iout,*) 'vector:',nroot,(v1(i),i=1,n)
            end if
c
c           ----- and form hc -----
c
            call rzero(v1,n)
c
            do 18 i=1,oldnum
               call iosys('read real "product:'//itoc(i)//
     #                    '" from bliu',-1,v2,0,' ')
               call saxpy(n,eigvec(i,nroot),v2,1,v1,1)
   18       continue
            call iosys('write real "hc:'//itoc(nroot)//'" to bliu',
     #                 n,v1,0,' ')
            if(debug) then
               write(iout,*) 'hc:',nroot,(v1(i),i=1,n)
            end if
            call iosys('read real "vector:'//itoc(nroot)//'" from bliu',
     $                 n,v2,0,' ')
c
c           ----- form (h-e)c -----
c
            do 17 i=1,n
               v1(i)=v1(i)-root(nroot)*v2(i)
 17         continue
            if(debug) then
               write(iout,*) '(h-e)c:root ',nroot,(v1(i),i=1,n)
            end if
c
            call iosys('write real "(h-e)c:'//itoc(nroot)//'" to bliu',
     #                 n,v1,0,' ')
c
c           ----- check on convergence, and print results -----
c
            sqcdif=sqrt(sdot(n,v1,1,v1,1))
c
            if (sqcdif.gt.maxdif) then
               maxdif=sqcdif
               numdif=nroot
            end if
c
            if (sqcdif.gt.cnvgnc) then
               pt=pt+1
               ptroot(pt)=nroot
            end if
c
            if (prnt) then
               write (iout,19) nroot,root(nroot)+efixed,sqcdif
   19          format (t10,i4,t24,f15.9,t50,2f15.9)
            end if
   20    continue
c
         ntodo=pt
c
c        ----- check if all roots are converged -----
c
         if (maxdif.lt.cnvgnc) then
            if (mxroot.eq.nroots) then
               status='converged'
               return
            else
c
c              ----- create vectors,etc for all of the roots
c
               mn=mxroot+1
               mx=min(nroots,mn+nattim-1)
c
c              ----- form the new "best" vectors -----
c     
               do 120 nroot=1,min(nroots,oldnum)
                  call rzero(v1,n)
                  do 116 i=1,oldnum
                     call iosys('read real "trial:'//itoc(i)//
     #                    '" from bliu',-1,v2,0,' ')
                     call saxpy(n,eigvec(i,nroot),v2,1,v1,1)
 116              continue
                  call iosys('write real "vector:'//itoc(nroot)//
     $                 '" to bliu',n,v1,0,' ')
c
c                 ----- and form hc -----
c
                  call rzero(v1,n)
c     
                  do 118 i=1,oldnum
                     call iosys('read real "product:'//itoc(i)//
     #                    '" from bliu',-1,v2,0,' ')
                     call saxpy(n,eigvec(i,nroot),v2,1,v1,1)
 118              continue
                  call iosys('write real "hc:'//itoc(nroot)//
     $                 '" to bliu',n,v1,0,' ')
                  call iosys('read real "vector:'//itoc(nroot)//
     $                 '" from bliu',n,v2,0,' ')
c
c                 ----- form (h-e)c -----
c
                  do 117 i=1,n
                     v1(i)=v1(i)-root(nroot)*v2(i)
 117              continue
c
                  call iosys('write real "(h-e)c:'//itoc(nroot)//
     $                 '" to bliu',n,v1,0,' ')
c
c                 ----- check on convergence, and print results -----
c
                  sqcdif=sqrt(sdot(n,v1,1,v1,1))
                  if (prnt) then
                     write (iout,19) nroot,root(nroot)+efixed,sqcdif
                  end if
                  if (sqcdif.gt.maxdif) then
                     maxdif=sqcdif
                     numdif=nroot
                  end if
 120           continue
c     
c              ----- set up the vectors and products for the 'solved' 
c                     roots, and reform the small matrix 
c
               do 140 nroot=1,nroots
                  call iosys('read real "vector:'//itoc(nroot)//
     $                 '" from bliu',n,v1,0,' ')
                  call iosys('write real "trial:'//itoc(nroot)//
     $                 '" to bliu',n,v1,0,' ')
                  call iosys('read real "hc:'//itoc(nroot)//
     $                 '" from bliu',n,v1,0,' ')
                  call iosys('write real "product:'//itoc(nroot)//
     $                 '" to bliu',n,v1,0,' ')
c
c                 ----- add the next row to the small matrix -----
c
                  ij=nroot*(nroot-1)/2
                  do 130 trial=1,nroot
                     call iosys('read real "trial:'//itoc(trial)//
     #                    '" from bliu',-1,v2,0,' ')
                     ij=ij+1
                     a(ij)=sdot(n,v1,1,v2,1)
 130              continue
 140           continue
c
c              ----- save the small matrix -----
c
               oldnum=nroots
               call iosys('write real "small matrix" to bliu',
     #                     maxvec*(maxvec+1)/2,a,0,' ')
               mnroot=mn
               mxroot=mx
               ntodo=mxroot-mnroot+1
               do 402 i=1,ntodo
                  ptroot(i)=mnroot+i-1
 402           continue
            end if
         end if
c
         added=0
         newtr=0
         troot=0
c---------------------------------------------------------- 'new trial'
      else if (oper.eq.'new trial') then
c
c        ----- create a new trial vector, start by forming hc/(e-hii)
c
         status=' '
   21    continue
         troot=troot+1
         nroot=ptroot(troot)
c
c        ----- check if we are done with new trial vectors -----
c
         if (troot.gt.ntodo) then
            status='done'
            if (newtr.gt.0) then
               return
            else
               nroot=numdif
            end if
         end if
c
         call iosys('read real "(h-e)c:'//itoc(nroot)//'" from bliu',
     #                                         -1,v1,0,' ')
c
c        ----- divide by the perturbation denominator -----
c
         eroot=root(nroot)
         call iosys('read real diagonals from bliu',-1,v2,0,' ')
         do 26 i=1,n
            v2(i)=eroot-v2(i)
            if (abs(v2(i)).lt.nrzero) v2(i)=1.0d+00
   26    continue
         do 27 i=1,n
            v1(i)=v1(i)/v2(i)
   27    continue
c
c        ----- and orthonormalize to previous vectors -----
c
         do 35 i=1,oldnum+newtr
            call iosys('read real "trial:'//itoc(i)//'" from bliu',
     #                 -1,v2,0,' ')
            ovrlap=sdot(n,v1,1,v2,1)
            call saxpy(n,-ovrlap,v2,1,v1,1)
            norm=sqrt(sdot(n,v1,1,v1,1))
            if (norm.lt.thresh.and.status.ne.'done') go to 21
   35    continue
c
         if (status.eq.'done') troot=ntodo+1
         call sscal(n,1.0d+00/norm,v1,1)
c
c        ----- second orthogonalization -----
c
         do 38 i=1,oldnum+newtr
            call iosys('read real "trial:'//itoc(i)//'" from bliu',
     #                 -1,v2,0,' ')
            ovrlap=sdot(n,v1,1,v2,1)
            call saxpy(n,-ovrlap,v2,1,v1,1)
            norm=sqrt(sdot(n,v1,1,v1,1))
   38    continue
         call sscal(n,1.0d+00/norm,v1,1)
c
c        ----- return this vector with an 'all okay' -----
c
         newtr=newtr+1
         call iosys('write real "trial:'//itoc(oldnum+newtr)//
     #              '" to bliu',n,v1,0,' ')
         status='trial'
c
c--------------------------------------------------------- 'get vector'
c
      else if (oper.eq.'get vector') then
         if (i3.le.0.or.i3.gt.nroots) then
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
         call iosys('read real "vector:'//itoc(i3)//'" from bliu',
     #              -1,v1,0,' ')
c----------------------------------------------------------- 'finish'
      else if (oper.eq.'finish') then
         call iosys('destroy bliu',0,0,0,' ')
c------------------------------------------------------------- errors
      else
         call lnkerr('bad call to bliu')
      end if
c
c
      return
      end

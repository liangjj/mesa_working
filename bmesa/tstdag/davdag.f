*deck @(#)davdag.f	1.5  7/30/91
      subroutine davdag(oper,status,d,v,hv,resid,n,maxvec,i3,
     1                  i4,i5,small,root,eigvec,cvn)
c
c***begin prologue     davdag
c***date written       960302   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords           eigenvalue, large matrices, davidson algorithm,
c                      liu algorithm, diagonalization
c***author             paul saxe (lanl), barry schneider(nsf)
c***purpose            iterative solution of an eigenvalue problem for
c                      one or a few of the lowest eigenvalues and
c                      vectors.
c***description
c
c    davdag is intended for extracting one or a few of the lowest
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
c   call davdag('initialize',print,d,thresh,r1,r1,n,
c                maxvec,nroots,iout,nattim,efixed,cnvgnc,r1)
c
c  on input:
c
c    print     character
c              'noprint' for no printing.
c
c    d         real (n)
c              the diagonal of the matrix.
c
c    thresh    real
c              the cutoff on the norm of the schmidt orthogonalized
c              potential trial vector. those smaller than "thresh" will
c              be thrown away.
c
c    r1        real variable which is not used in this call
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
c----------------------------------------------------------------------
c
c  'with vectors'
c
c     call davdag('with vectors',char,r1,r1,r1,r1,j1,j1,i3,j1,j1,
c                  r1,r1,r1)
c
c  on input:
c
c    i3       integer
c             number of new trial vectors entered.
c
c    r1       real and unused variable
c
c    j1       integer and unused variable
c
c----------------------------------------------------------------------
c
c  'solve'
c
c    call davdag('solve',status,d,v,hv,resid,n,maxvec,j1,j1,j1,
c                 small,root,eigvec,cvn)
c
c  on input:
c
c   d         real
c             diagonals of matrix
c
c   v         real
c             basis vectors which are both transformed
c             and augmented from iteration to iteration.
c
c   hv        real
c             hamiltonian times matrix.  as with the vectors, these are
c             transformed and augmented from iteration to iteration.
c             the subset of hvecs coming from the immediate previous
c             iteration do not contain the diagonals.
c
c   resid     real
c             used as scratch during the solve.
c
c   n         integer
c             the size of a side of the matrix, and the length of the
c             vectors.
c
c   maxvec    integer
c             the maximum number of expansion vectors permitted.
c
c   n         integer
c             size of matrix,
c
c   small     real (maxvec*(maxvec+1)/2)
c             the small hamiltonian matrix.  this is transformed and 
c             augmented from iteration to iteration.
c
c   j1        integer
c             unused
c
c  on output:
c
c   status    character
c             either 'continue' or 'converged'.
c
c   hv        real
c             hamiltonian times matrix.  as with the vectors, these are
c             transformed and augmented from iteration to iteration.
c             upon exit the hvecs from all previous iterations have been
c             transformed and in addition, the diagonals have been added.
c             the hvecs which are added at this iteration do not contain the 
c             effect of the diagonals.
c
c   resid     real
c             returns the unconverged residuals which are needed to construct
c             the next trial vectors.
c 
c   root      real (nroots)
c             the eigenvalues of small matrix.
c
c   eigvec    real (maxvec,nroots)
c             the eigenvectors of the effective eigenvalue problem.
c
c   j1        integer
c             number of unconverged residuals
c 
c----------------------------------------------------------------------
c
c  'new trial'
c
c     call davdag('new trial',status,trial,v1,n,maxvec,i1,i2,i3,r1,
c                  r1,root,eigvec)
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
c            either 'trial' if there is another trial vector, or 'done'
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
c     call davdag('get vector',status,vector,r1,n,maxvec,iroot,i1,i2,
c                  r2,r3,r4)
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
c     call davdag('finish',c1,r1,r2,i1,i2,i3,i4,i5,r3,r4,r5)
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
c***end prologue       davdag
c
c
      implicit integer (a-z)
c
      parameter (mxroot=500)
c
      real*8 d(n), v(n,*), hv(n,*), resid(n,*)
      real*8 small(maxvec,*), root(maxvec), eigvec(maxvec,*), cvn(*)
      real*8 sdot, maxdif, eroot, nrzero, norm, thresh, efixed, test
      real*8 ovrlap, cnvgnc
      character*(*) oper, status
      character*4 itoc
      character*3 truth
      character*80 title
      logical prnt, diag
      logical debug
      integer ptroot(mxroot)
c
      common /io/ inp,itape6
c
      data debug/.false./
      data nrzero / 1.0d-06 /
c
      save oldnum, added, totnum, nattim, numcnv
      save diag, cycle, maxdif, cnvgnc, efixed, nroots
      save troot, ntodo, newtr, thresh
      save nrzero, ptroot, prnt, iout, lwr, upr
c
c
c**********************************************************************c
c         decide what to do based on "oper" passed in
c**********************************************************************c
c
c                      'initialize'
c
      if (oper(1:10).eq.'initialize') then
c**********************************************************************c
c              save local copies of important variables
c**********************************************************************c
         nroots=i3
         iout=i4
         nattim=i5
         thresh=v(1,1)
         efixed=small(1,1)
         cnvgnc=root(1)
         prnt=status.ne.'noprint'
         diag=oper(12:).eq.'with diagonals'
c*********************************************************************c
c               initialize some local variables
c*********************************************************************c
         call izero(ptroot,mxroot)
         oldnum=0
         newtr=0
         added=0
         cycle=0
         totnum=0
         numcnv=0 
         lwr=1
         upr=nattim
c         
c                       'with vectors'
c
      else if (oper.eq.'check matrix size') then
c*********************************************************************c
c           check for exceeding maximum number of expansion vectors
c*********************************************************************c
c
         num2ad=i3
         if (oldnum+added+num2ad.gt.maxvec) then
             call lnkerr('davdag exceeded maximum number, '
     #                    //itoc(maxvec) //' of expansion vectors')
         end if
c*********************************************************************c
c                      update numbers
c*********************************************************************c
c
         added=added+num2ad
         totnum=totnum+num2ad
c         
c                            'solve'
c
      else if (oper.eq.'solve') then
         status='continue'
c*********************************************************************c
c            set up and solve the small set of equations
c            note that the number of expansion vectors
c            is equal to or greater than the number of roots
c*********************************************************************c
         if (oldnum.eq.0) then
             call rzero(small,maxvec*maxvec)
         end if
c
c*********************************************************************c
c        first, add product of trial and diagonals to product
c*********************************************************************c
c
         do 10 vector=oldnum+1,oldnum+added
            do 20 i=1,n
               hv(i,vector) = d(i)*v(i,vector) + hv(i,vector)
   20       continue
c*********************************************************************c
c                add the next row to the small matrix
c*********************************************************************c
c
            do 30 trial=1,vector
               small(vector,trial)=sdot(n,v(1,trial),1,hv(1,vector),1)
               small(trial,vector)=small(vector,trial)
   30       continue
   10    continue
c   
c*********************************************************************c
c                  print the the small matrix 
c*********************************************************************c
c
         oldnum=oldnum+added
         nnpnum=oldnum*(oldnum+1)/2
         if(debug) then
            title=' davidson:small matrix'
            call prntrm(title,small,oldnum,oldnum,maxvec,maxvec,iout)       
         endif
c
c*********************************************************************c
c                   solve the small eigenproblem
c*********************************************************************c
         call mmove(small,eigvec,oldnum,oldnum,maxvec,maxvec)
         call tred2(maxvec,oldnum,eigvec,root,resid,eigvec)
         call tql2(maxvec,oldnum,root,resid,eigvec,error)
         if(error.ne.0) then
            call lnkerr('davidson: error code from rsp')
         end if
c
c*********************************************************************c
c                   print results for this cycle
c*********************************************************************c
c
         cycle=cycle+1
         if (prnt) then
            write (iout,1) cycle, oldnum, totnum
            write(iout,2)
         end if
c
c*********************************************************************c
c               at this point we have a set of oldnum vectors
c               of which the lowest nroots vectors are what we 
c               want to converge.
c*********************************************************************c   
c
c
c*********************************************************************c
c               form the new "best" vectors.
c               copy them back to their proper vector slots.
c               at this point there are oldnum vectors which are        
c               the results of diagonalizing the small eigenproblem.
c               this number is not necessarily the number of roots
c               which we need.
c               then do the same for the products of h on the vectors.     
c*********************************************************************c
c
         call ebcxx(resid,v,eigvec,n,oldnum,oldnum,n,n,maxvec)
         call copy(resid,v,n*oldnum)
c
         call ebcxx(resid,hv,eigvec,n,oldnum,oldnum,n,n,maxvec)
         call copy(resid,hv,n*oldnum)
c
c*********************************************************************c
c                print results if requested
c*********************************************************************c
c
         if(debug) then
            title='new best vectors for cycle = '//itoc(cycle)
            call prntrm(title,v,n,oldnum,n,maxvec,iout)
         end if
c
         if(debug) then
            title='new best products for cycle = '//itoc(cycle)
            call prntrm(title,hv,n,oldnum,n,maxvec,iout)
         end if
c
c*********************************************************************c
c                we now form the residuals (h-e)c for the
c                all of the eigenvectors. 
c*********************************************************************c
c
         do 40 nroot=1,oldnum
            do 50 i=1,n
               resid(i,nroot)=hv(i,nroot)-root(nroot)*v(i,nroot)
 50         continue
            cvn(nroot)=sqrt(sdot(n,resid(1,nroot),1,resid(1,nroot),1))
 40      continue
         if(debug) then
            title=' (h-e)c '
            call prntrm(title,resid,n,oldnum,n,maxvec,iout)
         end if
c
c**********************************************************************c
c              check on convergence, and print results
c**********************************************************************c
         status='continue'
         do while(status.eq.'continue')
            pt=0
            ncnv=0
            maxdif=0.d0
            lower=lwr
            if(oldnum.ge.upr) then
               upper=upr
            else
               upper=oldnum
            endif               
            do 60 nroot=lower,upper
               if (cvn(nroot).gt.cnvgnc) then
                   pt=pt+1
                   ptroot(pt)=nroot
                   truth='no'
               else
                   ncnv=ncnv+1
                   truth='yes'
               end if
               if (cvn(nroot).gt.maxdif) then
                   maxdif=cvn(nroot)
               end if
               if (prnt) then
                  write (iout,3) nroot, root(nroot)+efixed, 
     1                                  cvn(nroot), truth
               endif
   60       continue
c
c**********************************************************************c
c          pt now tells us how may roots are unconverged and
c          the array ptroot has which roots they are.                  
c**********************************************************************c
c 
            ntodo=upper-lower+1-ncnv
            i3=ntodo
c
            if (prnt) then
               endfile iout
               backspace iout
            end if
            if(ntodo.eq.0) then
               numcnv=numcnv+upper-lower+1
               if(numcnv.eq.nroots) then
c**********************************************************************c 
c                 all roots converged
c**********************************************************************c      
                  status='converged'
                  return
               else
                  lwr=upr+1
                  upr=min(upr+nattim,nroots)
               endif
            else
               status='more'
            endif
         enddo               
         if (prnt) then
             endfile iout
             backspace iout
         end if
c
c**********************************************************************c
c               reform the small matrix in the new basis 
c**********************************************************************c
            do 70 nroot=1,oldnum
               do 80 trial=1,nroot
                  small(nroot,trial)
     1                               =
     2                           sdot(n,hv(1,nroot),1,v(1,trial),1)
                  small(trial,nroot) = small(nroot,trial)
c              write(iout,*) 'i',nroot,'j',trial,small(nroot,trial)
  80           continue
  70        continue
c
         added=0
         newtr=0
         troot=0
c         
c                      'new trials'
c
      else if (oper.eq.'new trials') then
c
c**********************************************************************c
c            create new trial vectors. start by forming hc/(e-hii)
c            for the unconverged roots. schmidt orthogonalize it to the
c            other vectors.  after the schmidt process check that the 
c            vector is not linearly dependent on the other vectors.
c            keep it or discard it according to the size of the norm.      
c**********************************************************************c
c
         status=' '
         do 90 nroot=1,ntodo
            pt=ptroot(nroot)                  
c         
c*********************************************************************c
c           form the trial vector by taking the residual for this
c           root and dividing by the perturbation denominator
c*********************************************************************c
c          
c            write(iout,*) 'input passed vector'
c            write(iout,*) (resid(ii,pt),ii=1,n)
            eroot=root(pt)
            do 100 i=1,n
               test=eroot-d(i)
               if (abs(test).ge.nrzero) then
                   resid(i,pt)=resid(i,pt)/test
               else     
                   resid(i,pt)=1.0d+00
               endif                
  100       continue
c
c
c*********************************************************************c
c          orthonormalize to previous vectors
c*********************************************************************c
c
            do 200 i=1,oldnum+newtr
               ovrlap=sdot(n,resid(1,pt),1,v(1,i),1)
               call saxpy(n,-ovrlap,v(1,i),1,resid(1,pt),1)
               norm=sqrt(sdot(n,resid(1,pt),1,resid(1,pt),1))
               if (norm.lt.thresh) go to 90
  200       continue
c
            call sscal(n,1.0d+00/norm,resid(1,pt),1)
c
c
c            do 300 i=1,oldnum+newtr
c               ovrlap=sdot(n,resid(1,pt),1,v(1,i),1)
c               call saxpy(n,-ovrlap,v(1,i),1,resid(1,pt),1)
c               norm=sqrt(sdot(n,resid(1,pt),1,resid(1,pt),1))
c  300       continue
c
c            call sscal(n,1.0d+00/norm,resid(1,pt),1)
c
c            write(iout,*) 'second'
c            write(iout,*) 'output passed vector'
c            write(iout,*) (resid(ii,pt),ii=1,n)
            newtr=newtr+1
            call copy(resid(1,pt),v(1,oldnum+newtr),n)
   90    continue
c         do 999 i=1,oldnum+newtr
c            do 998 j=1,i
c               ovrlap=sdot(n,v(1,i),1,v(1,j),1)
c               write(iout,22) i,j,ovrlap
c 998        continue
c 999     continue   
c 22      format(1x,'i=',i3,1x,'j=',i3,1x,'overlap=',e15.8)
         if(newtr.eq.0) then
            write(iout,4)
            status='none left'
         else
             status='done'                              
             write(iout,5) newtr
         endif 
         i3=newtr   
c
c                     'get vector'
c
      else if (oper.eq.'get vector') then
         if (i3.le.0.or.i3.gt.nroots) then
            status='error'
            return
         end if
c
         if (maxdif.ge.cnvgnc) then
            status='not converged'
         else
            status='ok'
         end if
c                          'finish'
c
      else if (oper.eq.'finish') then
         write(iout,99)
   99 format(/,1x,'eigenvalue calculation finished')         
      else                    
         call lnkerr('bad call to davidson')
      end if
c
c
      return
    1 format (/,5x,' cycle = ',i4,1x,' number of vectors = ',i4,1x,
     1                     ' total number of vectors = ',i4)
    2 format(/,5x,'     root     ',3x,'     energy    ',7x,
     1            '      rms    ',3x,'   converged   ')
    3 format(9x,i3,10x,e15.8,6x,e15.8,8x,a3)
    4 format(/,15x,'there are no more trial vectors available')
    5 format(/,15x,'adding ',i3,' new vectors to the space')
      end



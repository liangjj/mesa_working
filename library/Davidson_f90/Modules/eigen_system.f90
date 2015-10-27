!***eigen_system_driver
!***begin prologue     eigen_system_driverd
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           iterative eigenvalue
!***author             schneider, barry (nsf)
!***source
!***purpose            find the eigenvalues of a real symmetric matrix
!***                   having diagonal dominance using the davidson
!***                   algorithm.
!***
!***description        the davidson algorithm as programmed here is
!***                   capable of getting the lowest m roots or of
!***                   climbing the eigenvalue tree in a piecewise fashion
!***                   by getting the first m1 roots, then the next m2 roots
!***                   then the next.... roots until all are obtained.
!***                   it can often happen that roots will be missed,
!***                   especially if one tries to compute too many
!***                   at each pass.  so, the routine must be used
!***                   intelligently.  it is typically simple to see if a
!***                   root has been missed by inspection.  the method avoids
!***                   diagonalizing a large matrix for the higher roots by
!***                   explicitly orthogonalizing the guesses to the converged
!***                   eigenvectors of the lower roots.
!***
!***     variable      type       size             description
!***     _______       ____       ____             ___________
!***      ibuf          integer 2*lenbuf      buffer to hold non-zero values of
!***                                          hamiltonian indices
!***      rbuf          real    lenbuf        buffer to hold non-zero values of
!***                                          hamiltonian values
!***      diag          real    n             diagonals of hamiltonian

!***      vec           real    n*maxvec      davidson vectors
!***                                          can be same storage as trials
!***      hvec          real    n*maxvec      action of the hamiltonian on vec
!***      b, bwrk       real    maxvec*maxvec small hamiltonian matrix and a
!***                                          copy
!***      eigwrk        real    maxvec        eigenvalues of small matrix
!***      work          real    5*maxvec      scratch array
!***      svec          real    maxvec*maxvec small scratch matrix
!***      resid         real    n*maxvec      residual vectors
!***      eig           real    nroot         converged eigenvalues
!***      rep           real                  nuclear repulsion
!***      cnverg        real                  convergence criterion for a root
!***      thresh        real                  overlap tolerance for accepting
!***                                          a new davidson vector
!***      n             integer               matrix size
!***      nroot         integer               number of roots to get
!***      ntrial        integer               number of trial vectors
!***      nattim        integer               number of roots to converge at a
!***                                          pass
!***      maxit         integer               maximum size of davidson space
!***      maxvec        integer               maximum number of vectors
!***                                          available for storage
!***      lenbuf        integer               size of buffer for non-zero
!***                                          labels and matrix elements
!***      header        character             labels for reading buffers
!***      prnt          logical               print flags
!***references

!***routines called
!***end prologue       eigen_system_driverd

  SUBROUTINE eigen_system_driver_d(ibuf,rbuf,diag,eig,vec,hvec,resid,b,bwrk,  &
                                   eigwrk,work,svec,rep,cnverg,thresh,n,nroot,  &
                                   ntrial,nattim,maxit,maxvec,lenbuf, headr,prnt,code)
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size,0:maximum_number_of_davidson_vectors)       :: vec
  REAL*8, DIMENSION(matrix_size,0:maximum_number_of_davidson_vectors)       :: h_vec
  REAL*8, DIMENSION(matrix_size,0:maximum_number_of_davidson_vectors)       :: residual
REAL*8, INTENT(IN)                       :: eig(nroot)
REAL*8, INTENT(IN OUT)                   :: b(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: bwrk(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: eigwrk(maxvec)
REAL*8, INTENT(IN OUT)                   :: work(5*maxvec)
REAL*8, INTENT(IN OUT)                   :: svec(maxvec,*)
REAL*8, INTENT(IN)                       :: rep
REAL*8, INTENT(IN OUT)                   :: cnverg
REAL*8, INTENT(IN OUT)                   :: thresh
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nroot
INTEGER, INTENT(IN)                      :: ntrial
INTEGER, INTENT(IN)                      :: nattim
INTEGER, INTENT(IN)                      :: maxit
INTEGER, INTENT(IN)                      :: maxvec
INTEGER, INTENT(IN)                      :: lenbuf
CHARACTER (LEN=*), INTENT(IN OUT)        :: headr(3)
LOGICAL, INTENT(IN)                      :: prnt(11)
CHARACTER (LEN=*), INTENT(IN)            :: code
IMPLICIT INTEGER (a-z)

REAL*8  error, sdot, tmp
REAL*8 maxerr, zero, one, nrzero, eci
REAL*4 time(2), secnds, delta(10), total
LOGICAL :: incore
CHARACTER (LEN=5) :: itoc
CHARACTER (LEN=8) :: cntrl, pascod, precon
CHARACTER (LEN=80) :: title
#ifdef decpointer
INTEGER*8 pre
#END IF decpointer
#ifdef sgipointer
INTEGER*4 pre
#END IF sgipointer





DIMENSION  pascod(2)

COMMON/io/inp, iout
DATA zero, one/ 0.d0, 1.d0/
DATA nrzero / 1.0D-06 /
DATA incore/.false./
DATA precon /'diagonal' /

!
!     find out how many passes are needed to get all the roots
!
  number_of_passes = number_of_roots/number_of_roots_per_pass
  number_left = number_of_roots - number_of_passes * number_of_roots_per_pass
  IF(number_left /= 0) THEN
     number_of_passes = number_of_passes + 1
   ELSE
     number_left = number_of_roots_per_pass
  END IF
!
! Read in the Trial Solutions
!
  Call iosys('read real guess_eigenvectors from guess_vectors',             &
              guess_size*matrix_size,vec,0,' ')
  Call iosys('open davidson_vectors as new',0,0,0,'Davidson_Vectors.out')
  converged_davidson_vectors = 0
  root_0=0
  root_n=0
  number_2_do=number_of_roots_per_pass
  number_of_vectors = number_of_guess_vectors
  delta(:)=0.0
  DO  i_pass =1, number_of_passes 
      time(1)=secnds(0.0)
      IF(i_pass == number_of_passes ) THEN
         number_2_do = number_left
      END IF
      root_0 = root_n + 1
      root_n = root_n + number_2_do
      WRITE(iout,2) i_pass, root_0, root_n
!
!     Initialize the vectors for this pass
!  
!-----------------------------------------------------------------------c
!                                                                       c
!                    Initialization Section                             c
!                                                                       c
!-----------------------------------------------------------------------c
  
!                      This is the preparation phase.
!        We have some trial vectors for these roots and some, possibly none,       &
!        converged Davidson roots. Then we take the trial vectors and              &
!        orthogonalize them to the converged Davidson vectors.                     &
!        After that we take those and orthonormalize them to get the               &
!        initial set of trial vectors for this set of roots.
!
      first_vector = 0
      number_expected = number_of_vectors      
      Call Initialize(vec(:,0:number_expected - 1 ),residual,number_expected,      &
                      actual_number,converged_davidson_vectors)
!
!        OK, we have a set of actual_number /=0 vectors after the Schmidt Process
!        and we can proceed
!
      IF(prnt(2)) THEN
         title='initial vectors for pass = '//itoc(i_pass)
         CALL prntfm(title,vec(:,0:actual_number - 1),matrix_size,actual_number,   &
                     matrix_size,actual_number,iout)
      END IF
!
      last_vector = actual_number - 1
      number_in = last_vector - first_vector + 1
      it = 0
      dvd_error=1.d+10
      WRITE(iout,4) dvd_error
      convergence_control='continue'
      eigen_values = 0.d0
      working_eigen_values = 0.d0
      time(2)=secnds(0.0)
      delta(1)=delta(1) + time(2) - time(1)
      DO WHILE ( convergence_control == 'continue'.and.it <               &
                                         maximum_number_of_iterations )
!----------------------------------------------------------------------c
!                                                                      c
!                    Iteration Sequence                                c
!                                                                      c
!     iteration is continued until all of the roots are converged      c
!     or if convergence is not achieved some preset maximum number of  c
!     iterations are performed.                                        c
!                                                                      c
!----------------------------------------------------------------------c
!
!        initialize the effect of the hamiltonian on these vectors.
!        time(1) = secnds(0.0)
         call packed_symmetric_matrix_on_vector                                              &
                                       (h_buf_d,                                             &
                                        ibuf,                                                &
                                        diag_d,                                              &
                                        vec_d(:,first_vector:last_vector),                   &
                                        h_vectors_d(:,first_vector:last_vector),             &
                                        non_zero,number_in)
         time(2)=secnds(0.0)
         delta(2)=delta(2) + time(2) - time(1)
!  
!        Construct/Update and Diagonalize the Small Hamiltonian Matrix.  The small eigenvector
!        matrix is in h_mat_work.
!  
         time(1)=secnds(0.0)         
         Call Solve_Small_Eigenvalue_Problem(h_mat_d,h_mat_work_d,working_eigen_values,it,   &
                                             first_vector,last_vector)            
         time(2)=secnds(0.0)
         delta(3)=delta(3) + time(2) - time(1)
!
!        Transform Relevant Quantities To The New Basis Defined By The Diagonalization
!        Then Compute The Residuals and Check For Convergence.  The converged
!        Eigenvalues Are Stored And The Eigenpairs Written To Disk.  The
!        Unconverged Residuals Are Moved From Their Current Positions To The Leading
!        Positions In The List Of New Vectors.
    
         WRITE(iout,5) it, number_in
         time(1)=secnds(0.0)
         number_converged = 0
         number_unconverged = 0
        CALL Form_Residuals(vec_d,h_vectors_d,residual_d,svec,eig,eigwrk,b,bwrk,rep,  &
        davidson_convergence,maximum_error,n,number_in,number_2_do,number_converged,    &
        number_unconverged, maxvec,it,prnt(5),pascod)
    time(2)=secnds(0.0)
    delta(4) = delta(4) + time(2) - time(1)
    
    
!           check to see if all num2do roots are converged
!           if so, we are done and can print results and then proceed to
!           the next set of roots.
    
    IF(con == num2do) THEN
      
      nlft=MIN(nend-num2do,ntrial)
      remain=ntrial-nlft
      WRITE(iout,6)
      DO  i=1,num2do
        actual=i+converged_davidson_vectors
        WRITE(iout,7) actual, eig(actual)+rep
      END DO
      
!              begin next set of roots with the best available
!              vectors.  these are the the remaining vectors coming
!              from the diagonalization of the small matrix augmented
!              with the remaining trials.
      IF(actual < nroot) THEN
        WRITE(iout,8) nlft, remain
        IF(nlft > 0) THEN
          WRITE(iout,*) ' Doing the vectors moves'
          CALL copy(vec(1,num2do+1),vec,n*nlft)
          IF(remain > 0) THEN
            CALL iosys('read real trials from bliu',  &
                n*remain,vec(1,nlft+1),n*nlft,' ')
          END IF
        ELSE
          CALL iosys('read real trials from bliu',n*ntrial, vec,0,' ')
        END IF
      END IF
      numvec=ntrial
      convergence_control='finished'
    ELSE
      
!           all roots are not converged.  set the error to the largest
!           current error and either restart or continue the
!           iteration sequence.
      
      dvd_error=MIN(dvd_error,maximum_error)
      
!           how many new vectors could be added in principle
      
      numnew = maxvec - nend
      
!           how many will we add
      
      addvec = MIN(numnew,uncon)
      
!           check if the number of old plus new vectors will exceed
!           the maxvec limit to determine if a restart is needed.
      
      chkno = nend + addvec
      IF(chkno >= maxvec) THEN
        
!              the maximum number of vectors is exceeded.
        
        nend = MIN(2*num2do,nend)
        WRITE(iout,9) nend
        numnew = maxvec - nend
        addvec = MIN(numnew,uncon)
      END IF
      WRITE(iout,11) addvec
      
!              maximum number of vectors is still within the allowed
!              limits.  add vectors to the set from the unconverged
!              residuals and put them after the current vectors.
      
      IF(addvec /= 0) THEN
        oldnum = nend
        nbeg = nend + 1
        time(1)=secnds(0.0)
        CALL genvec(vec(1,nbeg),resid,diag,eigwrk,eigwrk,dum,  &
            dum,dum,n,dum,pre,addvec,iter, precon,dum,prnt(9))
        time(2)=secnds(0.0)
        delta(4) = delta(4) + time(2) - time(1)
        
!              orthonormalize the new trials to the converged roots
!              and then to the old vectors to get an additional
!              nout vectors.
        
        num=addvec
        nout=num
        time(1)=secnds(0.0)
        IF(converged_davidson_vectors /= 0) THEN
          CALL rdciv(resid,pascod(2),n,converged_davidson_vectors)
          CALL abschm(resid,vec(1,nbeg),thresh,n,converged_davidson_vectors,  &
              addvec,num,.true.,.false.)
        END IF
        nend=nend+num
        IF(num /= 0) THEN
          CALL gschmt(vec,thresh,n,nbeg,nend, nout,.true.,.false.)
        END IF
        IF(num == 0.OR.nout == 0) THEN
          WRITE(iout,12) num, nout
          RETURN
        END IF
        time(2)=secnds(0.0)
        delta(5) = delta(5) + time(2) - time(1)
        
        IF(nout == 0) THEN
          
!                 no more vectors write out unconverged results.
          
          WRITE(iout,13)
          RETURN
        END IF
        nend=nbeg+nout-1
        IF(prnt(10)) THEN
          CALL tstovl(vec,n,nend)
        END IF
      END IF
      
!           we now have either a reconstituted set of vectors or
!           an augmented set of vectors.
      
!              operate with hamiltonian on these vectors
      
      time(1)=secnds(0.0)
      title='h on initial vectors iteration-'//itoc(iter)
      CALL honv(ibuf,rbuf,diag,vec(1,nbeg),hvec(1,nbeg),  &
          n,nout,headr(1),lenbuf,ntot,incore, title,prnt(3))
      time(2)=secnds(0.0)
      delta(6) = delta(6) + time(2) - time(1)
      
!           update the curent small hamiltonian matrix.
      
      time(1)=secnds(0.0)
      CALL hsmall(b,bwrk,vec,hvec,n,nbeg,nend,maxvec, 'fill',.false.)
      time(2)=secnds(0.0)
      delta(7) = delta(7) + time(2) - time(1)
    END IF
  END DO
  IF(iter >= maxit) THEN
    WRITE(iout,14)
    RETURN
  END IF
  converged_davidson_vectors=converged_davidson_vectors+num2do
END DO
WRITE(iout,15)
DO  i=1,nroot
  eci=eig(i)+rep
  WRITE(iout,7) i, eci
END DO
total=delta(1) + delta(2) + delta(3) + delta(4) +  &
    delta(5) + delta(6) + delta(7)
WRITE(iout,16) ( delta(i), i=1,7), total
RETURN
1    FORMAT(/,1X,'davidson eigenvalue solver',/,10X,  &
    'number of roots               = ',i4,/,10X,  &
    'number of roots at a time     = ',i4,/,10X,  &
    'number of trials              = ',i4,/,10X,  &
    'maximum number of iterations  = ',i4,/,10X,  &
    'maximum number of vectors     = ',i4,/,10X,  &
    'convergence criterion for rms = ',e15.8)
2    FORMAT(/,5X,'pass = ',i3,/,5X, 'processing root = ',i4,' to root =',i4)
3    FORMAT(/,5X,'cannot even begin davidson calculation:',/,5X,  &
    'orthonormalization of initial vectors yields null ' 'set')
4    FORMAT(/,5X,'beginning davidson iterations:',/,5X,  &
    'initial error = 'e15.8)
5    FORMAT (/,5X,'cycle = ',i4,2X,'no. vectors = ',i4)
6    FORMAT(/,30X,'summary',/,17X,'root',22X,'energy')
7    FORMAT(15X,i4,20X,f15.8)
8    FORMAT(/,5X,'number of added vectors from unconverged roots = ',  &
    i5,/,5X, 'number of added vectors from original trial set = ',  &
    i5)
9    FORMAT(/,10X,'***** maximum number of vectors exceeded *****',  &
    /10X,'      contract back to ',1X,i4,' vectors')
11   FORMAT(/,5X,'number of vectors added = ',i4)
12   FORMAT(/,5X,'cannot continue',/,5X, 'num = ',i5,2X,'nout = ',i5,/,5X,  &
    'write out results and return')
13   FORMAT(/1X,'no more orthonormal vectors can be added.',/,5X,  &
    'quit and return unconverged results')
14   FORMAT(/,5X,'iteration limit exceeded.',/,5X, 'quit and return to main')
15   FORMAT(///,30X,'final summary',/,17X,'root',22X,'energy')
16   FORMAT(/,25X,'Time Summary for Eigen_system_driverd',  &
    /,5X,'Initialization         = ',e15.8,  &
    /,5X,'Small diagonalization  = ',e15.8,  &
    /,5X,'Residual Formation     = ',e15.8,  &
    /,5X,'Vector Generation      = ',e15.8,  &
    /,5X,'Schmidt Process        = ',e15.8,  &
    /,5X,'Matrix Vector Multiply = ',e15.8,  &
    /,5X,'Small Matrix Update    = ',e15.8,  &
    /,5X,'                 Total = ',e15.8)
END SUBROUTINE eigen_system_driverd
!***********************************************************************
!***********************************************************************

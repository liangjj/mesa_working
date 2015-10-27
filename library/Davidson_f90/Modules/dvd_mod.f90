!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  MODULE Davidson_Module
!deck Davidson_Modules
!**begin prologue     Davidson_Module
!**date written       030208   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       Davidson_Module
!***********************************************************************
!***********************************************************************
                      MODULE Davidson_Module
                      USE input_output
                      USE Iterative_Global
                      USE Pack_Global
                      IMPLICIT NONE
!***********************************************************************
!***********************************************************************
  INTEGER                                :: matrix_size
  INTEGER                                :: number_of_roots
  INTEGER                                :: number_of_right_hand_sides
  INTEGER                                :: number_of_roots_per_pass
  INTEGER                                :: maximum_number_of_iterations
  INTEGER                                :: maximum_number_of_davidson_vectors
  INTEGER                                :: guess_size
  INTEGER                                :: number_of_guess_vectors
  INTEGER                                :: begin
  INTEGER                                :: end
  REAL*8                                 :: davidson_convergence
  REAL*8                                 :: overlap_tolerance
  REAL*8                                 :: ene
  LOGICAL, DIMENSION(11)                 :: prdvd
  LOGICAL                                :: dvdall
  CHARACTER(LEN=80)                      :: action
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvddat
!***begin prologue     dvddat
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            data entry for davidson routine.
!***
!***references

!***routines called
!***end prologue       dvdddat

  SUBROUTINE dvddat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL                   :: dollar
  LOGICAL                   :: logkey
  INTEGER                   :: intkey
  REAL*8                    :: fpkey
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF( dollar ('$davidson_parameters',card,cpass,card,inp) ) THEN
      overlap_tolerance = fpkey(card,'overlap_tolerance',1.0D-10,' ')
      matrix_size = intkey(card,'matrix_size',2,' ')
      number_of_roots = intkey(card,'number_of_roots',1,' ')
      number_of_roots = intkey(card,'number_of_right_hand_sides',1,' ')
      number_of_roots_per_pass = intkey(card,'number_of_roots_per_pass',number_of_roots,' ')
      guess_size = intkey(card,'guess_size',10*number_of_roots_per_pass,' ')
      guess_size = min(matrix_size,guess_size)
      davidson_convergence = fpkey(card,'davidson_convergence',1.0D-08,' ')
      drop_tol = fpkey(card,'drop_tolerance',1.0d-10,' ')
      maximum_number_of_iterations = intkey(card,'maximum_number_of_iterations',1,' ')
      maximum_number_of_davidson_vectors = intkey(card,'maximum_number_of_davidson_vectors',1,' ')
      number_of_guess_vectors = intkey(card,'number_of_guess_vectors',guess_size,' ')
  END IF
  prdvd(1)=logkey(card,'print=trials',.false.,' ')
  prdvd(2)=logkey(card,'print=vectors',.false.,' ')
  prdvd(3)=logkey(card,'print=h-on-vectors',.false.,' ')
  prdvd(4)=logkey(card,'print=hamiltonian',.false.,' ')
  prdvd(5)=logkey(card,'print=iteration-information',.false.,' ')
  prdvd(6)=logkey(card,'print=residuals',.false.,' ')
  prdvd(7)=logkey(card,'print=transformed-vectors',.false.,' ')
  prdvd(8)=logkey(card,'print=transformed-h-on-vectors', .false.,' ')
  prdvd(9)=logkey(card,'print=new-trial-vectors',.false.,' ')
  prdvd(10)=logkey(card,'print=overlaps',.false.,' ')
  dvdall=logkey(card,'print=all',.false.,' ')
  prdvd(:)=.false.
  IF(dvdall) THEN
     prdvd(:)=.true.
  END IF
  WRITE(iout,1) number_of_roots, number_of_roots_per_pass, matrix_size, guess_size, number_of_guess_vectors,               &
                maximum_number_of_davidson_vectors, maximum_number_of_iterations, overlap_tolerance, davidson_convergence
1 FORMAT(/,15X,'iterative diagonalization information',/,/,5X,  &
                'number of roots                     = ',i4,/,5X,  &
                'number of roots at a time           = ',i4,/,5X,  &
                'size of matrix                      = ',i8,/,5X,  &
                'size of guess matrix                = ',i7,/,5X,  &
                'size of initial vectors             = ',i7,/,5X,  &
                'maximum number of davidson vectors  = ',i7,/,5X,  &
                'max. number of iterations           = ',i7,/,5X,  &
                'overlap tolerance                   = ',e15.8,/,5X,  &
                'convergence criterion               = ',e15.8)
  END SUBROUTINE dvddat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck pack_matrix
!***begin prologue     pack_matrix
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Pack the Matrix.
!***
!***references

!***routines called
!***end prologue       pack_matrix

  SUBROUTINE pack_matrix(matrix,matrix_size)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     pack all non-zero, non-diagonal elements
!
  number=0
  DO i=1,matrix_size
     DO j=1,i-1
        IF(abs(matrix(i,j)).gt.drop_tol) THEN
           number=number
           ibuf(1,nonzro) = i
           ibuf(2,nonzro) = j
           packed_matrix_d(number)=matrix(i,j)
        END IF
     END DO
  END DO 
!
!
  DO i=1,matrix_size
     matrix_diagonal_d(i)= matrix(i,i)
  END DO
  END SUBROUTINE pack_matrix
!***********************************************************************
!***********************************************************************
!deck lindvd.f
!***begin prologue     lindvd
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           iterative eigenvalue
!***author             schneider, barry (nsf)
!***source
!***purpose            iterative linear system solver using
!***                   Davidson algorithm.
!***
!***description        solve the linear set of equations
!***                    A |Psi> = |B>
!***                   where |B> is a vector.
!***references

!***routines called
!***end prologue       lindvd
  SUBROUTINE lindvd(t,list,cnverg,thresh,  &
    n,nrhs,maxit,maxvec,lenbuf,incore,header, prnt)
  INTEGER                                :: begin
  INTEGER                                :: end
  INTEGER                                :: out
  INTEGER                                :: iter
  REAL*8, DIMENSION(:,:)                 :: resid(n,maxvec)
  REAL*8, DIMENSION(:,:)                 :: t(n,*)
  REAL*8                                 :: error
  REAL*8                                 :: sdot
  REAL*8                                 :: max_err
  REAL*8                                 :: zero=0.d0
  REAL*8                                 :: one=1.d0
  WRITE(iout,1) nrhs, maxit, maxvec, cnverg
!-----------------------------------------------------------------------c
!                                                                       c
!                    Initialization Section                             c
!                                                                       c
!-----------------------------------------------------------------------c

!     initialize the set of input vectors as the
!     orthonormalized set of right hand side vectors

  vec_d(:,1) = rhs_d(:)  
  begin = 1
  end = 1
  CALL gschmt(vec_d,overlap_tolerance,matrix_size,begin,end,out,.true.,.false.)
  end=out
  IF(prnt(2)) THEN
    title='initial trial vectors'
    CALL prntfm(title,vec_d(:,begin),matrix_size,out,matrix_size,out,iout)
  END IF
!
!    Multiply the Hamiltonian on the starting vector.
!
  title='h on initial vectors'
  Call packed_symmetric_matrix_on_vector(packed_matrix_d, ibuf,                   &
                                         matrix_diagonal_d, vec_d(:,nbeg),        &
                                         h_vectors_d(:,nbeg), number)
!
!     initialize the small hamiltonian matrix and right hand side.
!
  action='initialize
  CALL bsmall
  IF(prnt(4)) THEN
     title='initial small matrix'
     CALL print_triangle(title,h_mat_work_tri_d,end,iout)
     title='initial small right hand side'
     CALL prntfm(title,small_rhs_work_d,end,1,                                    &
                 maximum_number_of_davidson_vectors,                              &
                 maximum_number_of_davidson_vectors,iout)
  END IF
!----------------------------------------------------------------------c
!                                                                      c
!                    Iteration Sequence                                c
!                                                                      c
!     iteration is continued until the solution is converged           c
!     or if convergence is not achieved some preset maximum number of  c
!     iterations are performed.                                        c
!                                                                      c
!----------------------------------------------------------------------c
  iter=0
  error=1.d+10
  WRITE(iout,3) error
  newrhs=nrhs
  totrhs=0
  DO WHILE ( error > davidson_convergence.AND.iter < maximum_number_of_iterations )
     iter = iter + 1
     WRITE(iout,4) iter, nend
!    Step 1:
!          get solutions of the small matrix.
!          small_rhs_work_d holds the initial matrix which is destroyed.
!
     CALL l_solve(h_mat_work_tri_d,small_rhs_work_d,ipvt)
      IF(prnt(5)) THEN
         title='solution of small matrix iteration = '//itoc(iter)
         CALL prntfm(title,small_rhs_work_d,end,1,                                &
                     maximum_number_of_davidson_vectors,1,iout)
  END IF
  
!     Step 2:
  
!           form the residuals and check for convergence.
!           t contains the all the solution vectors.
!           when lares is finished, the converged solutions are stored
!           on the disk and the unconverged residuals and solutions appear
!           first in the list.  the b and btmp matrices are updated to
!           reflect this and the current value of newrhs changed to the
!           new number.
  
  CALL lares(pvec,hpvec,b,btmp,rhs,energy,cnverg,resid,  &
      maxerr,t,list,n,nend,newrhs,con,uncon,maxvec, iter,prnt(5))
  newrhs=uncon
  totrhs=totrhs+con
  IF(totrhs == nrhs) THEN
    
!           all solutions are converged.  copy them in to rhs and quit
    
    RETURN
  ELSE
    
!           all solutions are not converged.  set the error to the largest
!           current error and continue the iteration sequence
    
    error=MIN(error,maxerr)
    
!           scale the residuals to get the next set of vectors to
!           be added to the krylov sequence.  leave them in resid
!           for the moment.
    
    CALL lavec(energy,diag,resid,n,uncon,iter,prnt(9))
    
!           orthogonalize the scaled residuals to the old vectors
!           to get an additional nout set of linearly independent
!           vectors.
    
    CALL abschm(pvec,resid,thresh,n,nend,uncon,nout, .true.,.false.)
    
!           orthonormalize these vectors and add them to the
!           existing set.
    
    CALL gschmt(resid,thresh,n,1,nout,nfinal,.true.,.false.)
    
!           how many new vectors could be added in principle
    
    numnew = maxvec - nend
    
!           how many will we add
    
    addvec = MIN(numnew,nfinal)
    WRITE(iout,5) nfinal, addvec
    
!           add them
    
    nbeg = nend + 1
    nend=nend+addvec
    IF(nend >= maxvec.AND.addvec == 0) THEN
      WRITE(iout,6)
      CALL lnkerr('vector space exceeded. quit')
    END IF
    CALL copy(resid,pvec(1,nbeg),n*addvec)
    IF(prnt(10)) THEN
      CALL tstovl(pvec,n,nend)
    END IF
    
!           we now have either a reconstituted set of vectors or
!           an augmented set of vectors.
    
!              operate with hamiltonian on these vectors
    
    title='h on initial vectors'
    CALL honv(ibuf,hbuf,diag,pvec(1,nbeg),hpvec(1,nbeg),  &
        n,addvec,header,lenbuf,ntot,incore, title,prnt(3))
    
!           update the curent small hamiltonian matrix and right hand sides.
    
    CALL bsmall(h,htmp,b,btmp,pvec,hpvec,rhs,energy,newrhs,  &
        n,nbeg,nend,maxvec,'fill')
  END IF
END DO
CALL iosys('write integer "size of davidson vector '//  &
    'space" to hamiltonian',1,nend,0,' ')
CALL iosys('write real "davidson vectors" to hamiltonian', nend*n,pvec,0,' ')
RETURN
1    FORMAT(/,1X,'davidson linear system solver',/,10X,  &
    'number of right hand sides    = ',i4,/,10X,  &
    'maximum number of iterations  = ',i4,/,10X,  &
    'maximum number of vectors     = ',i4,/,10X,  &
    'convergence criterion for rms = ',e15.8)
2    FORMAT(/,5X,'cannot even begin davidson calculation:',/,5X,  &
    'orthonormalization of initial vectors yields null' '  set')
3    FORMAT(/,5X,'beginning iterations:',/,5X, 'error initialized at ',e15.8)
4    FORMAT(/,5X,'iteration            = ',i5,  &
    /,5X,'size of vector space = ',i5)
5    FORMAT(/,5X,'number of linearly independent vectors after '  &
    'schmidt step = ',i4,/,5X,'number actually added ' '             = ',i4)
6    FORMAT(/,5X,'size of allowed vector space exceeded. quit')
END SUBROUTINE lindvd
!***********************************************************************
!***********************************************************************
!deck b_small.f
!***begin prologue     b_small
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       b_small
  SUBROUTINE b_small
  REAL*8                                 :: ddot
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: count
  IF(action == 'initialize') THEN
     count = 0
     DO i=1,END
        DO j=1,i
           count = count + 1
           h_mat_tri_d(count) = ddot(matrix_size,vec_d(:,i),1,h_vectors_d(:,j),1)
           h_mat_work_tri_d(count) = h_mat_tri_d(count)
        END DO
        small_rhs_d(i) = ddot(matrix_size,vec_d(:,i),1,rhs_d(:),1)
        small_rhs_work_d(i) = small_rhs_d(i)
     END DO
  ELSE IF(action == 'fill') THEN
     DO i=1,END
        DO  j=begin,END
            count = i*(i-1)/2 + j
            h_mat_tri_d(count) = ddot(matrix_size,vec_d(:,i),1,h_vectors_d(:,j),1)
            h_mat_work_tri_d(count) = h_mat_d(count)
        END DO
        small_rhs_d(i) = ddot(matrix_size,vec_d(:,i),1,rhs_d(:),1)
        small_rhs_work_d(i) = small_rhs_d(i)
     END DO
  END IF
END SUBROUTINE b_small
!***********************************************************************
!***********************************************************************
!deck l_solve
!***begin prologue     l_solve
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           linear system solve
!***author             schneider, barry (nsf)
!***source
!***purpose            driver for direct linear system solve.
!***
!***references

!***routines called
!***end prologue       l_solve
  SUBROUTINE l_solve(matrix,rhs,ipvt)
  REAL*8, DIMENSION(:)                   :: matrix
  REAL*8, DIMENSION(:)                   :: rhs
  INTEGER, DIMENSION(:)                  :: ipvt
  CALL dspsv('u',matrix_size,1,matrix,ipvt,rhs,matrix_size,info)
  IF(info /= 0) THEN
     CALL lnkerr('error from linear solve routine')
  END IF
END SUBROUTINE l_solve
END MODULE Davidson_Module

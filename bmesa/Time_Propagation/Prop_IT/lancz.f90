!deck lancz.f
!***begin prologue     lancz
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           eigenvalues, eigenvectors
!***author             schneider, barry (nsf)
!***source
!***purpose            lanczos diagonalization.
!***
!***description        b v  = (mat - a ) v   - b    v
!***                    j j           j   j-1   j-1  j-2
!***
!***
!***                   v   = matrix of dimension iter*n
!***                         on entry first vector must be filled with
!***                         the starting vector to prime the Lanczos
!***                         recursion.
!***                   mat = the matrix you are trying to tridiagonalize
!***                   a,b = the recursion coefficients
!***                   scr = scratch vector of length n
!**                    vec = eigenvectors of tridiagonal matrix iter*iter
!***                   n   = vector length
!***                   iter = number of iterations
!***                   note as with all simple Lanczos methods, one can
!***                   get linear dependence problems and multiple copies
!***                   of eigenvalues if iter is large.
!***references

!***routines called
!***end prologue       lancz

  SUBROUTINE lancz(v,a,b,vec,hvec,n,iter)
  USE io
  IMPLICIT NONE (a-z)
  COMPLEX*16, DIMENSION(n,0:iter)        :: v
  REAL*8, DIMENSION(1,iter)              :: a, b
  REAL*8, DIMENSION(iter,iter)           :: vec
  COMPLEX*16, DIMENSION(n)               :: hvec
  INTEGER                                :: n, iter
  COMPLEX*16                             :: cdotc
  REAL*8                                 :: anorm
  CHARACTER (LEN=80) :: title
!
!         normalize the input vector
!
  anorm=SQRT(cdotc( n,v(1,0),1,v(1,0),1) )
  anorm=1.d0/anorm
  CALL cvscal(v(1,0),v(1,0),anorm,n)
!
!         we now have a normalized first function
!
  IF (iter > 1) THEN
!  
!         form mat times the first function
!  
     CALL h_v_dvr(v(1,0),hvec,n)
  
!         calculate a(1)
  
     a(1)=cdotc(n,v(1,0),1,hvec,1)
  
!         form mat times the first function - a(1) times the first function
!         and store it in the next polynomial
  
     v(:,1) = hvec(:) - a(1) * v(:,0)
  
!         calculate b(1)
  
     b(1)=SQRT( cdotc(n,v(1,1),1,v(1,1),1) )
  
!         normalize the second polynomial
  
     CALL cvscal(v(1,0),v(1,0),1.d0/b(1),n)
  END IF
  IF (iter > 2) THEN
      DO i=2,iter
    
!            multiply the last calculated polynomial by mat
    
         CALL h_v_dvr(v(1,i-1),hvec,n)
    
!            orthogonalize it to the two previous polynomials
!            calculating a(i) as we go
    
         a(i)=cdotc(n,v(1,i-1),1,hvec,1)
         v(:,i)=hvec(:) - a(i)*v(:,i-1)
         v(:,i) = v(:,i) - b(i-1)*v(:,i-2)
    
!            calculate b(i)
    
         b(i) = SQRT( cdotc(n,v(1,i),1,v(1,i),1) )
    
!            normalize the polynomial and we are done
    
        CALL cvscal(v(1,i),v(1,i),1.d0/b(i),n)
      END DO
  END IF

!     diagonalize the tridiagonal matrix for the eigenvalues

  vec = 0.d0
  DO  i=1,iter
    vec(i,i)=1.d0
  END DO
  CALL copy(b,scr(2),iter-1)
  CALL imtql2(iter,iter,a(1),scr,vec,ierr)
  title='eigenvalues of lanczos polynomial'
  CALL prntrm(title,a(1),iter,1,iter,1,iout)
END SUBROUTINE lancz

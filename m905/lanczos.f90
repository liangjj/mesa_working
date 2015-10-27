!deck lanczos.f
!***begin prologue     lanczos
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
!***end prologue       lanczos

  SUBROUTINE lanczos(it)
  USE lanczos_global
  IMPLICIT NONE 
  COMPLEX*16                             :: cdotc
  INTEGER                                :: it
!
!    Compute b(it) and the next normalized polynomial
!    
     b(it)=SQRT(cdotc( n,v(1,it-1),1,v(1,it-1),1) )
     CALL cvscal(v(1,it-1),v(1,it-1),1.d0/b(it),n)
!
!    Perform matrix multiply
!
     CALL lanczos_mult(v(1,it-1),hvec,n)
!    
!    Calculate a(it)
!    
     a(it)=cdotc(n,v(1,it-1),1,hvec,1)
!
!    Calculate the right hand side of the recursion
!
     v(:,it) = hvec(:) - a(it)*v(:,it-1) 
     IF(it > 1) then
        v(:,it) = v(:,it) - b(it)*v(:,it-2)
     END IF
END SUBROUTINE lanczos

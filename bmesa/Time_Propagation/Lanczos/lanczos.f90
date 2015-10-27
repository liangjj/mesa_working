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
!***                   v    = matrix of dimension maxit*n3d
!***                          on entry first vector must be filled with
!***                          the starting vector to prime the Lanczos
!***                          recursion.
!***                   hvec =  vector of dimension n3d
!***                   a,b  = the recursion coefficients of dimension maxit
!**                    vec  = eigenvectors of tridiagonal matrix 
!***                          of dimension maxit*maxit
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
     b_lanczos(it)=SQRT(cdotc( n3d,vec(1,it-1),1,vec(1,it-1),1) )
     CALL cvscal(vec(1,it-1),vec(1,it-1),1.d0/b_lanczos(it),n3d)
!
!    Perform matrix multiply
!
     CALL h_v_dvr(vec(1,it-1),hvec,1)
!    
!    Calculate a_lanczos(it)
!    
     a_lanczos(it)=cdotc(n3d,vec(1,it-1),1,hvec,1)
!
!    Calculate the right hand side of the recursion
!
     vec(:,it) = hvec(:) - a_lanczos(it)*vec(:,it-1) 
     IF(it > 1) then
        vec(:,it) = vec(:,it) - b_lanczos(it)*vec(:,it-2)
     END IF
END SUBROUTINE lanczos

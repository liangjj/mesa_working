!deck c_vect_2d.f
!**begin prologue     c_vect_2d
!**date written       040607   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             3-dim
!**purpose            Sort eigenvectors according energy and place
!***                  the desired root into initial state.
!**references
!**routines called
!**end prologue       c_vect_2d
  SUBROUTINE c_vect_2d(wave_function,etemp,itemp,root)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                 :: root
  REAL*8, DIMENSION(nphy(2),nphy(1))      :: wave_function
  REAL*8, DIMENSION(nphy(2)*nphy(1))      :: etemp
  INTEGER, DIMENSION(nphy(2)*nphy(1),2)   :: itemp
  REAL*8                                  :: tmp
  INTEGER                                 :: i, j, k, i1, j1
  INTEGER                                 :: ii, count
!
!       The eig and ind arrays are destroyed by this routine.
!
!
!       Set up the index array
!
  count=0 
  do i=1,nphy(1)
     do j=1,nphy(2)
        count=count+1
        itemp(count,1)=i
        itemp(count,2)=j
     END DO
  END DO
!
!        Fill the eigenvalue array with all of the eigenvalues of the zeroth
!        order Hamiltonian.  These are simply the sum of the eigenvalues of
!        each dimension
!
  count=0
  DO  i=1,nphy(1)
      do j=1,nphy(2)
         count = count + 1
         etemp(count) = grid(1)%eigv_0(i) + grid(2)%eigv_0(j)
      END DO
  END DO
!
!        Sort the eigenvalues so that we can find the smallest
!
  DO  ii=2,n3d
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      DO  j=ii,n3d
          IF(etemp(j) < tmp) THEN
             k=j
             tmp=etemp(j)
          END IF
      END DO
      IF(k /= i) THEN
         itemp(i,1) = itemp(k,1)
         itemp(i,2) = itemp(k,2)
         itemp(k,1)=i1
         itemp(k,2)=j1
         etemp(k) = etemp(i)
         etemp(i) = tmp
      END IF
  END DO
!
!        Now that we know which eigenvector in each dimension fill the psi0
!        arrray with the appropriate eigenvector as the product.
!
  DO  i=1,nphy(1)
      DO j=1,nphy(2)
         wave_function(j,i) =  grid(1)%eigvec_0(i,itemp(root,1))      &
                                             *                        &
                               grid(2)%eigvec_0(j,itemp(root,2))
      END DO
  END DO
  energy=etemp(root)
  WRITE(iout,1) root, energy
1    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE c_vect_2d

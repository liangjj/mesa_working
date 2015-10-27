! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Psi0_fil_2}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck psi0_fil_3
!**begin prologue     psi0_fil_3
!**date written       030308   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           fill, sort
!**author             schneider, barry (nsf)
!**source             arnoldi
!**purpose            sort and fill etempenvector.
!**references
!**routines called
!**end prologue       psi0_fil_3
  SUBROUTINE psi0_fil_3(etemp,itemp,root)
  USE arnoldi_global
  USE dvr_Shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                 :: etemp
  INTEGER, DIMENSION(n3d,3)              :: itemp
  INTEGER                                :: root
  REAL*8                                 :: tmp
  INTEGER                                :: i, j, k, i1, j1, k1
  INTEGER                                :: ii, count
!
! Set up the index array
!
  count=0
  do i=1,nphy(1)
     do j=1,nphy(2)
        do k=1,nphy(3)
           count=count+1
           itemp(count,1)=i
           itemp(count,2)=j
           itemp(count,3)=k
        END DO
     END DO
  END DO 
!
! Fill the etempenvalue array with all of the etempenvalues of the zeroth
! order Hamiltonian.  These are simply the sum of the etempenvalues of
! each dimension
!
  count=0
  DO  i=1,nphy(1)
      do j=1,nphy(2)
         do k=1,nphy(3)
            count = count + 1
            etemp(count) = grid(1)%eigv_0(i) + &
                           grid(2)%eigv_0(j) + &
                           grid(3)%eigv_0(k)
         end do
      END DO
  END DO
!
! Sort the etempenvalues so that we can fitemp the smallest
!
  DO  ii=2,n3d
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      k1=itemp(i,3)
      DO  j=ii,n3d
          IF(etemp(j) < tmp) THEN
             k=j
             tmp=etemp(j)
          END IF
      END DO
      IF(k /= i) THEN
         itemp(i,1)=itemp(k,1)
         itemp(i,2)=itemp(k,2)
         itemp(i,3)=itemp(k,3)
         itemp(k,1)=i1
         itemp(k,2)=j1
         itemp(k,3)=k1
         etemp(k) = etemp(i)
         etemp(i) = tmp
      END IF
  END DO
!
! Now that we know which etempenvector in each dimension fill the psi0
! arrray with the appropriate etempenvector as the product.
!
  count = 0
  DO  i=1,nphy(1)
      do j=1,nphy(2)
         do k=1,nphy(3)
            count = count + 1
            psi0(count) = grid(1)%eigvec_0(i,itemp(root,1)) * &
                          grid(2)%eigvec_0(j,itemp(root,2)) * &
                          grid(3)%eigvec_0(k,itemp(root,2))
         end do
      end do
  END DO
END SUBROUTINE psi0_fil_3






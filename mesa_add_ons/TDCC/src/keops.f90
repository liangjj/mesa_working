module keops
  USE nrtype
  !************************************************************************************
  !************************************************************************************
  ! ----------------------------------------------------------------------
contains
  ! ----------------------------------------------------------------------
  subroutine KE_quarterdt_3D_xodd
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    !************************************************************************************
    !************************************************************************************
    !  Kinetic-energy operator action on the wave function for half time-step (dt/2)

    USE globalmod
    USE pfedvrmod
    implicit NONE

    INTEGER(I4B) :: i, j, k, m, n
    REAL(DP) :: ctemp

    !-----------------------------------------------------------------------------------
    !   Looping  over the parallelizing dimension (1)
    !-----------------------------------------------------------------------------------

    !------------------------------------performing on the odd elements---------------------------------------------------

    DO i=ngridstart(2), ngridstop(2)       
       DO j=1, nwave
          DO k=regstart(1), regstop(1), 2

             IF( k.eq.1 )THEN

                do m=2, num_fun(1,k)
                   ctemp=czero
                   do n=2, num_fun(1,k)
                      ctemp=ctemp+Kin_operator_quarterdt(1,k,m,n)*coeff3d(nindex(1,k,n),i,j)
                   end do
                   ctemp_mat3(nindex(1,k,m))=ctemp
                end do

                do n=2, num_fun(1,k)
                   coeff3d(nindex(1,k,n),i,j)=ctemp_mat3(nindex(1,k,n)) 
                end do

             ELSE

                do m=1, num_fun(1,k)
                   ctemp=czero
                   do n=1, num_fun(1,k)
                      ctemp=ctemp+Kin_operator_quarterdt(1,k,m,n)*coeff3d(nindex(1,k,n),i,j)
                   end do
                   ctemp_mat3(nindex(1,k,m))=ctemp
                end do

                do n=1, num_fun(1,k)
                   coeff3d(nindex(1,k,n),i,j)=ctemp_mat3(nindex(1,k,n)) 
                end do

             ENDIF


          END DO
       END DO
    END DO

    !------------------------------------------------------------

    return

  end subroutine KE_quarterdt_3D_xodd




  !************************************************************************************
  !************************************************************************************
  ! ----------------------------------------------------------------------
  ! ----------------------------------------------------------------------
  subroutine KE_halfdt_3D_xeven
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    !************************************************************************************
    !************************************************************************************
    !  Kinetic-energy operator action on the wave function for half time-step (dt/2)

    USE globalmod
    USE pfedvrmod
    implicit NONE

    INTEGER(I4B) :: i, j, k, m, n
    REAL(DP) :: ctemp

    !-----------------------------------------------------------------------------------	
    !   Looping  over the parallelizing dimension (1)
    !-----------------------------------------------------------------------------------	

    !------------------------------------performing on the even elements---------------------------------------------------

    DO i=ngridstart(2), ngridstop(2)       
       DO j=1, nwave
          DO k=regstart(1)+1, regstop(1), 2

             do m=1, num_fun(1,k)
                ctemp=czero
                do n=1, num_fun(1,k)
                   ctemp=ctemp+Kin_operator_halfdt(1,k,m,n)*coeff3d(nindex(1,k,n),i,j)
                end do
                ctemp_mat3(nindex(1,k,m))=ctemp
             end do

             do n=1, num_fun(1,k)
                coeff3d(nindex(1,k,n),i,j)=ctemp_mat3(nindex(1,k,n)) 
             end do

          END DO
       END DO
    END DO

    !------------------------------------------------------------

    return

  end subroutine KE_halfdt_3D_xeven


  !************************************************************************************
  !************************************************************************************
  ! ----------------------------------------------------------------------
  ! ----------------------------------------------------------------------
  subroutine KE_quarterdt_3D_yodd
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    !************************************************************************************
    !************************************************************************************
    !  Kinetic-energy operator action on the wave function for half time-step (dt/2)

    USE globalmod
    USE pfedvrmod
    implicit NONE

    INTEGER(I4B) :: i, j, k, m, n
    REAL(DP) :: ctemp

    !-----------------------------------------------------------------------------------	
    !   Looping  over the parallelizing dimension (2)
    !-----------------------------------------------------------------------------------	

    !------------------------------------performing on the odd elements---------------------------------------------------

    DO i=ngridstart(1), ngridstop(1)       
       DO j=1, nwave
          DO k=regstart(2), regstop(2), 2

             IF( k.eq.1 )THEN

                do m=2, num_fun(2,k)
                   ctemp=czero
                   do n=2, num_fun(2,k)
                      ctemp=ctemp+Kin_operator_quarterdt(2,k,m,n)*coeff3d(i,nindex(2,k,n),j)
                   end do
                   ctemp_mat3(nindex(2,k,m))=ctemp
                end do

                do n=2, num_fun(2,k)
                   coeff3d(i,nindex(2,k,n),j)=ctemp_mat3(nindex(2,k,n)) 
                end do

             ELSE

                do m=1, num_fun(2,k)
                   ctemp=czero
                   do n=1, num_fun(2,k)
                      ctemp=ctemp+Kin_operator_quarterdt(2,k,m,n)*coeff3d(i,nindex(2,k,n),j)
                   end do
                   ctemp_mat3(nindex(2,k,m))=ctemp
                end do

                do n=1, num_fun(2,k)
                   coeff3d(i,nindex(1,k,n),j)=ctemp_mat3(nindex(2,k,n)) 
                end do

             ENDIF


          END DO
       END DO
    END DO

    !------------------------------------------------------------

    return

  end subroutine KE_quarterdt_3D_yodd




  !************************************************************************************
  !************************************************************************************
  ! ----------------------------------------------------------------------
  ! ----------------------------------------------------------------------
  subroutine KE_halfdt_3D_yeven
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    !************************************************************************************
    !************************************************************************************
    !  Kinetic-energy operator action on the wave function for half time-step (dt/2)

    USE globalmod
    USE pfedvrmod
    implicit NONE

    INTEGER(I4B) :: i, j, k, m, n
    REAL(DP) :: ctemp

    !-----------------------------------------------------------------------------------	
    !   Looping  over the parallelizing dimension (1)
    !-----------------------------------------------------------------------------------	

    !------------------------------------performing on the even elements---------------------------------------------------

    DO i=ngridstart(1), ngridstop(1)       
       DO j=1, nwave
          DO k=regstart(2)+1, regstop(2), 2

             do m=1, num_fun(2,k)
                ctemp=czero
                do n=1, num_fun(2,k)
                   ctemp=ctemp+Kin_operator_halfdt(2,k,m,n)*coeff3d(i,nindex(2,k,n),j)
                end do
                ctemp_mat3(nindex(2,k,m))=ctemp
             end do

             do n=1, num_fun(2,k)
                coeff3d(i,nindex(2,k,n),j)=ctemp_mat3(nindex(2,k,n)) 
             end do

          END DO
       END DO
    END DO

    !------------------------------------------------------------

    return

  end subroutine KE_halfdt_3D_yeven

end module keops

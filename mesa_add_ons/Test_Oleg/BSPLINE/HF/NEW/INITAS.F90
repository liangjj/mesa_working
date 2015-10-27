!==================================================================
  SUBROUTINE initas
!==================================================================
!
!   Sets ( or Initializes ) the array in symmetric storage mode:
!
!       db1 --- matrix of integral <B_i,B'_j>
!       db2 --- matrix of integral <B_i,B"_j>
!       sb  --- matrix of integral <B_i,B_j>
!       r1  --- matrix of integral <B_i,r B_j>
!       rm1 --- matrix of integral <B_i,(1/r)B_j>
!       rm2 --- matrix of integral <B_i,(1/r^2)B_j>
!               where i=1,..,ns, j=1,...ks
!
!   SUBROUTINES called
!       mdb mrm
!
!   Calling sequence:
!       initas
!        /  \
!      mdb  mrm
!
!------------------------------------------------------------------
!
    USE spline_galerkin

    IMPLICIT NONE

    ! .. sets db1 --- matrix of integral <B_i,B'_j>
    CALL mdb(1, db1)  
 
    ! .. sets db2 --- matrix of integral <B_i,B"_j>
    CALL mdb(2, db2)  
 
    ! .. sets sb  --- matrix of integral <B_i,B_j>
    CALL mrm(0, sb)  

    ! .. sets r1  --- matrix of integral <B_i,r B_j>
    CALL mrm(1, r1)  

    ! .. sets rm1 --- matrix of integral <B_i,(1/r)B_j>
    CALL mrm(-1, rm1)  
 
    ! .. sets rm2 --- matrix of integral <B_i,(1/r^2)B_j>
    CALL mrm(-2, rm2)  

  END SUBROUTINE initas

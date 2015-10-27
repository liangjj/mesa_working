!==================================================================== 
      Real*8 FUNCTION F_3J (J1,M1,J2,M2,J,M,clebsch)
!==================================================================== 
! 
!     determines the 3J or Clebsh-Gordon coefficients 
!     (the momenta are used in (2J+1)-representation) 
! 
!     Call:  F_3J 
! 
!     determines the value of the 3j-symbols without direct using of 
!     factorials. The following expression for the 3j-symbols is used: 
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977) 
! 
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} * 
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ] 
!                         SUM(z) {   (-1)^z  / 
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! * 
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] } 
! 
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ] 
! 
!     If we introduce the auxiliary values a(i) and b(i) 
!     (see below the text of program) then 
! 
!     3j =         (-1) ^ Sum[a(i)] 
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] } 
!                  Sum(z) { (-1)^z  / 
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] } 
! 
!     (below the moments are used in (2J+1)-representation) 
! 
!-------------------------------------------------------------------- 
!-------------------------------------------------------------------- 
      Implicit None 
      INTEGER, intent(in) :: J1, M1, J2, M2, J, M 
      INTEGER             :: J_a, M_a, J_b, M_b, J_c, M_c 
      LOGICAL             :: clebsch 
      Real*8, EXTERNAL    :: Z_3j 
      J_a = J1 + J1 + 1
      M_a = M1 + M1 + 1
      J_b = J2 + J2 + 1
      M_b = M2 + M2 + 1
      J_c = J + J + 1
      M_c = M + M + 1
      IF (clebsch) THEN
          F_3J = (-1)**((J_a-J_b+M_c-1)/2)*sqrt(DBLE(J_c))*   & 
                  Z_3J(J_a,M_a,J_b,M_b,J_c,-M_c+2) 
      ELSE
          F_3J = Z_3J(J_a,M_a,J_b,M_b,J_c,M_c)
      END IF
      END FUNCTION F_3J
!-------------------------------------------------------------------- 
!-------------------------------------------------------------------- 

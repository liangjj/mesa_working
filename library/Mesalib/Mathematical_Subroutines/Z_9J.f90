!======================================================================
      FUNCTION Z_9J (j1,j2,j3,j4,j5,j6,j7,j8,j9)
!======================================================================
!
!     determination of 9j-symbols     
!                                    
! {j1 j2 j3}                           {j1 j2 j3} {j4 j5 j6} {j7 j8 j9}        
! {j4 j5 j6} = SUM(j) (-1)^(2j) (2j+1) 
! {j7 j8 j9}                           {j6 j9 j } {j2 j  j8} {j  j1 j4}
!
!
!     (the momenta are used in (2J+1)-representation)
!
!----------------------------------------------------------------------

      Implicit none

      Integer*4, intent(in)                :: j1,j2,j3,j4,j5,j6,j7,j8,j9
      Integer*4                            :: j,i1,i2
      Real*8, External                     :: Z_6J
      Real*8                               :: Z_9J
      
      i1 = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))+1
	  i2 = min(j1+j9,j2+j6,j4+j8)-1

	  Z_9J = 0.d0; if(i1.gt.i2) Return;  if(mod(i2-i1,2).ne.0) Return

      Do j = i1,i2,2
       Z_9J = Z_9J + j * (-1)**(j-1) * &
                     Z_6J(j1,j2,j3,j6,j9,j ) * &
                     Z_6J(j4,j5,j6,j2,j ,j8) * &
                     Z_6J(j7,j8,j9,j ,j1,j4) 
      End do
       
      End FUNCTION Z_9J

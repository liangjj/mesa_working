!====================================================================
     FUNCTION Z_6JJ(j1,j2,j3,j4,j5,j6)
!====================================================================

      IMPLICIT NONE

      Integer*4, intent(in)             :: j1,j2,j3,j4,j5,j6
      Real*8                            :: Z_6JJ
      Real*8, External                  :: Z_6J

      Z_6JJ = Z_6J(j1+j1+1,j2+j2+1,j3+j3+1,j4+j4+1,j5+j5+1,j6+j6+1)

      End FUNCTION Z_6JJ

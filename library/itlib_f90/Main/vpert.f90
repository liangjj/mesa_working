!deck vpert.f
!***begin prologue     vpert
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            potential matrix
!***
!***description        calculate the dependent potential
!***                   matrix elements in the dvr representation.
!***
!***references

!***routines called
!***end prologue       vpert
  SUBROUTINE vpert(key,vword,zeroit,prnt)
  USE io
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  CHARACTER (LEN=*)                      :: key
  INTEGER                                :: vword
  LOGICAL                                :: zeroit, prnt
  REAL*8, DIMENSION(3,3)                 :: a, b, d, len, dif
  REAL*8                                 :: fpkey, depth
  CHARACTER (LEN=32)                     :: TYPE
  CHARACTER (LEN=80)                     :: chrkey, title
  CHARACTER (LEN=1)                      :: itoc
  CHARACTER (LEN=24), DIMENSION(4)       :: phr
  LOGICAL                                :: dollar, logkey
  INTEGER                                :: i, j
  ALLOCATE (v_pert(nphy(4)))
  IF(zeroit) THEN
     v_pert=0.d0
  END IF
  IF(spdim == 1) THEN
     RETURN
  END IF
  IF(dollar(key,card,cpass,inp) ) THEN
     TYPE=chrkey(card,'potential','none',' ')
     prnt=logkey(card,'print=potential',.false.,' ')
     WRITE(iout,1) TYPE
  END IF
  IF(TYPE == 'none') THEN
     RETURN
  ELSE IF(TYPE == 'exponential') THEN
     DO  i=1,spdim
         phr(1)='a'//itoc(i)
         phr(2)='b'//itoc(i)
         DO  j=1,spdim
             phr(3)=phr(1)(1:2)//itoc(j)
             phr(4)=phr(2)(1:2)//itoc(j)
             a(i,j)=fpkey(card,phr(3)(1:3),0.d0,' ')
             b(i,j)=fpkey(card,phr(4)(1:3),0.d0,' ')
         END DO
     END DO
     DO  i=1,spdim
         DO  j=1,i
             IF(a(i,j) /= a(j,i).OR.b(i,j) /= b(j,i)) THEN
                CALL lnkerr('error in off-diagonal exponential '//  &
                            'parameters')
             END IF
         END DO
     END DO
     IF(spdim == 2) THEN
        CALL vexp12(v_pert,grid(1)%pt,grid(2)%pt,a(1,2),b(1,2),  &
                    nphy(1),nphy(2),.false.)
        IF(prnt) THEN
           title='2d exponential'
        END IF
     END IF
     IF(spdim == 3) THEN
        CALL vexp123(v_pert,grid(1)%pt,grid(2)%pt,grid(3)%pt,     &
                     a(1,2),a(1,3),a(2,3), b(1,2),b(1,3),b(2,3),  &
                     nphy(1),nphy(2),nphy(3),.false.)
     END IF
  ELSE
     CALL lnkerr('error in potential type')
  END IF
1    FORMAT(/,5X,'interaction potential = ',a32)
END SUBROUTINE vpert















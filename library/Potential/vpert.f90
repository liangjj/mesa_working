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
  SUBROUTINE vpert(x1,x2,x3,v1,v2,v3,v4,n,dim,key,vword,zeroit,prnt)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: dim
  INTEGER, DIMENSION(dim+1)              :: n 
  REAL*8, DIMENSION(n(1))                :: x1
  REAL*8, DIMENSION(n(2))                :: x2
  REAL*8, DIMENSION(n(3))                :: x3
  REAL*8, DIMENSION(n(1))                :: v1
  REAL*8, DIMENSION(n(2))                :: v2
  REAL*8, DIMENSION(n(3))                :: v3
  REAL*8, DIMENSION(n(4))                :: v4
  CHARACTER (LEN=*)                      :: key
  INTEGER                                :: vword
  LOGICAL                                :: zeroit, prnt
  REAL*8, DIMENSION(3,3)                 :: a, b, d, len, dif
  REAL*8                                 :: fpkey, depth
  CHARACTER (LEN=32)                     :: TYPE
  CHARACTER (LEN=1600)                   :: card
  CHARACTER (LEN=80)                     :: cpass, chrkey, title
  CHARACTER (LEN=1)                      :: itoc
  CHARACTER (LEN=24), DIMENSION(4)       :: phr
  LOGICAL                                :: dollar, logkey
  INTEGER                                :: i, j
  IF(zeroit) THEN
     v4=0.d0
  END IF
  IF(dim == 1) THEN
     RETURN
  END IF
  IF(dollar(key,card,cpass,inp) ) THEN
     TYPE=chrkey(card,'potential','none',' ')
     prnt=logkey(card,'print=potential',.false.,' ')
     WRITE(iout,1) TYPE
  END IF
  IF(TYPE == 'none') THEN
     RETURN
  ELSE IF(TYPE == 'well') THEN
     IF(dim >= 2) THEN
        dif(1,2)=ABS( x2(n(2)) - x1(n(1)) )
        dif(2,1)=dif(1,2)
     END IF
     IF(dim == 3) THEN
        dif(1,3)=ABS( x1(n(1)) - x3(n(3)) )
        dif(3,1)=dif(1,3)
        dif(2,3)=ABS( x2(n(2)) - x3(n(3)) )
        dif(3,2)=dif(2,3)
     END IF
     DO  i=1,dim
         phr(1)='d'//itoc(i)
         phr(2)='l'//itoc(i)
         DO  j=1,dim
             phr(3)=phr(1)(1:2)//itoc(j)
             phr(4)=phr(2)(1:2)//itoc(j)
             d(i,j)=fpkey(card,phr(3)(1:3),0.d0,' ')
             len(i,j)=fpkey(card,phr(4)(1:3),dif(i,j),' ')
         END DO
     END DO
     DO  i=1,dim
         DO  j=1,i
             IF(d(i,j) /= d(j,i).OR.LEN(i,j) /= LEN(j,i)) THEN
                CALL lnkerr('error in off_diagonal well parameters')
             END IF
         END DO
     END DO
     IF(dim == 2) THEN
        CALL vwel12(v4,x1,x2,LEN(2,1),d(2,1),n(1),n(2),.false.)
     END IF
     IF(dim == 3) THEN
        CALL vwel123(v4,x1,x2,x3,LEN(2,1),LEN(3,1),LEN(3,2),  &
                     d(2,1),d(3,1),d(3,2),n(1),n(2),n(3),.false.)
     END IF
     IF(prnt) THEN
        title='well'
        CALL prntrm(title,v4,n(dim+1),1,n(dim+1),1,iout)
     END IF
  ELSE IF(TYPE == 'exponential') THEN
     DO  i=1,dim
         phr(1)='a'//itoc(i)
         phr(2)='b'//itoc(i)
         DO  j=1,dim
             phr(3)=phr(1)(1:2)//itoc(j)
             phr(4)=phr(2)(1:2)//itoc(j)
             a(i,j)=fpkey(card,phr(3)(1:3),0.d0,' ')
             b(i,j)=fpkey(card,phr(4)(1:3),0.d0,' ')
         END DO
     END DO
     DO  i=1,dim
         DO  j=1,i
             IF(a(i,j) /= a(j,i).OR.b(i,j) /= b(j,i)) THEN
                CALL lnkerr('error in off_diagonal exponential '// 'parameters')
             END IF
         END DO
     END DO
     IF(dim == 2) THEN
        CALL vexp12(v4,x1,x2,a(1,2),b(1,2),n(1),n(2),.false.)
        IF(prnt) THEN
           title='2d exponential'
        END IF
     END IF
     IF(dim == 3) THEN
        CALL vexp123(v4,x1,x2,x3,a(1,2),a(1,3),a(2,3), b(1,2),b(1,3),b(2,3),  &
                     n(1),n(2),n(3),.false.)
     END IF
  ELSE
     CALL lnkerr('error in potential type')
  END IF
1    FORMAT(/,5X,'interaction potential = ',a32)
END SUBROUTINE vpert















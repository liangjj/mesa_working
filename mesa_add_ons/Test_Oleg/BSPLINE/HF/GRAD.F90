!=======================================================================
  FUNCTION grad(i,j,ier) 
!=======================================================================
!  The grad function computes the following directly
!          <P(J)^D + L(I)/R ^P(I)> WITH L(I) > L(J)
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j
    INTEGER, INTENT(OUT) :: ier
    REAL(KIND=8) :: grad

    LOGICAL :: ixch
    REAL(KIND=8), DIMENSION(ns,ks) :: g, a
      
    IF ( IABS(L(I) - L(J)) .NE. 1) THEN
      ier = 1
      grad = 0
    else
      LL = MAX0(L(I),L(J))
      i1 = i 
      i2 = j
      ixch = .false.
      IF ( L(I) < L(J) ) THEN
	i1 = j
	i2 = i
	ixch = .true.
      end if

      !  form the sums

      last = min0(max(i1),max(i2))
      grad = SUM( p(2:last,i1)*p(2:last,i2)*rm1(2:last,ks)
      grad = ll*grad
      do jp = 2,ks
        ! .. j' < jj
	lastjp = last-jp+1
        grad = grad  &       
	    +  SUM( g(2:lastjp,ks-jp+1)*                                     &
                          (p(2:lastjp,i2)*p(2:lastjp-jp+1,i1) -              &
	                   p(2:lastjp,i1)*p(2:lastjp-jp+1,i2)) )             &
	    + ll*SUM( rm1(2:lastjp,ks-jp+1)*                                 &
			  (p(2:lastjp,i2)*p(2:lastjp-jp+1,i1)+               &
                           p(2:lastjp,i1)*p(2:lastjp-jp+1,i2)) )
      end do
      IF ( ixch ) grad = -grad
    end if
    END FUNCTION grad

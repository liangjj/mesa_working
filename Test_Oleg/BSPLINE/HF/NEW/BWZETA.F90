!=======================================================================
   FUNCTION bwzeta(i1)
!=======================================================================
!
!  Computes the nuclear spin-orbit parameter and the
!  using the formula derived by blume and watson.
!
!
!----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i1
     REAL(KIND=8) :: bwzeta

      dimension ss(3)
!
!
      zeta = fine*z*quadr(i1,i1,-3)
      lb = l(i1)
      do i = 1,nwf
        if (i .eq. i1) cycle
        la = l(i)
        zeta = zeta -sum(i)*sn(i1, i, i1, i, 0)
        if (sum(i) .ne. 4*l(i)+2) go to 10
          call bwint(la,lb)
          ke1 = 2
          if (la .ne. lb) ke1 = iabs(la-lb)
          ip = 0
          do k = ke1,la+lb,2
            ip = ip+1
            zeta = zeta+coefn2(ip)*sn(i1, i, i, i1, k-2)
     :                 +coefnk(ip)*sn(i, i1, i1, i, k)
     :                 +coefvk(ip)*(vk(i1,i,i,i1,k-1)-vk(i,i1,i1,i,k-1))
          end do
      end do
      zeta = 2.d0*zeta
      c= sum(i1)
      if (c .ne. 1.d0) then
         ss(1) = sn(i1,i1,i1,i1,0)
         c = c + c - 3.d0
         zeta = zeta - c*ss(1)
         if (lb .eq. 2) then
            ss(2) = sn(i1,i1,i1,i1,2)
            zeta = zeta + ss(2)*6.d0/7.d0
         else if (lb .eq. 3) then
            ss(2) = sn(i1,i1,i1,i1,2)
            ss(3) = sn(i1,i1,i1,i1,4)
            zeta = zeta + ss(2) + ss(3)/2.2d0
         end if
      end if
      bwzeta = zeta
 
  CONTAINS
  !=====================================================================
    SUBROUTINE bwint(lc,lo)
  !=====================================================================
  !
  !   Compute the spin-orbit parameter using the Blume and Watson method
  !
  !----------------------------------------------------------------------
  !
  
    IMPLICIT NONE
      COMMON/blume/coefn2(4),coefnk(4),coefvk(4)
  !
  ! ... lc is the l-value of the filled subshell, lo is the l-value
  !     of the partially-filled subshell.
  !
      if(lc > 3 .or. lo >  4) then
      write(iscw,'(A,I3,A,I3)') 'Incorrect calling of bwint with lc =', &`j
         lc, 'and lo =',lo
         return
      end if
      lc1 = lc + 1

      go to (10,20,30,40), lc1
   10 go to (11,12,13,14), lo
  !
  ! ... s-p
  !
   11 coefnk(1) = 1.d0
      coefn2(1) = -2.d0
      coefvk(1) = 1.d0
      return
  !
  ! ... s-d
  !
   12 coefnk(1) = 6.d0/5.d0
      coefn2(1) = -9.d0/5.d0
      coefvk(1) = 3.d0/5.d0
      return
  !
  ! ... s-f
  !
   13 coefnk(1) = 9.d0/7.d0
      coefn2(1) = -12.d0/7.d0
      coefvk(1) = 3.d0/7.d0
      return
  !
  ! ... s-g
  !
   14 coefnk(1) = 4.d0/3.d0
      coefn2(1) = -5.d0/3.d0
      coefvk(1) = 1.d0/3.d0
      return
   20 go to (21,22,23,24), lo
  !
  ! ... p-p
  !
   21 coefnk(1) = 0.d0
      coefn2(1) = 3.d0
      coefvk(1) = 9.d0/5.d0
      return
  !
  ! ... p-d
  !
   22 coefnk(1) = 3.d0/7.d0
      coefnk(2) = 36.d0/35.d0
      coefn2(1) = -12.d0/5.d0
      coefn2(2) = 0.d0
      coefvk(1) = 3.d0/5.d0
      coefvk(2) = 36.d0/35.d0
      return
  !
  ! ... p-f
  !
   23 coefnk(1) = 1.d0/7.d0
      coefnk(2) = 10.d0/7.d0
      coefn2(1) = -18.d0/7.d0
      coefn2(2) = 0.d0
      coefvk(1) = 18.d0/35.d0
      coefvk(2) = 5.d0/7.d0
      return
  !
  ! ... p-g
  !
   24 coefnk(1) = 5.d0/77.d0
      coefnk(2) = 18.d0/11.d0
      coefn2(1) = -18.d0/7.d0
      coefn2(2) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 6.d0/11.d0
      return
   30 go to (31,32,33,34), lo
  !
  ! ... d-p
  !
   31 coefnk(1) = 59.d0/7.d0
      coefnk(2) = -18.d0/7.d0
      coefn2(1) = -4.d0
      coefn2(2) = 0.d0
      coefvk(1) = -1.d0
      coefvk(2) = 18.d0/7.d0
      return
  !
  ! ... d-d
  !
   32 coefnk(1) = 6.d0/7.d0
      coefnk(2) = 0.d0
      coefn2(1) = 3.d0
      coefn2(2) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 10.d0/7.d0
      return
  !
  ! ... d-f
  !
   33 coefnk(1) = 9.d0/7.d0
      coefnk(2) = -13.d0/77.d0
      coefnk(3) = 75.d0/77.d0
      coefn2(1) = -18.d0/7.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 3.d0/7.d0
      coefvk(3) = 75.d0/77.d0
      return
  !
  ! ... d-g
  !
   34 coefnk(1) = 741.d0/693.d0
      coefnk(2) = -215.d0/429.d0
      coefnk(3) = 210.d0/143.d0
      coefn2(1) = -3.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 255.d0/693.d0
      coefvk(3) = 105.d0/143.d0
      return
   40 go to (41,42,43,44), lo
  !
  ! ... f-p
  !
   41 coefnk(1) = 52.d0/3.d0
      coefnk(2) = -20.d0/3.d0
      coefn2(1) = -9.d0
      coefn2(2) = 0.d0
      coefvk(1) = -9.d0/5.d0
      coefvk(2) = 10.d0/3.d0
      return
  !
  ! ... f-d
  !
   42 coefnk(1) = 5.d0
      coefnk(2) = 142.d0/55.d0
      coefnk(3) = -20.d0/11.d0
      coefn2(1) = -18.d0/5.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = -3.d0/5.d0
      coefvk(2) = 2.d0/5.d0
      coefvk(3) = 20.d0/11.d0
      return
  !
  ! ... f-f
  !
   43 coefnk(1) = 1.d0
      coefnk(2) = 5.d0/11.d0
      coefnk(3) = 0.d0
      coefn2(1) = 3.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = 1.d0/5.d0
      coefvk(2) = 5.d0/11.d0
      coefvk(3) = 175.d0/143.d0
      return
  !
  ! ... f-g
  !
   44 coefnk(1) = 53.d0/33.d0
      coefnk(2) = 57.d0/143.d0
      coefnk(3) = -115.d0/429.d0
      coefnk(4) = 392.d0/429.d0
      coefn2(1) = -8.d0/3.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefn2(4) = 0.d0
      coefvk(1) = 1.d0/3.d0
      coefvk(2) = 3.d0/11.d0
      coefvk(3) = 57.d0/143.d0
      coefvk(4) = 392.d0/429.d0
      return

    END FUNCTION bwint
  
  END SUBROUTINE bwzeta

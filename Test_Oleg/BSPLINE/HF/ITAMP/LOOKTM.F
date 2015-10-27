.main

*     ------------------------------------------------------------------
*                       l o o k - t m
*     ------------------------------------------------------------------
*
*     add the deviations to the average energy for a partially filled
*       p- or d- shell
*
      SUBROUTINE looktm( l, sl, sen, q, ip, nsl)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
 
      CHARACTER sl*2, sen*1
      INTEGER iptr(5)
      CHARACTER*3 terms(51)
      data    iptr/6,11,19,35,51/
      data        terms/'3p2','1d2','1s0','4s3','2d3','2p1',
*             .. d2 and d3 terms
     :      '3f2','3p2','1g2','1d2','1s0','4f3','4p3','2h3','2g3',
     :      '2f3','2d1','2d3','2p3',
*            ... d4 terms ...
     :       '5d4','3h4','3g4','3f2','3f4','3d4','3p2','3p4',
     :       '1i4','1g2','1g4','1f4','1d2','1d4','1s0','1s4',
*            ... d5 terms ...
     :       '6s5','4g5','4f3','4d5','4p3','2i5','2h3','2g3',
     :       '2g5','2f3','2f5','2d1','2d3','2d5','2p3','2s5'/
 
*
*  --- search for a partially unfilled p- or d-shell
*
      n = q
      if (n .gt. 2*l+1) n = 4*l+2 - n
      ip = 0
      nsl = 0
      if (n .gt. 1  .and. l .le. 2) then
         if (l .eq. 1) then
            ibegin = 1
            iend = 6
         else
            ibegin = iptr(n-1) + 1
            iend = iptr(n)
         end if
1        i = ibegin
10       if (sl .eq. terms(i)(1:2)) then
            if (sen .eq. ' ' .or. sen .eq. terms(i)(3:3)) then
	       nsl = nsl + 1
	       ip = i
            end if
         end if
         i = i+1
         if (i .le. iend) go to 10
      else if ( n .eq. 1 .and. sl(1:1) .eq. '2') then
	 nsl = 1
      end if
      return
      end

*
*     ------------------------------------------------------------------
*                       d e v
*     ------------------------------------------------------------------
*
*     add the deviations to the average energy for a partially filled
*       p- or d- shell
*
      SUBROUTINE dev(iel, l, q, i, done)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
     
      LOGICAL done
      INTEGER f2pp(6), f2dd(45), f4dd(45)
      data    f2pp/-3,3,12,-9,0,6/
      data    f2dd/-58,77,50,-13,140,
*             ... d3 coefficients
     :       -93,42,-12,-57,123,105,69,-12,
*             ... d4 coefficients
     :        -105,-69,-24,66,12,39,21,57,
     :        -51,30,48,84,219,111,210,138,
*             ... d5 coefficients
     :        -175,-85,23,-22,-112,-76,-58,167,
     :        23,-85,59,140,104,86,320,113/
      data    f4dd/5,-70,15,50,140,
*             ... d3 coefficients
     :        -30,-105,30,55,-45,105,-15,30,
*             ... d4 coefficients
     :        -105,15,-10,45,-30,-45,70,-55,
     :        75,135,20, 0,30,-15,210,-30,
*             ... d5 coefficients
     :        -175,-50,-40,-85,35,50,110,-15,
     :        -5,125,-25,140,20,-40,-100,-55/
     
      done = .true.
      n = q
      if (n .gt. 2*l+1) n = 4*l+2 - n
      if (n .gt. 1) then
         if (l .eq. 1) then
            call add(2*f2pp(i)/25.d0,2,iel,iel,.true.)
         else if (l .eq. 2) then
            i = i-6
            call add(2*f2dd(i)/441.d0,2,iel,iel,.true.)
            call add(2*f4dd(i)/441.d0,4,iel,iel,.true.)
	 else
	    done = .false.
         end if
      end if
      return
      end

*
*     ------------------------------------------------------------------
*               n o d e c
*     ------------------------------------------------------------------
*
*      counts the number of nodes of the function pde(j) in the range
*   j = 40,...,m-10.   the node counting procedure counts the local max
*   and min values.   only nodes between sufficiently large max and
*   min values are counted.
*
*
      INTEGER FUNCTION nodec(pde,ns,m)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION pde(*)
*
*   ***** find max |pde(j)|
*
      m = idamax(ns,pde,1)
      pdem = pde(m)
*
*   *****  count the number of local max or min's
*
      ncc = 0
      sign = 0.d0
      diff1 = pde(4) - pde(3)
      do 2 j = 4,ns-1
        diff2 = pde(j+1) - pde(j)
        if (diff2*diff1 .lt. 0.d0 .and. diff1 .ne. 0.d0) then
*
*   *****  a max or min has been found.   test if it is
*          sufficiently large
*
          if ( dabs(pde(j))/pdem .gt. 0.05d0 ) then
*
            if (sign .eq. 0.d0 ) then
*              .. this is the first significant maximum
               m = j
	    else
*           .. this is not the first, test whether a sign change has
*              occurred since the last significant max or min
              if (pde(j)*sign .lt. 0.d0 ) ncc = ncc + 1
	    end if
*           .. reset for the next anti-node
            sign = pde(j)
          end if
	end if
        diff1 = diff2
2     continue
      nodec = ncc
      end

*deck thtw
	subroutine thtw(p,q,e,e2,d,x1,s)	! sturm sequence
	integer p ,q ,s
	real*8 e(*) ,e2(*) ,d(*) ,x1 ,u ,v ,machep
	data machep/1.0d0/
	if (machep .eq. 1.0d0) then
	   do while (1.0d0 + machep .gt. 1.0d0)
	      machep = 0.5d0*machep
	   end do
           machep = 2.0d0*machep
	end if
	s = p - 1
        u = 1.0d0
	do i = p ,q
	   if ( u .ne. 0.0d0 ) then
	      v = e2(i) / u
	   else 
	      if ( e2(i) .eq. 0.0d0 ) then
		 v = 0.0d0
	      else
		 v = abs(e(i)) / machep
	      end if
	   end if
           u = d(i) - x1 - v
           if (u .lt. 0.0d0) s = s + 1
	end do
	return
	end	

*deck @(#)prtrte.f	5.1  11/6/94
      subroutine prtrte(route,iout)
c***begin prologue     prtrte.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           print, route
c***author             martin, richard (lanl)
c***source             @(#)prtrte.f	5.1   11/6/94
c***purpose            prints the nonstandard route.
c***description
c     call prtrte(route,iout)
c       route   the nonstandard route.
c       iout    the output file.
c
c***references
c***routines called    none
c***end prologue       prtrte.f
      implicit none
c     --- input variables -----
      integer iout
c     --- input arrays (unmodified) ---
      character*(*) route
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer rtbeg,pos,j
c
 1000 format(' route:')
 1010 format(5x,128a1,(/,(8x,121a1)))
c
c     --- write the route in its nonstandard form.
c         the end of a route segment is delimited by a semicolon.
      write(iout,1000)
      rtbeg=1
   10    pos=index(route(rtbeg:),';')
         if(pos.eq.0) then
         else
            write(iout,1010) (route(j:j),j=rtbeg,rtbeg+pos-1)
            rtbeg=rtbeg+pos
            goto 10
         endif
c
c
      return
      end

*deck %W%  %G%
      subroutine chain(next)
c***begin prologue     chain
c***date written       870207  (yymmdd)
c***revision date      870708  (yymmdd)
c
c     8 july 1987      pws at lanl
c         changing to a cray-ctss version. the name now has a 'c'
c         appended.
c
c***keywords           chain, links, exit
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c***purpose            chain invkes the next link for mesa
c                      on the ctss cray's, this is done by writing the name
c                      of the next link to the controller (terminal)
c                      the ccl shell takes care of the rest.
c***description
c
c  n.b.   this routine cannot call either lnkerr or write to 'iout'.
c         the output file has already been closed at this point.
c
c***references
c***routines called
c***end prologue       chain
c
      character*(*) next,temp*8
      logical ioinq
      logical exist
c
c     ----- check for l998 which indicates that the job is over -----
c
      if (next.eq.'m998') stop
c
      temp=next
      call locase(temp,temp)
      exist=ioinq(temp,'exist')
      do 53 i=1,len(temp)
         if (temp(i:i).eq.' ') temp(i:i)=char(0)
   53 continue
      if (.not.exist) then
        call access(error,'dn'l,temp)
      endif
      call lgo(temp)
c
c
      return
      end

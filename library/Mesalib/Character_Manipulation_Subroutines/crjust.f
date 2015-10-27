*deck @(#)crjust.f	5.1  11/6/94
      subroutine crjust(instr,outstr)
c***begin prologue     crjust
c***date written       850601  (yymmdd)
c***revision date      871025  (yymmdd)
c
c   25 october 1987    rlm at lanl
c      fixing case of blank string by 'returning'.
c
c***keywords           character, right justify
c***author             martin, richard (lanl)
c***source             @(#)crjust.f	5.1   11/6/94
c***purpose            right justifies one string into another.
c***description
c                      call crjust(instr,outstr)
c                        instr   input character string.
c                        outstr  output character string.
c
c                      the two strings may be the same.
c***references
c***routines called    (none)
c***end prologue       crjust
      implicit integer(a-z)
      character*(*) instr,outstr
c
c
      lenin=len(instr)
      lenout=len(outstr)
c     find the length of instr excluding trailing blanks
      li=cskipb(instr,' ')
      if(li.eq.0) then
c
c        return if a blank string was sent.
c
      else if(li.eq.lenout) then
         outstr=instr
      else if(li.lt.lenout) then
         outstr(lenout+1-li:lenout)=instr(1:li)
         outstr(1:lenout-li)=' '
      else
         outstr=instr(li+1-lenout:li)
      endif
c
c
      return
      end

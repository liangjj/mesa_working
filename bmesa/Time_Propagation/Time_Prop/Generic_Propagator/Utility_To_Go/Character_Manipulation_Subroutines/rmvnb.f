*deck @(#)rmvnb.f	5.1  11/6/94
      subroutine rmvnb(instr,outstr)
c***begin prologue     rmvnb
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           nulls, blanks
c***author             martin, richard (lanl)
c***source             @(#)rmvnb.f	5.1   11/6/94
c***purpose            removes nulls and blanks from a character string.
c***description
c     call rmvnb(instr,outstr)
c     instr   the input character string.
c     outstr  the output character string.
c     the two strings may be the same if desired.
c***references
c***routines called    (none)
c***end prologue       rmvnb
      implicit integer(a-z)
      character*(*) instr,outstr
      character blank*1, null*1
c
      data blank/' '/, null/' '/
      save blank,null
c
c
      lin=len(instr)
      lout=len(outstr)
      j=0
      do 10 i=1,lin
         if(instr(i:i).ne.blank.and.instr(i:i).ne.null) then
            j=j+1
            if(j.le.lout) outstr(j:j)=instr(i:i)
         endif
 10   continue
      if(j.lt.lout) outstr(j+1:)=blank
c
c
      return
      end

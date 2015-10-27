*deck @(#)captlz.f	2.1  10/10/91
      subroutine captlz(instr, outstr)
c***begin prologue     captlz
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, string, capitalize
c***author             martin, richard (lanl)
c***source             @(#)captlz.f	2.1   10/10/91
c***purpose            returns an upper-case copy of an input string.
c***description
c                      call captlz(instr,outstr)
c                        instr   input character string.
c                        outstr  output character string.
c
c                      the two strings may be the same.
c
c***references
c***routines called    (none)
c***end prologue       captlz
      implicit integer (a-z)
      character instr*(*), outstr*(*)
      character lcalph*26, ucalph*26
      data lcalph/'abcdefghijklmnopqrstuvwxyz'/
      data ucalph/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
c
      lenstr=len(instr)
      if(lenstr.gt.0) then
         do 100 i=1,lenstr
            cur=index(lcalph,instr(i:i))
            if(cur.gt.0) then
               outstr(i:i)=ucalph(cur:cur)
            else
               outstr(i:i)=instr(i:i)
            endif
  100    continue
      endif
 
c
c
      return
      end

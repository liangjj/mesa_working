*deck @(#)streqc.f	5.1  11/6/94
      function streqc(instr1, instr2)
c***begin prologue     streqc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           string, equality
c***author             martin, richard (lanl)
c***source             @(#)streqc.f	5.1   11/6/94
c***purpose            convert to lower-case and compare two strings.
c***description
c                      streqc is a logical function used as:
c                        test=streqc(instr1,instr2)
c                          instr1    character string 1.
c                          instr2    character string 2.
c
c                      the strings are first converted to lower-case.
c                      streqc returns .true. if the two decapitalized strings
c                      are the same, .false. otherwise.
c                      the shorter of the two strings is blank filled.
c
c***references
c***routines called    locase(chr)
c***end prologue       streqc
      implicit integer (a-z)
      logical streqc
      character instr1*(*), instr2*(*)
      character*1 char1, char2
c
c
      streqc=.false.
      len1=len(instr1)
      len2=len(instr2)
      if(len1.eq.0.or.len2.eq.0) return
      do 10 i=1,max(len1,len2)
         if(i.gt.len1) then
            char1=' '
         else
            char1=instr1(i:i)
         endif
         if(i.gt.len2) then
            char2=' '
         else
            char2=instr2(i:i)
         endif
         call locase(char1,char1)
         call locase(char2,char2)
         if(char1.ne.char2)  return
   10 continue
         streqc=.true.
c
c
      return
      end

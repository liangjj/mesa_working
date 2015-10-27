*deck @(#)getfld.f	5.1  11/6/94
      subroutine getfld(string,cursor,result,found,delims,action,iout)
c***begin prologue     getfld
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           characters, fields, search
c***author             martin, richard (lanl)
c***source             @(#)getfld.f	5.1   11/6/94
c***purpose            returns a field within a string.  the field is anything
c                      delimited by two characters supplied by the call.
c***description
c                      call getfld(string,cursor,result,found,delims,
c                                  action,iout)
c
c                      module to search for a field within a string.
c                      a field is anything contained within two delimiting
c                      characters.
c
c                         string:  the string to be searched.
c
c                         cursor:  the input string is searched beginning
c                                  at cursor+1.  on output, cursor points
c                                  to the position of the end delimiter.
c
c                         result:  the output field. initialized to blanks.
c
c                         found:   the value of the end delimiter, or 'eor'
c                                  to signify that no complete field was
c                                  found.
c
c                         delims:  a string containing the delimiters.
c                                  the first character is the field initiator,
c                                  the second the field terminator.
c                                  a wild card symbol, '@', signifies that
c                                  any character is acceptable in the field.
c                                  thus the combination '@/' returns
c                                  everything up to (or including) the next
c                                  slash, and '/@' returns everything
c                                  after(or including) the next slash.
c
c
c                         action:  what to do with the delimiters.
c                                  'keep':    retain them in result. a wild
c                                             card delimiter is automatically
c                                             kept.
c                                  'discard': omit them in result.
c
c                         iout:    the unit number for posting messages.
c
c***references
c***routines called    (none)
c***end prologue       getfld
      implicit integer(a-z)
      character*(*) string,result,found,delims,action
      character*1 wildc,blank,id,ed
c
      data blank/' '/, wildc/'@'/
      save blank,wildc
c
c
      result=' '
      found='eor'
      lenstr=len(string)
      id=delims(1:1)
      ed=delims(2:2)
c
c     find the beginning of the field.
      if(id.eq.wildc) then
         bfld=cursor+1
         ecur=bfld
      else
         ecur=index(string(cursor+1:),id)
         bfld=cursor+ecur
      endif
      if(ecur.eq.0) return
c
c     find the end of the field.
      if(ed.eq.wildc) then
         efld=lenstr
         ecur=efld
      else
         ecur=index(string(bfld+1:),ed)
         efld=bfld+ecur
      endif
      if(ecur.eq.0) return
c
c     copy the field into the result string.
      cursor=efld
      found=ed
      if(id.ne.wildc.and.action.eq.'discard') bfld=bfld+1
      if(ed.ne.wildc.and.action.eq.'discard') efld=efld-1
      if(bfld.gt.efld) then
         efld=bfld
      endif
      result=string(bfld:efld)
c
c
      return
      end

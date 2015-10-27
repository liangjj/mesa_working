*deck @(#)dollar.f	5.1  11/6/94
      function dollar(key,rtcard,card,inp)
c***begin prologue     dollar
c***date written       870204  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           route
c***author             martin, richard (lanl)
c***source
c***purpose            locates a directive section in the input file.
c***description
c     test=dollar(key,ncards,rtcard,card,inp)
c       key     the section name to search for.
c       rtcard  an array containing the individual cards of the section.
c       card    a character scratch array.
c       inp     the input file name.
c
c     dollar is a logical function which searches through the
c     input deck for a section which begins with 'key'.
c     if the section is found, dollar=.true., and all information between
c     the key and an end delimiter ('$') is returned in rtcard.
c***references
c***routines called    positn(chr), lnkerr(mdutil)
c***end prologue       dollar
      implicit integer(a-z)
      logical dollar
      character*(*) key, rtcard, card
      character itoc*4, tmpkey*8
      logical positn
c
 1000 format(a)
c
c
      tmpkey=key
      if(positn(key,card,inp)) then
         dollar=.true.
         lencrd=len(card)
         lenrt=len(rtcard)
         mxcard=lenrt/lencrd
         ecur=0
   10    read(inp,1000,end=20) card
            if(index(card,'$').eq.0) then
               if((ecur+lencrd).gt.lenrt)
     $            call lnkerr(tmpkey//' section is too long.'
     $               //'the maximum is '//itoc(mxcard)//' cards.')
               rtcard(ecur+1:)=card
               ecur=ecur+lencrd
               goto 10
            else
            endif
      else
         dollar=.false.
         rtcard=' '
      endif
      return
c
c
   20 call lnkerr('no terminus found for '//tmpkey)
c
c
      return
      end

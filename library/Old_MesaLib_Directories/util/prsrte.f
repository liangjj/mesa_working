*deck @(#)prsrte.f	5.1  11/6/94
      subroutine prsrte(action,route,ov,ops,links,nlnks,jump,iout)
c***begin prologue     prsrte
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           parse, route
c***author             martin, richard (lanl)
c***source
c***purpose            module to parse or create a route card.
c***description
c                      mesa is guided in it's machinations by a construct
c                      called the 'route'.  the route is a series of cards,
c                      one for each 'overlay' to be executed.  an overlay is
c                      a related sequence of programs.  for example, the codes
c                      which relate to integrals are all in overlay 3. overlay
c                      3 consist of 4 links: l301 (basis set,ecp data)
c                                            l302 (one-electron integrals)
c                                            l312 (two-electron integrals)
c                                            l330 (two-electron integral sort)
c                      the route card for an overlay defines the sequence in
c                      which the links are to run, the options in force,
c                      and a directive which tells it what to do after all the
c                      links on this card have been executed.
c                      an individual route card is terminated by a semicolon,
c                      and has the general structure:
c                        overlay/options,/links(jump);
c                      for example, the atypical route:
c                        1/inau,/1;
c                        :opt
c                        2//1(:endopt);
c                        2//2;
c                        3//1,2,12,30;
c                        4//1;
c                        5/maxcyc=60,/1(:opt);
c                        :endopt
c                        20//1;
c                      directs the program to execute l1 (reading the
c                      coordinates in atomic units), l201,l202,
c                      l301,l302,l312,l330, l401 , l501 (limiting the scf to a
c                      maximum of 60 cycles),and finally l2001.
c                      if the jump is missing, the controller proceeds to the
c                      next overlay card, otherwise it skips to the card label
c                      equal to the jump.  the jump in overlay 5 above causes
c                      the program to execute l201 again after completing
c                      l501.  in this case the jumps are used to control an
c                      optimization.  l201 will override the jump to :endopt
c                      until convergence is reached via a call chainx(1),
c                      at which point mesa will execute l2001.
c
c                      call prsrte(action,route,ov,ops,links,nlnks,jump,iout)
c                        action   character string describing what to do.
c                                 'make'    construct a route card.
c                                 'crack'   intepret a route card.
c                        route    the route card, a character variable.
c                        ov       the overlay number.
c                        ops      the options character string.
c                        links    array of links to execute on this card.
c                        nlnks    the dimension of links.
c                        jump     the jump to take at the conclusion of this
c                                 card (character*16).
c                        iout     unit number of the output file.
c
c***references
c***routines called    itoc(chr),pakstr(chr),getfld(chr),lnkerr(mdutil),
c                      ctoi(chr),ffnext(chr)
c***end prologue       prsrte
      implicit integer(a-z)
      character*(*) action,route,ops,jump
      character cov*16,cln*16,itoc*16,found*16,ffnext*16
      character*1 blank,slash,lparen,rparen,semic,comma
      character*80 lnkstr
      dimension links(*)
c
      data blank/' '/, slash/'/'/, lparen/'('/, rparen/')'/, semic/';'/
      data comma/','/
      save blank,slash,lparen,rparen,semic,comma
c
c      module to parse a route card.
c
c
      if(action.eq.'make') then
         route=itoc(ov)
         icur=index(route,blank)
         ecur=index(ops,blank)
         route(icur:)=ops(:ecur)
            do 100 i=1,nlnks
               icur=index(route,blank)
               cln=itoc(links(i))
               ecur=index(cln,blank)
               route(icur:)=cln(:ecur-1)//comma
  100       continue
c        write over the comma after the last link.
         icur=index(route,blank)-1
         if(jump.ne.blank) then
            ecur=index(jump,blank)
            route(icur:)=lparen//jump(:ecur-1)//rparen
            icur=index(route,blank)
         endif
         route(icur:)=semic
c
c
      else if(action.eq.'crack') then
         call pakstr(route,lenrt)
c
c        get the overlay number.
         icur=0
         call getfld(route(1:lenrt),icur,cov,found,'@/','discard',iout)
         if(found.eq.'eor') then
            call lnkerr('no overlay found on the route card.')
         else
            ov=ctoi(cov)
         endif
c
c        get the options.
         icur=icur-1
         call getfld(route(1:lenrt),icur,ops,found,'//','discard',iout)
         if(found.eq.'eor') then
            call lnkerr('no options field found on the route card.')
         endif
c
c        get the jump.
         ilnks=icur
         call getfld(route(1:lenrt),icur,jump,found,'()','discard',
     $                  iout)
         if(found.eq.'eor') then
            jump=blank
         endif
c
c        get the links.
         if(found.eq.'eor') then
            call getfld(route(1:lenrt),ilnks,lnkstr,found,'@;',
     $                  'discard',iout)
         else
            call getfld(route(1:lenrt),ilnks,lnkstr,found,'@(',
     $                  'discard',iout)
         endif
         nlnks=0
         icur=0
  200       found=ffnext(lnkstr,icur,start,end)
            if(found.eq.'integer') then
               nlnks=nlnks+1
               links(nlnks)=ctoi(lnkstr(start:end))
               goto 200
            endif
c        check for proper syntax, just to be a hard guy.
         if(route(lenrt:lenrt).ne.semic)
     $      call lnkerr('route card must end with a semicolon.')
      endif
c
c
      return
      end

*deck @(#)chainx.f	5.1  11/6/94
      subroutine chainx(jmpovr)
c***begin prologue     chainx
c***date written       850601  (yymmdd)
c***revision date      910311  (yymmdd)
c
c    11 march 1991     rlm at lanl
c           getting rid of the /clocks/ common block.
c     2 feb   1987     rlm at lanl
c           requiring all route cards, including labels, to end with
c           a semicolon.
c     1 december 1986  pws at lanl
c           changing the name of the next links from lnnnc to mnnnc for
c           the mesa system.
c
c***keywords           chain, links, exit
c***author             martin, richard (lanl)
c***source             mdutil
c***purpose            the exit routine for mesa.  it determines which
c                      routine to run next, and sends this information
c                      to the controller.
c***description
c                      call chainx(jmpovr)
c
c                      this routine is called at the end of each program in
c                      mesa.  it's purpose is to determine which program
c                      is to run next, and to pass this information to the
c                      controller.  the sequence of programs to run is stored
c                      in route, a character string containing the route in
c                      nonstandard form.  this is documented in more detail
c                      in 'link1'.
c
c                        'chainx' is called with an argument, 'jmpovr',
c                      which allows the exiting program some control over
c                      which program runs next.
c                       jmpovr=0 ... no modification of sequence stored in
c                                    'route'.
c
c                       jmpovr=1 ... the route stored in 'route' may have
c                                    loops or skips encoded in it in the
c                                    form of forward or backward jumps,
c                                    causing control to skip to some program
c                                    other than the one stored sequentially
c                                    next in 'route'.  specifying a value
c                                    of one for 'jmpovr' causes these to be
c                                    ignored.  this is used primarily to exit
c                                    a loop at the end an optimization.
c
c
c
c***references
c
c***iosys i/o          mxcore         integer     read     1
c                      /rtstat/       integer     written
c                      "route jump"     character   written  8
c                      /clocks/       integer     written
c                      /route/        character   read     (400)
c
c***routines called    timing(ctss), dattim(mdutil), prsrte(util),
c                      msgtty(ctss), close(ctss), closep(ctss)
c***end prologue       chainx
      implicit integer(a-z)
      parameter (mxlnk=200,mxcard=50)
      integer nenter(mxlnk),lnknos(mxlnk),coruse(mxlnk)
      character*4096 newops
      character*16 chrkey,newinp
      character*80 card,route(mxcard)
      character today*24,nxtl*16,cov*16,cln*16,itoc*16,jump*8
      logical dollar
      real*8 tstop(3),elapsd(3),tstart,times(3,mxlnk)
c
      common/lnkinf/tstart(3)
      common/rtstat/ncards,curcrd,curlnk,ov,links(mxlnk),nlnks
      common/io/inp,iout
c
 1000 format(' <<<< leave link ',i4,' on ',a24,
     $       ', mxcore=',i7,',elapsed cpu/io/mem seconds:',
     $        3f8.1,' >>>>')
 1030 format(5x,80a1)
 1020 format(' new options in force:')
c
c
c     read in the common blocks which oversee the management of
c     mesa.
c        '/route/' contains the route generated in link1 coded in the
c        nonstandard form.  each "card" of route is a character
c        string delimited by semicolons.
c
c        '/rtstat/' contains information about the link just completed.
c        ncards is the total number of route cards, and curcrd is the
c        card just run.  similarly, nlnks is the total number of links
c        to execute on this card,and curlnk points to the link just run.
c        ov,links, and jump are the card parameters.
c
c        '/clocks/' contains global information regarding link timings
c        and memory usage.
c        nenter is the number of times we've encountered this link ,
c        lnknos contains the link number, coruse the memory used, and
c        times(i,.) contains the cpu,io,and memory charges, respectively.
c
c
c
c     update the timing information.
      call iosys('read integer mxcore from rwf',1,mxcore,0,' ')
      call timing(tstop(1),tstop(2),tstop(3))
c     call trakio(0,-1)
c
c     retrieve the link management information.
      call iosys('read integer "route status" from rwf',-1,ncards,0,' ')
      call iosys('read character "route jump" from rwf',8,0,0,jump)
      call iosys('read integer "links:nenter" from rwf',
     $     -1,nenter,0,' ')
      call iosys('read integer "links:lnknos" from rwf',
     $     -1,lnknos,0,' ')
      call iosys('read integer "links:coruse" from rwf',
     $     -1,coruse,0,' ')
      call iosys('read real "links:times" from rwf',
     $     -1,times,0,' ')
      call iosys('read character route from rwf',-1,0,0,route)
c
c     update the summary information.
      nchain=100*ov+links(curlnk)
c     find available space for summary in case this is a new link.
      spot=1
      do 10 i=1,mxlnk
         if(lnknos(i).ne.0) then
            spot=i+1
         endif
   10 continue
c     has this link been executed previously?
      do 20 i=1,mxlnk
         if(nchain.eq.lnknos(i)) spot=i
   20 continue
      nenter(spot)=nenter(spot)+1
      lnknos(spot)=nchain
      coruse(spot)=max(mxcore,coruse(spot))
      do 30 i=1,3
         elapsd(i)=tstop(i)-tstart(i)
         times(i,spot)=times(i,spot)+elapsd(i)
   30 continue
c
c     print the program statistics, if requested.
      prttc=0
      if(prttc.eq.1) then
         call dattim(today)
         write(iout,1000) nchain,today,mxcore,(elapsd(i),i=1,3)
      endif
c
c     what do we do next?
      newcrd=curcrd
      newlnk=curlnk+1
c     have we exhausted the links to be run on this card?
      if(newlnk.gt.nlnks) then
         newlnk=1
         if(jump.eq.' '.or.jmpovr.eq.1) then
            newcrd=newcrd+1
         else
            call pakstr(jump,lenjmp)
            newcrd=rtelbl(ncards,route,jump(1:lenjmp))
            if(newcrd.gt.ncards) call lnkerr('could not find the'
     $                          //' jump label in the route:'//jump)
            newcrd=newcrd+1
         endif
c
c        have found a new card.  is it also a label?
   40    rtstrt=cskipf(route(newcrd),' ')
         if(route(newcrd)(rtstrt:rtstrt).eq.':') then
c           two possibilities: a simple label or a goto.
            if(route(newcrd)(rtstrt:rtstrt+4).eq.':goto') then
c              accept everything until the card terminator as the jump.
               rtend=index(route(newcrd),';')
               jump=route(newcrd)(rtstrt+5:rtend-1)
               call pakstr(jump,lenjmp)
               newcrd=rtelbl(ncards,route,jump(1:lenjmp))
               if(newcrd.gt.ncards) call lnkerr('could not find the'
     $                          //' jump label in the route:'//jump)
            else
               newcrd=newcrd+1
            endif
            goto 40
         endif
      endif
c
c     perhaps there are no more cards to be run.
      if(newcrd.gt.ncards) then
c        set the program number to l998, an exit flag to the
c        controller.
         ov=9
         links(newlnk)=98
c        must someday check for 'link1" stuff.
      else
         call prsrte('crack',route(newcrd),ov,newops,links,nlnks,jump,
     $                iout)
         call iosys('write character "route jump" to rwf',8,0,0,jump)
c        if this is a new card, see if new options or memory have been
c        specified.
         if(newcrd.ne.curcrd) then
            newinp='$'//chrkey(newops,'ops','ops',' ')
            if(newinp.ne.'$ops ') then
               if(dollar(newinp,newops,card,inp)) then
                  call locase(newops,newops)
                  call pakstr(newops,lops)
                  write(iout,1020)
                  write(iout,1030) (newops(i:i),i=1,lops)
                  mxlops=len(newops)
                  if(lops.gt.mxlops)
     $               call lnkerr('new option string too long.')
                  call iosys('write character options to rwf',mxlops,
     $                        0,0,newops)
               endif
            endif
            newsiz=intkey(newops,'siz',0,' ')
            call iosys('write integer newsiz to rwf',1,newsiz,0,' ')
         endif
      endif
c
c     write the management blocks back out.
      curcrd=newcrd
      curlnk=newlnk
      call iosys('write integer "route status" to rwf',-1,ncards,0,' ')
      call iosys('write integer "links:nenter" to rwf',
     $      mxlnk,nenter,0,' ')
      call iosys('write integer "links:lnknos" to rwf',
     $      mxlnk,lnknos,0,' ')
      call iosys('write integer "links:coruse" to rwf',
     $      mxlnk,coruse,0,' ')
      call iosys('write real "links:times" to rwf',
     $      3*mxlnk,times,0,' ')
c
c     delete the link volatile files.
c     delete the overlay volatile if changing overlays.
c
c     close the internally managed files.
      call iosys('close all',0,0,0,' ')
c
c     pack the message and send it to the controller.
      cov=itoc(ov)
      cln=itoc(links(newlnk))
      lencov=cskipb(cov,' ')
      lencln=cskipb(cln,' ')
      if(lencln.eq.1) then
         nxtl='m'//cov(1:lencov)//'0'//cln(1:lencln)
      else
         nxtl='m'//cov(1:lencov)//cln(1:lencln)
      endif
      lenl=cskipb(nxtl,' ')
      call chain(nxtl(1:lenl))
c
c     close the external communication files.
c
      close (inp)
      close (iout)
c
c
      return
      end

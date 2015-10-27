*deck @(#)pm1.f	5.1  11/6/94
      program m1
c***begin prologue     pm1.f
c***date written       850601   (yymmdd)
c***revision date      11/6/94
c
c        22 may      1992    rlm at lanl
c              updating commentary
c        21 may      1991    rlm at lanl
c              making sure the route is lower case when it leaves.
c        11 march    1991    rlm at lanl
c              getting rid of the /clocks/ common block.
c        24 february 1989 by bhl at llnl
c              option to read maxsiz from ops
c        27 may 1987 by pws at lanl
c              merging with the sun compatible version.
c        4 february 1987 by rlm at lanl
c              drastically changed the nature of this link.
c              the options string is now common to all links.
c              the standard routes are read from the dat file.
c              all nonstandard route cards must end with semicolons,
c              even labels.
c        10 october 1986 by rlm at lanl
c              added one-electron property integrals(m1902) to route and
c              options.
c        24 july 1986 by pws at lanl
c              added mulliken population (m1951) to route and options
c        1 december 1986   pws at lanl
c              changing 'names' to a character array
c              fixing opens to character unit names.
c
c***keywords           m1, link 1, $route, input
c***author             martin, richard (lanl)
c***source             @(#)pm1.f	5.1   11/6/94
c***purpose            this link parses the basic route information
c
c***description
c     m1 reads the $route, $nonstd, and $title sections of the input deck.
c     in case the route is to be inferred from commands in the route section
c     it generates a proper route sequence for the links to follow.
c     it also handles intializing the read-write file in case
c     of a restart ($restart, not yet implemented).
c
c***references
c
c***routines called
c
c***end prologue       m1
      USE Set_Mesa_Parameters,    ONLY : buflen, maxcor, maxsmall,
     $                                   maxbig, maxnbfkohn, 
     $                                   maxprimkohn, maxchan, maxltop, 
     $                                   maxlmtop
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer mxcard,mxlnk
      parameter (mxcard=50,mxlnk=200)
c
      integer ncards,curcrd,curlnk,ov,links,nlnks
      integer maxsiz,lenrt,lenops,ecur,cskipb,i
      integer inp,iout
      real*8 tstart
      character*4096 route,restart,ops,refops
      character*4096 title
      character*80 card
      character rttyp*16
      character*128 names(50)
      logical dollar
c
      common/rtstat/ncards,curcrd,curlnk,ov,links(mxlnk),nlnks
      common/lnkinf/tstart(3)
      common/io/inp,iout
c
c
 1010 format(' title:')
 1020 format(5x,80a1)
 1030 format(' options:')
 1040 format(  10x,'buflen      = ',i10,1x, 
     $             'maxcor      = ',i10,1x,
     $             'maxsmall    = ',i10,/,10x,
     $             'maxbig      = ',i10,1x,
     $             'maxnbfkohn  = ',i10,1x,
     $             'maxprimkohn = ',i10,/,10x,
     $             'maxchan     = ',i10,1x,
     $             'maxltop     = 'i10,1x,
     $             'maxlmtop    = 'i10)
c 
      call m1init(names,maxsiz)
      write(iout,1040) buflen, maxcor, maxsmall, maxbig, maxnbfkohn,
     $                 maxprimkohn, maxchan, maxltop, maxlmtop
c
c     --- open the external files(input,output) ---

c
c     --- determine the type of run ---
c         dollar will return the contents of an input section if it
c         finds it.  for example, the nonstandard route comes back
c         in the array route if it is out there.
      if(dollar('$nonstd',route,card,inp)) then
         rttyp='nonstd'
      else if(dollar('$restart',restart,card,inp)) then
         rttyp='restart'
      else
         rttyp='generate'
      endif
c
c     --- retrieve the options string ---
      if(dollar('$route',ops,card,inp)) then
      else
         call lnkerr('no $route section found.')
      endif
c
c     --- this section yet to be implemented. ---
c         the intent is to provide a "reference options" list
c         which will map acceptable abbreviations onto the full
c         option list.
c     if(dollar('$refops',refops,cards,names(4)) then
c     else
c        temp=names(4)
c        call lnkerr('could not find the reference options list on:'
c    $              //temp)
c     endif
c     call iosys('write character refops on rwf',0,0,0,refops)
c     call opschk(ops,refops,iout)
c
c     --- generate the route from the options, if necessary ---
c         decapitalize the input strings which will be parsed.
      call locase(ops,ops)
      call locase(refops,refops)
      call locase(restart,restart)
      if(rttyp.eq.'generate') then
         call makrte(route,ops,iout,names(4))
      endif
      call locase(route,route)
c
c     --- initialize the rwf ---
c         note that the maximum memory allowed to the job may be specified in
c         the route input and will override the execute line replacement
c         if present.
      call inirwf(rttyp,names,maxsiz)
c
c     --- initialize some physical constants which will be used by all links,
c         and store the system common blocks.
      call setcom(rttyp,route,lenrt,ops,lenops)
c
c     --- read the title section ---
      if(dollar('$title',title,card,inp)) then
         ecur=cskipb(title,' ')
      else
         title=' '
         ecur=1
      endif
      ecur=(ecur+7)/8*8
      call iosys('write character title to rwf',ecur,0,0,title)
      write(iout,1010)
      write(iout,1020) (title(i:i),i=1,ecur)
c
c     --- print the route ---
      call prtrte(route,iout)
c
c     --- print the options string ---
      write(iout,1030)
      write(iout,1020) (ops(i:i),i=1,lenops)
c
c     exit (gracefully?).
      if(rttyp.eq.'restart') then
         call chainx(200)
      else
         call chainx(0)
      endif
c
c
      call exit
      end

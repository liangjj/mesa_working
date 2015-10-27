*deck @(#)manmem.f          
      subroutine manmem(need,pointer,ngot,name,fail)
c***begin prologue     getmem
c***date written       981110  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           memory management, dynamic memory, core adjustment
c***author             schneider b. (nsf)
c***source             @(#)manmem.f	1.1   9/6/91
c***description
c                      call manmem(need,pointer,ngot,name,fail)
c
c                      module to manage dynamic memory.
c                      input arguments:
c                         need   ... additional number of 4 byte words 
c                                    requested.
c                                ... if need is negative, the memory is freed
c                                ... and the pointer variable released.
c                         pointer... pointer to first word of variable.
c                         ngot   ... number of words actually gotten.
c                         name   ... the name of the routine requesting memory.
c                         fail   ... what to do if you cannot get the required
c                                    amount.  if .eq.0 abort.  if .ne.0 get what
c                                    you can.
c***references
c
c***iosys i/o          
c                      mxcore         integer     written   1
c
c***routines called    lnkerr(mdutil),
c***end prologue       manmem
c
      implicit integer(a-z)
c
c     note that the parameter maxcor refers to the maximum number of
c     integer words available for use.  it is set by the superuser and
c     depends on the core availability on the particular machine running
c     the mesa suite.
c
c     the parameter below is set for a machine with 512MBy of core
c
      parameter (maxcor=100000000, mincor=30000000)
c
      character*(*) name
      logical called
      common/io/inp,iout
c
      data called/.true./
      save mxcore, mncore, called
      pointer(pointer,a(1))
c
c     check if this is the first call and initialize the variables.
c
      if(called) then
c  
c        this is the first call.  set mxcore to maxcore and reset called
c        to false so that any subsequent call to the routine will understand.
c
         mxcore = maxcor
         mncore = mincor
         called=.false.
      endif
c
c     try and get the memory request
c
      if(need.eq.0) then
         ngot = mincor
      elseif(need.gt.mxcore) then
         if(fail.eq.0) then
c
c           user requested abort if need exceeds availability
c
            write(iout,1000) name, need, mxcore
            call lnkerr('user requested abort: need exceeds '
     $                     //'availability.')
         else
c
c           user requests getting as much as possible
c
            ngot = mxcore - 30000000
         endif
      else
            ngot = need
      endif
      pointer=malloc(4*ngot)
      if(pointer.ne.0) then
c
c        we can get the words. tell the user and decrement mxcore       
c
         write(iout,1000) name, ngot, mxcore
         mxcore = mxcore - ngot
         mncore = mncore - ngot
         call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
      else
c
c        we cannot get the number of requested words.  quit.
         write(iout,1000) name, ngot, mxcore
         call lnkerr('user requested abort: need exceeds '
     $                  //'availability.')
      endif
      return
 1000 format(/,1x,'link            = ',a8,/,1x,
     1          'words requested = ',i10,/,1x,
     2          'words available = ',i10)
      end

*deck @(#)getmem.f          
      subroutine getmem(need,pointer,ngot,name,fail)
c***begin prologue     getmem
c***date written       981110  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           memory management, dynamic memory, core adjustment
c***author             schneider b. (nsf)
c***source             @(#)getmem.f	1.1   9/6/91
c***description
c                      call getmem(need,pointer,ngot,name,fail)
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
c                                    amount.  if .eq.0 abort.  if .ne.0 get 
c                                    what you can.
c***references
c
c***iosys i/o          
c                      mxcore         integer     written   1
c
c***routines called    lnkerr(mdutil),
c***end prologue       getmem
c
      implicit none
c
c     note that the parameter mincore and maxcor refer to the 
c     minimum and maximum number of integer words available for use.  
c     it is set by the superuser and depends on the core availability
c     on the particular machine running the mesa suite.
c
c     mincor is set at 1GB and maxcor at 3GB of core on the Dec Alpha wit
c     4GB of memory.
c
      integer need, ngot, fail
      integer inp, iout
      integer maxcor, mxcore
      integer mincor, mncore
      logical prnt, called
      character*(*) name
      integer a
      parameter (maxcor=800000000)
      parameter (mincor=300000000)
      common/io/inp,iout
      data called/.true./
      save mxcore, mncore, called
      pointer(pointer,a(1))
      prnt=.false.
c
c     check if this is the "first" call and initialize the variables.
c
      if(called) then
         if(prnt) then
            write(iout,1) name
         endif
         if(need.lt.0) then
            call lnkerr('quit. need lt. 0 on first call')
         endif
c  
c        this is the first call.
c
         mxcore = maxcor
         mncore = mincor
         call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
         call iosys('write integer mncore on rwf',1,mncore,0,' ')
         called = .false.
      endif
      if(need.ge.0) then
         call iosys('read integer mxcore from rwf',1,mxcore,0,' ')   
         call iosys('read integer mncore from rwf',1,mncore,0,' ')   
         if(need.le.mxcore) then
c           
c           we can get the request.
c 
            ngot=need
            pointer=malloc(4*ngot)
            if(pointer.eq.0) then
               call lnkerr('malloc failure. zero pointer')
            endif
c            if(prnt) then
               write(iout,2) name, need, mxcore, ngot
c            endif
            mxcore=mxcore-ngot
            mncore=mncore-ngot
            call iosys('write integer mxcore to rwf',1,mxcore,0,' ')    
            call iosys('write integer mncore to rwf',1,mncore,0,' ')    
         else
c         
c           the request is too big.  either we quit or get what we can.
c           
            if(fail.eq.0) then
c
c              user requested abort.
c
               write(iout,2) name, need, mxcore, fail
               call lnkerr('user requested abort: need exceeds '
     $                        //'availability.')
            else
c
c              get what you can
c
               ngot = mxcore - 10000
               pointer=malloc(4*ngot)
               if(pointer.eq.0) then
                  call lnkerr('malloc failure. zero pointer')
               endif
c	       if(prnt) then
                  write(iout,2) name, need, mxcore, ngot
c	       endif	  
               mxcore=mxcore-ngot
               mncore=mncore-ngot
               call iosys('write integer mxcore to rwf',1,mxcore,0,' ')    
               call iosys('write integer mncore to rwf',1,mncore,0,' ')    
            endif
         endif
      elseif(need.lt.0) then      
c         if(prnt) then
            write(iout,3) name, abs(need)
c         endif
c
c        free the memory and restore mxcore to its proper value
         call free(pointer)
         mxcore = mxcore - need
         mncore = mncore - need
         call iosys('write integer mxcore to rwf',1,mxcore,0,' ')      
         call iosys('write integer mncore to rwf',1,mncore,0,' ')      
      endif
      return
 1    format(/,1x,'first call to memory in link = ',a8)
 2    format(/,1x,'link            = ',a8,/,1x,
     1            'words requested = ',i10,/,1x,
     2            'words available = ',i10,/,1x,
     3            'words gotten    = ',i10)
 3    format(/,1x,'free memory link = ',a8,/,1x,
     1            'words released   = ',i10)
      end











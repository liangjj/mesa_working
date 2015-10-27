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
c                                    amount.  if .eq.0 abort.  if .ne.0 get 
c                                    what you can.
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
      common/io/inp,iout
c
      data called/.true./
      save mxcore, mncore, called
      logical prnt
      pointer(pointer,a(1))
      prnt=.true.
c
c     check if this is the "first" call and initialize the variables.
c
      if(need.eq.0) then
         if(prnt) then
            write(iout,1) name
         endif
c  
c        this is the first call.  set mxcore to maxcore.
c
         mxcore = maxcor
         mncore = mincor
         call iosys('write integer mxcore on rwf',1,mxcore,0,' ')         
         call iosys('write integer mncore on rwf',1,mncore,0,' ')
c         
      elseif(need.gt.0) then
c
c        find out how much we have left
c
         call iosys('read integer mxcore from rwf',1,mxcore,0,' ')         
         call iosys('read integer mncore from rwf',1,mncore,0,' ')         
c
c     try and get the memory request
c
         if(need.gt.mxcore) then
c
c           check if its ok to get the amximum available
c           
            if(fail.ne.0) then
               need = mxcore - 10000
               ngot=need
            else
c
c              user requested abort if need exceeds availability
c
               write(iout,1000) name, need, mxcore
               call lnkerr('user requested abort: need exceeds '
     $                        //'availability.')
            endif
         else
            ngot = need
c
c           in principle we can get the words. call malloc, check the pointer
c           and if the pointer is non-zero continue.  otherwise abort.
c
            pointer=malloc(4*ngot)
            if(pointer.eq.0) then
               write(iout,1001) name, need, mxcore, pointer
               call lnkerr('malloc pointer is zero: '
     $                     //'cannot get requested memory. quit')
            else
               if(prnt) then
                  write(iout,1000) name, need, mxcore
               endif
               mxcore = mxcore - ngot
               mncore = mncore - ngot                                        
               call iosys('write integer mxcore on rwf',1,mxcore,0,' ')      
               call iosys('write integer mncore on rwf',1,mncore,0,' ')
            endif               
         endif
      elseif(need.lt.0) then
         if(prnt) then
            write(iout,2) name, abs(need)
         endif
c
c           free the memory and restore mxcore and mncore to their proper
c           values
         call free(pointer)
         mxcore = mxcore - need
         mncore = mncore - need
         call iosys('write integer mxcore on rwf',1,mxcore,0,' ')      
         call iosys('write integer mncore on rwf',1,mncore,0,' ')
      else
         call lnkerr('manmem called with bad value of need')
      endif
      return
 1    format(/,1x,'initialize link = ',a8)
 1000 format(/,1x,'link            = ',a8,/,1x,
     1            'words requested = ',i10,/,1x,
     2            'words available = ',i10)
 1001 format(/,1x,'link            = ',a8,/,1x,
     1            'words requested = ',i10,/,1x,
     2            'words available = ',i10,/,1x,
     3            'pointer         = ',i1)
 2    format(/,1x,'free memory link = ',a8,/,1x,
     1            'words released   = ',i10)
      end











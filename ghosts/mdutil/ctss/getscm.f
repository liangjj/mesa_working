*deck %W%  %G%
      subroutine getscm(need,core,ngot,name,fail)
c***begin prologue     getscm
c***date written       890114  (yymmdd)
c***revision date      891023  (yymmdd)
c        
c                      october 23,1989  rlm at lanl
c                      the codes now request all the memory available
c                      at the beginning of execution.  a common block
c                      /memory/ has been added to communicate with the
c                      major routines.  
c***keywords           memory adjustment, dynamic memory, core adjustment
c***author             lengsfield, byron (llnl)
c***source             mdutil
c***description
c                      call getscm(need,core,ngot,name,fail)
c
c                      module to increment blank common length.
c                      input arguments:
c                         need ... additional number of words requested.
c                         core ... starting address of memory region.
c                         name ... the name of the routine requesting memory.
c                                  used for posting output messages.
c                         fail ... what to do if you can't have that much.
c                                  .eq.0  ... abort if request exceeds availabilty.
c                                  .ne.0  ... get all you can.
c
c                      output arguments.
c                         ngot ... the number of words obtained.
c
c                      special cases:
c                         if need=-1, get as much core as possible.
c                         if need=0,  return possible expansion relative to 'core'
c                                     in 'ngot'.
c
c***references
c
c***iosys i/o          maxsiz         integer     read      1
c                      mxcore         integer     written   1
c
c***routines called    getufl(ctss), getfl(ctss), lnkerr(mdutil),
c                      memadj(ctss)
c***end prologue       getscm
      implicit integer(a-z)
c
      real*8 start
      common //start(2)
      common /memory/ ioff
      common/io/inp,iout
      character*(*) name
      logical called
      data called/.false./
      save maxsiz,mxcore,called,mxaval
c
c
 1000 format(1x,'  requested ',i7,' words of memory;',/,
     $       ' memadj  error condition.',i7)
c
c     recover maxsiz, the maximum working space allowed by the user.
      if(.not.called) then
         called=.true.
         call iosys('read integer maxsiz from rwf',1,maxsiz,0,' ')
c        find the start of blank common.
         lstrt=loc(start(1))
c        adjust the memory up by maxsiz.
         call memadj(maxsiz,error)  
         if(error.ne.0) then
           write(iout,1000) maxsiz,error
           call lnkerr(' specified too much core in call to getscm ')
         end if
c        the memory space will be assumed to begin at the beginning of blank
c        common.
         iaddr=lstrt
         ioff=iaddr-lstrt+1
         ngot=maxsiz
         mxaval=maxsiz
         mxcore=0
         return
      endif
c
c     ----- determine the maximum amount of additional memory possible -----
c
      blcom=loc(core)
      if (need.eq.0) then
         ngot=mxaval-(blcom-iaddr)
         return
      else if(need.eq.-1) then
         ngot=mxaval-(blcom-iaddr)
         return
      else
         ngot=mxaval-(blcom-iaddr)
      endif
c
      if(need.gt.ngot) then
            nwant=mxaval+need-ngot
            write(iout,1000) (name(i:i),i=1,2), nwant, mxaval
            call lnkerr('user requested abort: need exceeds '
     $                  //'availability.')
      endif
c
c     have parameter mxcore keep track of what was requested.
      mxcore=max(need,mxcore)
      call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
c
c
      return
      end

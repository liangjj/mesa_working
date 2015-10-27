*deck @(#)nugetscm.f	5.1  11/6/94
      subroutine getscm(need,core,ngot,name,fail)
c***begin prologue     getscm
c***date written       850601  (yymmdd)
c***revision date      900131  (yymmdd)
c
c  31 january  1990   rlm at lanl
c      modifying to request a large chunk of the memory at the beginning
c      of execution.  a common block /memory/ has been added to communicate
c      with the main programs. this is done in order to provide compatibility
c      with the developments at llnl. 
c
c  21 may 1986 pws at lanl
c              fixing bug by subtracting 'curavl' from 'mxaval' if need=-1
c              and in the check for going over the top of core. this should
c              prevent getting more core than have asked for.
c
c***keywords           memory adjustment, dynamic memory, core adjustment
c***author             martin, richard and saxe, paul (lanl)
c***source             @(#)nugetscm.f	5.1   11/6/94
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
c                                  .eq.0  ... abort if request exceeds availabil
c                                  .ne.0  ... get all you can.
c
c                      output arguments.
c                         ngot ... the number of words obtained.
c
c                      special cases:
c                         if need=-1, get as much core as possible.
c                         if need=0,  return possible expansion relative to 'cor
c                                     in 'ngot'.
c
c***references
c
c***iosys i/o          maxsiz         integer     read      1
c                      mxcore         integer     written   1
c
c***routines called    lnkerr(mdutil),
c***end prologue       getscm
c
      implicit integer(a-z)
c
c
      character*(*) name
      logical called
      dimension core(*)  
c
      common // a(2)
      common/memory/ioff
      common/io/inp,iout
c
      data called/.false./
      save maxsiz,mxcore,called,mxaval,iaddr
c
 1000 format(1x,16a1,'requested ',i7,' words of memory;',
     $       i7,' available.')
c
c
c     recover maxsiz, the maximum working space specified by the user.
c     the working space will begin at the beginning of blank common.
      if(.not.called) then
         called=.true.
         call iosys('read integer maxsiz from rwf',1,maxsiz,0,0)
c        find the start of blank common.
         lstrt=loc(a(1))
c        adjust the memory up by maxsiz.
         call malloc()
         if(error) then


         endif
c
c        return the offset address through common /memory/.
         iaddr=lstrt
         ioff=iaddr-lstrt+1
         ngot=maxsiz
         mxaval=maxsiz
         mxcore=0
         return
      endif
c
c     determine the maximum amount of additional memory possible.
c
      blcom=loc(core)
      if (need.eq.0) then
         ngot=mxaval-(blcom-iaddr)
         return
      else if (need.eq.-1) then
         ngot=mxaval-(blcom-iaddr)
         return
      else
         ngot=mxaval-(blcom-iaddr)
      endif
c
c
      if(need.eq.ngot) then
         nwant=mxaval+need-ngot
         write(iout,1000) (name(i:i),i=1,2), nwant, mxaval
         call lnkerr('user requested abort: need exceeds '
     $               //'availability.'
      endif
c
c
      have parameter mxcore keep track of what was requested.
      mxcore=max(need,mxcore)
      call iosys('write integer mxcore on rwf',1,mxcore,0,0)
c
c
      return
      end

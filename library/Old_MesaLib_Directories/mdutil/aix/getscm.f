*deck @(#)getscm.f	5.1  11/6/94
      subroutine getscm(need,core,ngot,name,fail)
c***begin prologue     getscm
c***date written       850601  (yymmdd)
c***revision date      870208  (yymmdd)
c
c   8 february 1987   pws at lanl
c      for the sun 3/50 and 3/160 under bsd 4.2 unix, making a large
c      blank common and letting paging do the job.
c
c  21 may 1986 pws at lanl
c              fixing bug by subtracting 'curavl' from 'mxaval' if need=-1
c              and in the check for going over the top of core. this should
c              prevent getting more core than have asked for.
c
c***keywords           memory adjustment, dynamic memory, core adjustment
c***author             martin, richard and saxe, paul (lanl)
c***source             @(#)getscm.f	5.1   11/6/94
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
c     note that this parameter refers to the maximum number of
c     integer words available for use.
      parameter (maxcor=5000000)
c
      character*(*) name
      logical called
      dimension core(*)  
c
      common // a(maxcor)
      common/memory/ioff
      common/io/inp,iout
c
      data called/.false./
      save maxsiz,mxcore,called
c
 1000 format(1x,16a1,'requested ',i7,' words of memory;',
     $       i7,' available.')
c
      ioff=1
      mxaval=(loc(a)-loc(core))/4+maxcor-1
      mxaval=maxcor
      curavl=mxaval
c
      if (need.eq.0) then
         ngot=mxaval-1000
         return
      end if
c
c
      nwant=need-curavl
c
      if(need.eq.-1) nwant=mxaval-curavl-1000
      if(nwant.gt.mxaval-curavl) then
         if(fail.eq.0) then
            write(iout,1000) (name(i:i),i=1,8), nwant, mxaval
            call lnkerr('user requested abort: need exceeds '
     $                  //'availability.')
         else
            nwant=mxaval
         endif
      endif
c
c     get the memory.
      ngot=curavl
      if(nwant.le.0) return
c
c
      mxcore=maxcor
      mxcore=loc(core(need))/4
      ngot=need
      call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
c
c
      return
      end

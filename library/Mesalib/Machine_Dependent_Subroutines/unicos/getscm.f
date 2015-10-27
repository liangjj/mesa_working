*deck @(#)getscm.f	5.1  11/6/94
      subroutine getscm(need,core,ngot,name,fail)
c***begin prologue     getscm
c***date written       890114  (yymmdd)
c***revision date      yymmdd  (yymmdd)
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
c                         as currently implemented, these two return
c                         identical results.
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
      parameter (maxcor=1)
c
      character*(*) name
      character*4 itoc
      logical called
c
      real*8 z
      common //z(maxcor)
      common /memory/ ioff
      common/io/inp,iout
c
c
      data called/.false./
      save maxsiz,mxcore,called,mxaval
c
c
 1000 format(1x,'  requested ',i7,' words of memory;',/,
     $       ' hpalloc  error condition.',i7)
c
c
c
c     determine the number of words available for user expansion.
c     get the starting address of the expansion region.
      blcom=loc(core)
c
c     recover maxsiz, the maximum field length allowed by the user.
      if(.not.called) then
         called=.true.
         call iosys('read integer maxsiz from rwf',1,maxsiz,0,0)
         lstrt=loc(z(1))
         call hpalloc(iaddr,maxsiz,ierr,iabort)
         if(ierr.ne.0) then
           write(iout,1000) maxsiz,ierr
           call lnkerr('hpalloc error:'//itoc(ierr))
         end if
         ioff=iaddr-lstrt+1
         ngot=maxsiz
         mxaval=maxsiz
         nused=0
         return
      endif
c
c     ----- determine the maximum amount of memory possible -----
c
      aval=mxaval-(loc(core)-loc(z(ioff)))
      if(need.eq.-1) then
        ngot=aval-10
        nused=ngot
      else if (need.eq.0) then
         ngot=aval-10
         nused=0
      else 
         if(need.gt.aval) then
            if(fail.eq.0) then
               write(iout,1000) (name(i:i),i=1,8)
               call lnkerr('user requestd abort: need exceeds'
     $                    //'availability.')
            else
               ngot=aval
               nused=mxaval
            endif
         else
            ngot=need
            nused=(loc(core)-loc(z(ioff)))+need
         endif
      end if
c
c
      mxcore=max(mxcore,nused)
      call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
c
c
      return
      end

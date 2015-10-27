*deck %W%  %G%
      subroutine getscm(need,core,ngot,name,fail)
c***begin prologue     getscm
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           memory adjustment, dynamic memory, core adjustment
c***author             martin, richard and saxe, paul (lanl)
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
      parameter (bufspace=70000)
c
      integer fcb(1)
      character*(*) name
      logical called
      common/io/inp,iout
      data called/.false./
      save maxsiz,mxcore
c
 1000 format(1x,16a1,'requested ',i7,' words of memory;',
     $       i7,' available.')
c
c     determine the number of words available for user expansion.
c     get the starting address of the expansion region.
c
      offset=1-loc(fcb)
c
      fl=and(fcb(offset+65),mask(128-24))
      hlm=and(shiftr(fcb(offset+65),24),mask(128-24))
      mfl=and(shiftr(fcb(offset+67),24),mask(128-24))-bufspace
c
      blcom=loc(core)
c     recover maxsiz, the maximum field length allowed by the user.
      if(.not.called) then
         called=.true.
         call iosys('read integer maxsiz from rwf',1,maxsiz,0,' ')
      endif
      mfl=min(mfl,maxsiz)
c
c     ----- determine the maximum amount of memory possible -----
c
      call memory('total',curavl)
      mxaval=curavl+mfl-fl+hlm-blcom
c
      if (need.eq.0) then
         ngot=mxaval-1000
         return
      end if
c
c
      nwant=need-(hlm-blcom)
      if(need.eq.-1) nwant=mxaval-(hlm-blcom)
      if(nwant.gt.mxaval-(hlm-blcom)) then
         if(fail.eq.0) then
            write(iout,1000) (name(i:i),i=1,8), nwant, mxaval-
     #                        (hlm-blcom)
            call lnkerr('user requested abort: need exceeds '
     $                  //'availability.')
         else
            nwant=mxaval-(hlm-blcom)
         endif
      endif
c
c     get the memory.
c
      if (nwant.lt.-5120000.or.nwant.gt.0) then
         call memory('uc',nwant)
      end if
c
      fl=and(fcb(offset+65),mask(128-24))
      hlm=and(shiftr(fcb(offset+65),24),mask(128-24))
      mxcore=fl
      ngot=hlm-blcom
      call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
c
c
      return
      end

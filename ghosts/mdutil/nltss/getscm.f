*deck %W%  %G%
c..bhl deck getscm
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
      real start
      common //start(2)
      common /memory/ ioff
      pointer(iaddr,b(1))
c
c
      character*(*) name
      logical called
      common/io/inp,iout
      data called/.false./
      save maxsiz,mxcore,called,mxaval
c
c
 1000 format(1x,'  requested ',i7,' words of memory;',/,
     $       ' hpalloc  error condition.',i7)
c
c     determine the number of words available for user expansion.
c     get the starting address of the expansion region.
c
c
c
      blcom=loc(core)
c     recover maxsiz, the maximum field length allowed by the user.
      if(.not.called) then
         called=.true.
         call iosys('read integer maxsiz from rwf',1,maxsiz,0,' ')
         lstrt=loc(start(1))
c         call link("unit69=terminal//")
c         write(69,*)' calling getspace maxsiz ',maxsiz
c..unicos
c..unicos         call hpalloc(iaddr,maxsiz,ierr,iabort)
c..unicos
c.ctss
         call getspace(iaddr,maxsiz)
c.ctss
         if(ierr.ne.0) then
           write(iout,1000) maxsiz,ierr
           call lnkerr(' bug call to hpalloc in getscm ')
         end if
         ioff=iaddr-lstrt+1
         ngot=maxsiz
         mxaval=maxsiz
         return
      endif
c
c     ----- determine the maximum amount of memory possible -----
c
      if (need.eq.0) then
         ngot=mxaval-(blcom-iaddr)
         return
      end if
c
c
      if(need.eq.-1) then
        ngot=mxaval-blcom+iaddr
        return
      end if
c
         ngot=mxaval-(blcom-iaddr)
c
      if(need.gt.ngot) then
             nwant=mxaval+need-ngot
            write(iout,1000) (name(i:i),i=1,2), nwant, mxaval
            call lnkerr('user requested abort: need exceeds '
     $                  //'availability.')
      endif
c
c     get the memory.
c
c      if (nwant.gt.0) then
c       call lnkerr(' bug in getscm')
c      end if
c
      call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
c
c
      return
      end

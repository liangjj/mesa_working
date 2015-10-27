*deck @(#)getscm.f	5.1  11/6/94
      subroutine getscm(need,core,ngot,name,fail)
c***begin prologue     getscm
c***date written       850601  (yymmdd)
c***revision date      910311  (yymmdd)
c
c  11 march    1991   rlm at lanl
c      calling falloc to get a block of core at job beginning.
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
c                      the way it works on the stardent is that the first
c                      call grabs a block of free memory determined
c                      by maxsiz, a parameter that is set by the user.
c                      this parameter is assumed to be in real*8 words.
c                      an offset relative to the beginning of
c                      blank common is returned, which represents
c                      the first word.  this is passed to the main
c                      programs in /memory/. 
c                      subsequent calls cannot increase
c                      the size of the block available, and simply
c                      check to see that it is enough, or return
c                      how much core is available.
c
c
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
c                         if need=0,  return possible expansion relative 
c                                     to 'core' in 'ngot'.
c
c                         as currently implemented, these two cases
c                         return identical results.
c
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
      parameter (maxcor=1)
c
      character*(*) name
      logical called
      dimension core(*)  
c
      real*8 z
      common // z(maxcor)
      common/memory/ioff
      common/io/inp,iout
c
      data called/.false./
      save called,mxaval,mxcore
c
 1000 format(1x,16a1,'requested ',i7,' words of memory;',
     $       i7,' available.')
c     
c     recover maxsiz, the field length requested by the user.
c     we assume this is in real*8 words. 
c     the routine returns information on length, however, in terms
c     of 4-byte words.
c     
      if (.not.called) then
         called=.true.
c         call iosys('read integer maxsiz from rwf',1,maxsiz,0,' ')
	maxsiz=4000000
c
c        request a block of maxsiz 8-byte words.  
c        the address of this block is addr. if can be accessed by
c        z(offset+1)
         call falloc(maxsiz,wptbyt(1),0,z(1),addr,offset)
c
         if (addr.eq.0) then
            write(iout,*) 'problem getting ',maxsiz,' real*8 words'
            call plnkerr(' falloc could not get the core requested.',
     $        666)
         end if
         ioff=offset+1
         ngot=wptoin(maxsiz)
         mxaval=wptoin(maxsiz)
         nused=0
         return
      end if
c
c
      aval=mxaval-((loc(core)-loc(z(ioff)))/itobyt(1))
      if (need.eq.-1) then
         ngot=aval-10
         nused=ngot
      else if (need.eq.0) then
         ngot=aval-10
         nused=0
      else
         if(need.gt.aval) then
            if(fail.eq.0) then
               write(iout,1000) (name(i:i),i=1,8), need, aval
               call plnkerr('user requested abort: need exceeds '
     $                     //'availability.',667)
            else
               ngot=aval
               nused=mxaval
            end if
         else
            ngot=need
            nused=(loc(core(need))-loc(z(ioff)))/itobyt(1)
         end if
      end if
c
c
      mxcore=max(mxcore,nused)
c      call iosys('write integer mxcore on rwf',1,mxcore,0,' ')
c
c
      return
      end

*deck @(#)traktm.f	1.1  7/23/91
      subroutine traktm(name,mode)
c***begin prologue     traktm
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           timing, sample, track
c***author             saxe, paul (lanl)
c***source             @(#)traktm.f	1.1   7/23/91
c***purpose            samples cpu,io,mem charges for a link.  useful in
c                      code development.
c***description
c                        name    the name of the subroutine (hollerith).
c                        mode    action to be taken.
c
c***references
c***routines called    lnkerr(mdutil), charges(ctss), timing(ctss)
c***end prologue       traktm
c
      parameter(mxstak=100,mxrtn=100,nslots=6)
c
      logical dead
      integer stak(nslots,mxstak),stats(nslots,mxrtn)
      real*8 cpusys,io,mem,cpu,sys,iomem,oldcpu,oldsys,oldmem,oldio
      real*8 cpu1,io1,mem1,sys1
c..bhl llnl
c  calling timeused instead of timing
c..bhl llnl
      integer pt
c
      common /io/     inp,iout
c
      data dead /.true./
c
c
      go to (100,200,300,400),mode+2
      call lnkerr('illegal mode in traktm')
c
  100 continue
      if (dead) return
c     fcpu=1.0/1.1
c     fsys=1.0/1.1
c     fio=1.0/0.25
c     fmem=1.0/0.25
c..llnl
      fcpu=1.0
      fsys=1.0
      fio=1.0
      fmem=1.0
c..llnl
      totcpu=0.0
      totsys=0.0
      totio=0.0
      totmem=0.0
      do 110 i=1,nrtn
         totcpu=totcpu+float(stats(3,i))*fcpu/1.0e+06
         totsys=totsys+float(stats(4,i))*fsys/1.0e+06
         totio =totio +float(stats(5,i))*fio/1.0e+06
         totmem=totmem+float(stats(6,i))*fmem/1.0e+06
  110 continue
      if (totcpu.eq.0.0) totcpu=1.0
      if (totsys.eq.0.0) totsys=1.0
      if (totio.eq.0.0) totio=1.0
      if (totmem.eq.0.0)totmem=1.0
      write (iout,101)
  101 format (//,'    ###### timing information #####',/,
     #        ' routine     calls',t25,'cpu (s)',t42,'sys (s)',
     #        t59,'io (s)',t76,'mem (mwd-s)',/)
      do 103 i=1,nrtn
         cpu=float(stats(3,i))*fcpu/1.0e+06
         sys=float(stats(4,i))*fsys/1.0e+06
         io =float(stats(5,i))*fio/1.0e+06
         mem=float(stats(6,i))*fmem/1.0e+06
         pcpu=cpu/totcpu*100.0e+00
         psys=sys/totsys*100.0e+00
         pio=io/totio*100.0e+00
         pmem=mem/totmem*100.0e+00
         write (iout,102) stats(1,i),stats(2,i),cpu,pcpu,sys,psys,
     #                 io,pio,mem,pmem
  102    format (1x,a6,i10,4(f10.4,1x,f5.1,'%'))
  103 continue
      write (iout,104) totcpu,totsys,totio,totmem
  104 format (' totals:',9x,4(f10.4,7x))
      return
c
  200 continue
      dead=.false.
      do 202 j=mxrtn
         do 201 i=nslots
            stats(i,j)=0
  201    continue
  202 continue
      do 204 j=1,mxstak
         do 203 i=1,nslots
            stak(i,j)=0
  203    continue
  204 continue
      dead=.false.
      pt=0
      nrtn=0
c     call charges(cpusys,io1,mem1)
c     call timing(cpu1,sys1,iomem)
c..llnl
      call timeused(mcpu,mio,msys)
      cpu1=float(mcpu)
      sys1=float(msys)
      io1=float(mio)
c..llnl
      return
c
  300 continue
      if (dead) return
c
c     ----- start timing of a routine -----
c
c     call charges(cpusys,io,mem)
c     call timing(cpu,sys,iomem)
      call timeused(mcpu,mio,msys)
      cpu=float(mcpu)
      sys=float(msys)
      io=float(mio)
      if (pt.gt.0) then
         stak(3,pt)=stak(3,pt)+(cpu-oldcpu)*1.0e+06
         stak(4,pt)=stak(4,pt)+(sys-oldsys)*1.0e+06
         stak(5,pt)=stak(5,pt)+(io-oldio)*1.0e+06
         stak(6,pt)=stak(6,pt)+(mem-oldmem)*1.0e+06
      end if
      pt=pt+1
      if (pt.gt.mxstak) then
         write (iout,301) mxstak
  301    format (/,' ##### to many levels for stak in timing:',i5,/)
         dead=.true.
      end if
      stak(1,pt)=name
      stak(3,pt)=0
      stak(4,pt)=0
      stak(5,pt)=0
      stak(6,pt)=0
c     call charges(cpusys,oldio,oldmem)
c     call timing(oldcpu,oldsys,iomem)
      call timeused(mcpu,mio,msys)
      oldcpu=float(mcpu)
      oldsys=float(msys)
      oldio=float(mio)
      return
c
  400 continue
      if (dead) return
c
c     ----- end of timing of routine -----
c
c     call charges(cpusys,io,mem)
c     call timing(cpu,sys,iomem)
      call timeused(mcpu,mio,msys)
      cpu=float(mcpu)
      sys=float(msys)
      io=float(mio)
      if (name.ne.stak(1,pt)) then
         write (iout,401) name,stak(1,pt)
  401    format (/,' #### error in timing routines, wrong routine',2a10,
     #           /)
         dead=.true.
      end if
c
      stak(3,pt)=stak(3,pt)+(cpu-oldcpu)*1.0e+06
      stak(4,pt)=stak(4,pt)+(sys-oldsys)*1.0e+06
      stak(5,pt)=stak(5,pt)+(io-oldio)*1.0e+06
      stak(6,pt)=stak(6,pt)+(mem-oldmem)*1.0e+06
c
      do 402 i=1,nrtn
         if (name.eq.stats(1,i)) go to 420
  402 continue
      nrtn=nrtn+1
      if (nrtn.gt.mxrtn) then
         write (iout,403) mxrtn
  403    format (/,' ##### error in timing: too many routines:',i6,/)
         dead=.true.
      end if
      stats(1,nrtn)=name
      stats(2,nrtn)=1
      stats(3,nrtn)=stak(3,pt)
      stats(4,nrtn)=stak(4,pt)
      stats(5,nrtn)=stak(5,pt)
      stats(6,nrtn)=stak(6,pt)
      go to 430
c
  420 continue
      stats(2,i)=stats(2,i)+1
      stats(3,i)=stats(3,i)+stak(3,pt)
      stats(4,i)=stats(4,i)+stak(4,pt)
      stats(5,i)=stats(5,i)+stak(5,pt)
      stats(6,i)=stats(6,i)+stak(6,pt)
c
  430 continue
      pt=pt-1
c     call charges(cpusys,oldio,oldmem)
c     call timing(oldcpu,oldsys,iomem)
      call timeused(mcpu,mio,msys)
      oldcpu=float(mcpu)
      oldsys=float(msys)
      oldio=float(mio)
      return
      end

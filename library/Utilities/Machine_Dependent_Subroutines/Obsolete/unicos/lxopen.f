*deck @(#)lxopen.f	5.1  11/6/94
      subroutine lxopen(names,mchsiz,maxsiz)
c***begin prologue     lxopen
c***date written       850601    (yymmdd)
c***revision date      910618    (yymmdd)
c
c   18 june     1991  rlm at lanl
c        splitting 'ints' into several independent units
c    1 december 1986  pws at lanl
c        1. changing mesadat file name to mesadat for mesa system
c        2. changing 'names' to a character array
c
c   21 august 1986  modified by pws at lanl
c       including the 'int' file in name replacement.
c
c***author             martin, richard (lanl)
c***source             @(#)lxopen.f	5.1   11/6/94
c***purpose            helps start-up a link.
c***description
c                      lxopen sets up default file numbers and names,
c                      processes file replacement and memory management
c                      requests, and opens the input and output
c                      files, the latter with position.
c
c                      call lxopen(names,mchsiz,maxsiz)
c                        names   hollerith file names.
c                        mchsiz  the maximum memory available on this machine.
c                        maxsiz  user defined maximum memory available to
c
c                      the correspondence in names:
c               names(1)         'inp'
c               names(2)         'out'
c               names(3)         'chk'
c               names(4)         'dat'
c               names(5)         'rwf'
c               names(6)         'rint'
c               names(7)         'int'
c               names(8)         'tint'
c               names(9)         'gint'
c               names(10)        'rdint'
c               names(11)        'dint'
c               names(12)        'zint'
c               names(13)        'ham'
c               names(14)        'moden'
c               names(15)        'aoden'
c               names(16)        'saoden'
c               names(17)        'gden'
c               names(18)        'fci'     
c     reserved for memory requests
c               names(19)        'siz'
c
c     new files added by bis for kohn work.
c               names(20)        'kohn'
c               names(21)        'kohndt'
c               names(22)        'grid'
c               names(23)        'orbs'
c               names(24)        'vstat'
c               names(25)        'ylms'
c               names(26)        'bessel'
c               names(27)        'knints'
c               names(28)        'tmat'
c               names(29)        'blktmt'
c
c***routines called    comand(mdutil), assign(ctss), openp(ctss),
c                      ctoi(chr), timing(ctss).
c***end prologue       lxopen
c
      implicit integer(a-z)
      real*8 tstart
      character*(*) names(*)
      character*8 siz
      common/lnkinf/tstart(3)
      common/io/inp,iout
      data nunits/28/
      save nunits
c
c
c     default the external file names.
      names(1)='inp'
      names(2)='out'
      names(3)='chk'
      names(4)='dat'
      names(5)='rwf'
      names(6)='rint'
      names(7)='int'
      names(8)= 'tint'
      names(9)= 'gint'
      names(10)= 'rdint'
      names(11)= 'dint'
      names(12)= 'zint'
      names(13)= 'ham'
      names(14)= 'moden'
      names(15)= 'aoden'
      names(16)= 'saoden'
      names(17)= 'gden'
      names(18)= 'fci'     
      names(19)= 'siz'
c     new files added by bis for kohn work.
      names(20)= 'kohn'
      names(21)= 'kohndt'
      names(22)= 'grid'
      names(23)= 'orbs'
      names(24)= 'vstat'
      names(25)= 'ylms'
      names(26)= 'bessel'
      names(27)= 'knints'
      names(28)= 'tmat'
      names(29)= 'blktmt'
c
c     parse the execute line message.
c     this associates an external file name('myinput'), the one known to
c     the system,  with the internal file name('inp').
c     it also searches for the the user-defined maximum field length.
c
      call comand(nunits+1,names)
c
c     ----- plug in defaults if files unchanged -----
c
      if (names(1).eq.'inp') names(1)='mesa.inp'
      if (names(2).eq.'out') names(2)='mesa.out'
      if (names(3).eq.'chk') names(3)='t/chk'
      if (names(4).eq.'dat') names(4)='/usr/tmp/mesa/mesa.dat'
      if (names(5).eq.'rwf') names(5)='t/rwf'
      if (names(6).eq.'rint') names(6)='t/rint'
      if (names(7).eq.'int') names(7)='t/int'
      if (names(8).eq.'tint') names(8)='t/tint'
      if (names(9).eq.'gint') names(9)='t/gint'
      if (names(10).eq.'rdint') names(10)='t/rdint'
      if (names(11).eq.'dint') names(11)='t/dint'
      if (names(12).eq.'zint') names(12)='t/zint'
      if (names(13).eq.'ham') names(13)='t/ham'
      if (names(14).eq.'moden') names(14)='t/moden'
      if (names(15).eq.'aoden') names(15)='t/aoden'
      if (names(16).eq.'saoden') names(16)='t/saoden'
      if (names(17).eq.'gden') names(17)='t/gden'
      if (names(18).eq.'fci') names(18)='t/fci'
c
      if (names(20).eq.'kohn') names(20)='t/kohn'
      if (names(21).eq.'kohndt') names(21)='t/kohndt'
      if (names(22).eq.'grid') names(22)='t/grid'
      if (names(23).eq.'orbs') names(23)='t/orbs'
      if (names(24).eq.'vstat') names(24)='t/vstat'
      if (names(25).eq.'ylms') names(25)='t/ylms'
      if (names(26).eq.'bessel') names(26)='t/bessel'
      if (names(27).eq.'knints') names(27)='t/knints'
      if (names(28).eq.'tmat') names(28)='t/tmat'
      if (names(29).eq.'blktmt') names(29)='t/blktmt'
c
c     initialize /io/.
c
      inp=8
      iout=9
c
c     assign the the input and output files.(output opened for appending)
c
      open (unit=inp,file=names(1),access='sequential',
     #      form='formatted',err=99,status='old')
      open (unit=iout,file=names(2),access='sequential',
     #      form='formatted',err=100,status='unknown')
      call wind(iout)
c
c     siz is the (user-defined) maximum field length allowed.
c     the default is to the full machine capability.
c
      siz=names(19)
      mchsiz=8000000
      if(siz.eq.'siz') then
         maxsiz=mchsiz
      else
         maxsiz=ctoi(siz)
         if(maxsiz.gt.mchsiz) maxsiz=mchsiz
      endif
c
c     initialize timing data.
      call timing(tstart(1),tstart(2),tstart(3))
c
c
      return
c
c     ----- errors opening the input and output files -----
c
   99 continue
         call abort('link 1 could not open the input file')
  100 continue
         call abort('link 1 could not open the output file')
  101 continue
         call lnkerr('link 1 could not open the basis-set data file')
c
c
      stop
      end

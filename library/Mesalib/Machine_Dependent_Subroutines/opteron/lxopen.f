*deck @(#)lxopen.f	5.2  4/18/95
      subroutine lxopen(names,mchsiz,maxsiz)
c***begin prologue     lxopen
c***date written       850601    (yymmdd)
c***revision date      930111   (yymmdd)
c
c   11 january  1993  rlm at lanl
c      modifying to send message to standard error when problem
c      occurs initializing file.
c    1 october  1992  rlm at lanl
c        modifying to pick up the default directory for files from
c        an environment variable
c        splitting 'ints' into several independent units
c        opening output file with 'append' access
c    1 december 1986  pws at lanl
c        1. changing mesadat file name to mesadat for mesa system
c        2. changing 'names' to a character array
c
c   21 august 1986  modified by pws at lanl
c       including the 'int' file in name replacement.
c
c***author             martin, richard (lanl)
c***source             @(#)lxopen.f	5.2   4/18/95
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
c               names(30)        'optint'
c               names(31)        'atomci'
c               names(32)        'lamdat'
c
c      new file added for bec work
c
c               names(33)        'bec'        
c
c      new file added by bis for hamiltonian manipulation
c
c               names(34)        'hconfig' 
c
c      new file added by bis for rmatrix data manipulation
c
c               names(35)        'rmtrx' 
c
c***routines called    comand(mdutil), assign(ctss), openp(ctss),
c                      ctoi(chr), timing(ctss).
c***end prologue       lxopen
c
      implicit integer(a-z)
      real*8 tstart
      character*(*) names(*)
      character*8 siz
      character*128 tmp,home
      character*8 itoc,pid
      integer stderr
c
      common/lnkinf/tstart(3)
      common/io/inp,iout
      data nunits/35/
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
      names(30)= 'optint'
      names(31)= 'atomci'
      names(32)= 'lamdat'
c     new files added by bis for bec and hamiltonian manipulation.      
      names(33)= 'bec'
      names(34)= 'hconfig'
      names(35)= 'rmtrx'
c
c     parse the execute line message.
c     this associates an external file name('myinput'), the one known to
c     the system,  with the internal file name('inp').
c     it also searches for the the user-defined maximum field length.
c
      call comand(nunits+1,names)

c
c     ----- plug in defaults if files unchanged -----
c           first get the default library names and 
c           append the process id
c           to keep one jobs files separate from another
c
      call getenv('MESA_TMP',tmp)
      ltmp=cskipb(tmp,' ')
      call getenv('MESA_HOME',home)
      lhome=cskipb(home,' ')
c  
c     someday we may want to append the process id to the name.
c     it can be done as follows.
c     id=getpid()
c     pid=itoc(id)
c     if (names(5).eq.'rwf') names(5)=tmp(1:ltmp)//'/rwf_'//pid
c
      if (names(1).eq.'inp') names(1)='mesa.inp'
      if (names(2).eq.'out') names(2)='mesa.out'
c
c     initialize /io/.
c
      inp=8
      iout=9
      ierr=stderr()
c
c     assign the the input and output files.(output opened 
c                                            for appending)
c
      open (unit=inp,file=names(1),access='sequential',
     #      form='formatted',err=99,status='old')
      open (unit=iout,file=names(2),access='append',
     #      form='formatted',err=100,status='unknown')
c
      if (names(3).eq.'chk') names(3)='mesa.chk'
      if (names(4).eq.'dat') names(4)=home(1:lhome)//'/mesa.dat'
      if (names(5).eq.'rwf') names(5)=tmp(1:ltmp)//'/rwf'
      if (names(6).eq.'rint') names(6)=tmp(1:ltmp)//'/rint'
c      if (names(6).eq.'rint') names(6)='tmp/mesa.rint'
      if (names(7).eq.'int') names(7)=tmp(1:ltmp)//'/int'
      if (names(8).eq.'tint') names(8)=tmp(1:ltmp)//'/tint'
      if (names(9).eq.'gint') names(9)=tmp(1:ltmp)//'/gint'
      if (names(10).eq.'rdint') names(10)=tmp(1:ltmp)//'/rdint'
      if (names(11).eq.'dint') names(11)=tmp(1:ltmp)//'/dint'
      if (names(12).eq.'zint') names(12)=tmp(1:ltmp)//'/zint'
      if (names(13).eq.'ham') names(13)=tmp(1:ltmp)//'/ham'
      if (names(14).eq.'moden') names(14)=tmp(1:ltmp)//'/moden'
      if (names(15).eq.'aoden') names(15)=tmp(1:ltmp)//'/aoden'
      if (names(16).eq.'saoden') names(16)=tmp(1:ltmp)//'/saoden'
      if (names(17).eq.'gden') names(17)=tmp(1:ltmp)//'/gden'
      if (names(18).eq.'fci') names(18)=tmp(1:ltmp)//'/fci'
c
      if (names(20).eq.'kohn') names(20)=tmp(1:ltmp)//'/kohn'
      if (names(21).eq.'kohndt') names(21)=tmp(1:ltmp)//'/kohndt'
      if (names(22).eq.'grid') names(22)=tmp(1:ltmp)//'/grid'
      if (names(23).eq.'orbs') names(23)=tmp(1:ltmp)//'/orbs'
      if (names(24).eq.'vstat') names(24)=tmp(1:ltmp)//'/vstat'
      if (names(25).eq.'ylms') names(25)=tmp(1:ltmp)//'/ylms'
      if (names(26).eq.'bessel') names(26)=tmp(1:ltmp)//'/bessel'
      if (names(27).eq.'knints') names(27)=tmp(1:ltmp)//'/knints'
      if (names(28).eq.'tmat') names(28)=tmp(1:ltmp)//'/tmat'
      if (names(29).eq.'blktmt') names(29)=tmp(1:ltmp)//'/blktmt'
      if (names(30).eq.'optint') names(30)=tmp(1:ltmp)//'/optint'
      if (names(31).eq.'atomci') names(31)=tmp(1:ltmp)//'/atomci'
      if (names(32).eq.'lamdat') names(32)=tmp(1:ltmp)//'/lamdat'
      if (names(33).eq.'bec')    names(33)=tmp(1:ltmp)//'/bec'
      if (names(34).eq.'hconfig') names(34)=tmp(1:ltmp)//'/hconfig'
      if (names(35).eq.'rmtrx')   names(35)=tmp(1:ltmp)//'/rmtrx'
c
c     siz is the (user-defined) maximum field length allowed.
c     the default is to the full machine capability.
c
      siz=names(19)
      mchsiz=200000000
      maxsiz=mchsiz
c     if(siz.eq.'siz') then
c        maxsiz=mchsiz
c     else
c        maxsiz=ctoi(siz)
c        if(maxsiz.gt.mchsiz) maxsiz=mchsiz
c     endif
c
c     initialize timing data.
      call timing(tstart(1),tstart(2),tstart(3))
c      write(iout,*) (names(i),i=1,28)
c
c
      return
c
c     ----- errors opening the input and output files -----
c
   99 continue
         write(ierr,*) 'link 1 could not open the input file:',
     $                 names(1)
         call abort()
  100 continue
         write(ierr,*) 'link 1 could not open the output file:',
     $                 names(2)
         call abort()
  101 continue
         write(ierr,*) 'link 1 could not open the data file:',
     $                 names(4)
         call abort()
c
c
      stop
      end

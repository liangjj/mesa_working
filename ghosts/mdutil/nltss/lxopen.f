*deck %W%  %G%
      subroutine lxopen(names,mchsiz,maxsiz)
c***begin prologue     lxopen
c***date written       850601    (yymmdd)
c***revision date      861201    (yymmdd)
c
c    1 december 1986  pws at lanl
c        1. changing mesadat file name to mesadat for mesa system
c        2. changing 'names' to a character array
c
c   21 august 1986  modified by pws at lanl
c       including the 'int' file in name replacement.
c
c***author             martin, richard (lanl)
c***source             %W%   %G%
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
c                      the corresponedence in names:
c                        names(1) ... 'inp'
c                        names(2) ... 'out'
c                        names(3) ... 'rwf'
c                        names(4) ... 'chk'
c                        names(5) ... 'dat'
c                        names(6) ... 'int'
c                        names(7) ... 'siz'
c
c***routines called    comand(mdutil), assign(ctss), openp(ctss),
c                      ctoi(chr), timing(ctss).
c***end prologue       lxopen
c
      implicit integer(a-z)
      real*8 tstart
      character*(*) names(*)
      character*8 siz
      character*8 userno,account,dropfil,suffix
c..bhl
      character*3 ans
      logical query
      character*8 name
c..bhl
      common/lnkinf/tstart(3)
      common/io/inp,iout
      data nunits/6/,runsiz/1000000/
c
      call dropfile(0)
c
c     default the external file names.
      names(1)=' inp'
      names(2)=' out'
      names(3)=' rwf'
      names(4)=' chk'
      names(5)='mesadat'
      names(6)=' int'
      names(7)='siz'
c
c     parse the execute line message.
c     this associates an external file name('myinput'), the one known to
c     the system,  with the internal file name('inp').
c     it also searches for the the user-defined maximum field length.
c
      call comand(nunits+1,names)
c
c     ----- plug in default names if needed -----
c
c
      siz=names(7)
c..bhl
      call userinfo(userno,account,dropfil,suffix)
      names(1)(1:1)=suffix(1:1)
      names(2)(1:1)=suffix(1:1)
      names(3)(1:1)=suffix(1:1)
      names(4)(1:1)=suffix(1:1)
      names(6)(1:1)=suffix(1:1)
c
c     initialize /io/.
c
      inp=8
      iout=9
c
c     assign the the input and output files.(output opened for appending)
c
c..bhl
      call link("unit6=terminal//")
      open (unit=inp,file=names(1),access='sequential',
     #      form='formatted',err=99,status='old')
      open (unit=iout,file=names(2),access='sequential',
     #      form='formatted',err=100,status='unknown')
      call wind (iout)
cps      open (unit=7,file=names(5),access='sequential',
cps     #      form='formatted',err=101,status='old')
c
c     siz is the (user-defined) maximum field length allowed.
c     the default is to the full machine capability.
c
c     call getflmx(mchsiz,idum)
 
      call izmaxfl(mchsiz)
c
c
      if(siz.eq.'siz') then
         call iosys('does maxsiz exist on rwf',0,0,0,ans)
         if(ans.eq.'no') then
          maxsiz=runsiz
         else
          call iosys('read integer maxsiz from rwf',1,maxsiz,0,' ')
         end if
      else
         maxsiz=ctoi(siz)
         if(maxsiz.gt.mchsiz) maxsiz=mchsiz
      endif
c
c     initialize timing data.
c
c
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

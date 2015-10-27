*deck @(#)liopen.f	5.1  11/6/94
      subroutine liopen(names,maxsiz)
c***begin prologue     liopen
c***date written       850601  (yymmdd)
c***revision date      910618  (yymmdd)
c
c     18 june     1991  rlm at lanl
c        splitting 'ints' file into several independent units.
c      1 december 1986  pws at lanl
c         switching 'names' to a character array and changing the open
c            of the read-write file
c
c     21 august 1986  modified by pws at lanl
c       adding 'int' file to name replacement stuff
c
c***keywords           internal files, input/output
c***author             martin, richard (lanl)
c***source
c***purpose            opens the internal files (rwf)
c***description
c                      module to open the rwf, initialize the global
c                      common area on the rwf, and save file names and memory
c
c                      call liopen(names,maxsiz)
c
c               names    names of the files known to the
c                        outside world.  the correspondence is:
c
c               names(1)     'inp'
c               names(2)     'out'
c               names(3)     'chk'
c               names(4)     'dat'
c               names(5)     'rwf'
c               names(6)     'rint'
c               names(7)     'int'
c               names(8)     'tint'
c               names(9)     'gint'
c               names(10)    'rdint'
c               names(11)    'dint'
c               names(12)    'zint'
c               names(13)    'ham'
c               names(14)    'moden'
c               names(15)    'aoden'
c               names(16)    'saoden'
c               names(17)    'gden'
c               names(18)    'fci'     
c               names(19)    'siz'
c
c     new files added by bis for kohn work.
c               names(20)    'kohn'
c               names(21)    'kohndt'
c               names(22)    'grid'
c               names(23)    'orbs'
c               names(24)    'vstat'
c               names(25)    'ylms'
c               names(26)    'bessel'
c               names(27)    'knints'
c               names(28)    'tmat'
c               names(29)    'blktmt'
c               names(30)    'optint'
c               names(31)    'atomci'
c               names(32)    'lamdat'
c               names(33)    'bec' 
c               names(34)    'hconfig'
c               names(35)    'rmtrx'
c
c***references
c
c***iosys i/o
c                      maxsiz         integer     written  1
c                      mxcore         integer     written  1
c
c***routines called    iosys(io)
c***end prologue       liopen
      implicit integer(a-z)
      character*(*) names(*)
      data mxcore/0/
      save mxcore
c
c
      call iosys('open rwf as old',1000000,0,0,names(5))
c
c     initialize the global common area on the rwf, and save file names,
c     memory and timing information.
c
      call iosys('write character "input filename" to rwf',
     $     0,0,0,names(1))
      call iosys('write character "output filename" to rwf',
     $     0,0,0,names(2))
      call iosys('write character "checkpoint filename" to rwf',
     $     0,0,0,names(3))
      call iosys('write character "data filename" to rwf',
     $     0,0,0,names(4))
      call iosys('write character "read-write filename" to rwf',
     $     0,0,0,names(5))
      call iosys('write character "raw integral filename" to rwf',
     $     0,0,0,names(6))
      call iosys('write character "integral filename" to rwf',
     $           0,0,0,names(7))
      call iosys('write character "transformed integral filename"'
     $         //' to rwf',0,0,0,names(8))
      call iosys('write character "guga integral filename" to rwf',
     $            0,0,0,names(9))
      call iosys('write character "raw derivative integral filename"'
     $            //' to rwf',0,0,0,names(10))
      call iosys('write character "derivative integral filename"'
     $            //' to rwf',0,0,0,names(11))
      call iosys('write character "zeroed integral filename" to rwf',
     $            0,0,0,names(12))
      call iosys('write character "hamiltonian filename" to rwf',
     $            0,0,0,names(13))
      call iosys('write character "mo density filename" to rwf',
     $            0,0,0,names(14))
      call iosys('write character "ao density filename" to rwf',
     $            0,0,0,names(15))
      call iosys('write character "sao density filename" to rwf',
     $            0,0,0,names(16))
      call iosys('write character "guga density filename" to rwf',
     $            0,0,0,names(17))
      call iosys('write character "ci formula filename" to rwf',
     $            0,0,0,names(18))
c
c
      call iosys('write character "kohn filename" to rwf',0,
     $            0,0,names(20))
      call iosys('write character "kohn data filename" to rwf',
     $            0,0,0,names(21))
      call iosys('write character "grid filename" to rwf',
     $            0,0,0,names(22))
      call iosys('write character "orbital filename" to rwf',
     $            0,0,0,names(23))
      call iosys('write character "potential filename" to rwf',
     $            0,0,0,names(24))
      call iosys('write character "spherical harmonic filename" to rwf',
     $            0,0,0,names(25))
      call iosys('write character "bessel function filename" to rwf',
     $            0,0,0,names(26))
      call iosys('write character "kohn integral filename" to rwf',
     $            0,0,0,names(27))
      call iosys('write character "kohn tmatrix filename" to rwf',
     $            0,0,0,names(28))
      call iosys('write character "kohn full tmatrix filename" to rwf',
     $            0,0,0,names(29))
      call iosys('write character "optical potential filename" to rwf',
     $            0,0,0,names(30))
      call iosys('write character "atomic ci filename" to rwf',
     $            0,0,0,names(31))
      call iosys('write character "linear algebraic filename" to rwf',
     $            0,0,0,names(32))  
      call iosys('write character "bec filename" to rwf',
     $            0,0,0,names(33))     
      call iosys('write character "hamiltonian manipulation '//
     $           'filename" to rwf',0,0,0,names(34))     
      call iosys('write character "r-matrix filename" to rwf',
     $                              0,0,0,names(35))     
c
c     check to see if the maximum size has been altered by
c     the options string.
c
      call iosys('read integer newsiz from rwf',1,newsiz,0,' ')
c
c     maxsiz is the maximum field length allowed to the job.
      maxsiz=max(maxsiz,newsiz)
      call iosys('write integer maxsiz to rwf',1,maxsiz,0,' ')
      call iosys('write integer mxcore to rwf',1,mxcore,0,' ')
c
c
      return
      end

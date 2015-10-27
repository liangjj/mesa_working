*deck @(#)inirwf.f	5.1  11/6/94
      subroutine inirwf(rttyp,names,maxsiz)
c***begin prologue     inirwf.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c
c    6 june     1991  rlm at lanl
c        splitting the 'ints' file into several subfiles.
c    1 december 1986  pws at lanl
c        changing 'names' and opens to character variables
c
c***keywords           rwf, i/o
c***author             martin, richard (lanl)
c***source             @(#)inirwf.f	5.1   11/6/94
c***purpose            prepares the rwf for a run.
c***description
c     call inirwf(rttyp,names,maxsiz)
c       rttyp   the type of route.
c       names   the names of the files known to the outside world.
c               element          external
c               names(1)         'inp'
c               names(2)         'out'
c               names(3)         'chk'
c               names(4)         'dat'
c               names(5)         'rwf'
c               names(6)         'rint'
c               names(7)         'int'
c               names(8)=        'tint'
c               names(9)=        'gint'
c               names(10)=       'rdint'
c               names(11)=       'dint'
c               names(12)=       'zint'
c               names(13)=       'ham'
c               names(14)=       'moden'
c               names(15)=       'aoden'
c               names(16)=       'saoden'
c               names(17)=       'gden'
c               names(18)=       'fci'     
c               names(19)=       'siz'
c
c     new files added by bis for other purposes.
c               names(20)=       'kohn'
c               names(21)=       'kohndt'
c               names(22)=       'grid'
c               names(23)=       'orbs'
c               names(24)=       'vstat'
c               names(25)=       'ylms'
c               names(26)=       'bessel'
c               names(27)=       'knints'
c               names(28)=       'tmat'
c               names(29)=       'blktmt'
c               names(30)        'optint'
c               names(31)        'atomci'
c               names(32)        'fedvr'
c               names(33)        'bec'        
c               names(34)        'hconfig' 
c               names(35)        'rmtrx' 
c               names(36)        'tdse' 
c       maxsiz  the maximum allowable field length.
c
c     this routine prepares the rwf for a run.  the action taken depends
c     on rttyp.
c***references
c***routines called    iosys(io), inicon(io), usrnam(mdutil)
c***end prologue       inirwf.f
      implicit none
c     --- input variables -----
      integer maxsiz
      character*(*) rttyp
c     --- input arrays (unmodified) ---
      character*(*) names(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout,mxcore
      character namusr*16
c
      data mxcore/0/
      save mxcore
c
      common /io/ inp,iout
c
 1000 format(5x,'user defined maxsiz:',i10)
c
c     --- open the read-write file.
      if(rttyp.eq.'restart') then
         call iosys('open rwf as old',500000,0,0,names(5))
      else
         call iosys('open rwf as new on ssd',500000,0,0,names(5))
      endif
c
c     --- initialize the global common area on the rwf, and save the
c         file names,memory and timing information.
      call iosys('write character "input filename" to rwf',
     $            0,0,0,names(1))
      call iosys('write character "output filename" to rwf',
     $            0,0,0,names(2))
      call iosys('write character "checkpoint filename" to rwf',
     $            0,0,0,names(3))
      call iosys('write character "data filename" to rwf',
     $            0,0,0,names(4))
      call iosys('write character "read-write filename" to rwf',0,0,0,
     $            names(5))
      call iosys('write character "raw integral filename" to rwf',0,0,0,
     $            names(6))
      call iosys('write character "integral filename" to rwf',0,0,
     $            0,names(7))
      call iosys('write character "transformed integral filename" to '//
     $           'rwf',0,0,0,names(8))
      call iosys('write character "guga integral filename" to rwf',0,
     $            0,0,names(9))
      call iosys('write character "raw derivative integral filename"'
     $         //' to rwf',0,0,0,names(10))
      call iosys('write character "derivative integral filename"'
     $         //' to rwf',0,0,0,names(11))
      call iosys('write character "zeroed integral filename" to rwf',
     $            0,0,0,names(12))
      call iosys('write character "hamiltonian filename" to rwf',0,
     $            0,0,names(13))
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
      call iosys('write character "fedvr filename" to rwf',
     $            0,0,0,names(32))
      call iosys('write character "bec filename" to rwf',
     $            0,0,0,names(33))     
      call iosys('write character "hamiltonian manipulation '//
     $           'filename" to rwf',0,0,0,names(34))     
      call iosys('write character "r-matrix filename" to rwf',0,0,0,
     $            names(35))    
      call iosys('write character "tdse filename" to rwf',
     $            0,0,0,names(36))
c
      write(iout,1000) maxsiz
c
      call iosys('write integer maxsiz to rwf',1,maxsiz,0,' ')
      call iosys('write integer mxcore to rwf',1,mxcore,0,' ')
c
c     --- write the user name on the rwf.
      call usrnam(namusr)
      call iosys('write character "user name" to rwf',0,0,0,namusr)
c
c     --- set the print flag to 'normal'.  it is reset in the 
c         intermediate stages of optimizations, where you don't really
c         want to see much.
      call iosys('write character "print flag" on rwf',0,0,0,'normal')
c
c
      return
      end

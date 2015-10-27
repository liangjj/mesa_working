*deck %W%  %G%
      subroutine pm903(rcore,icore)
c***begin prologue     pm903
c***date written       850101  (yymmdd) 
c***revision date      910605  (yymmdd)    
c
c  15  january, 1996   rlm at lanl
c      adding option to compute and print geminals
c   5  june,    1991   rlm at lanl
c      changes to mn903 calling list. 
c   11 january, 1991   bhl at llnl
c      unicos fixes.
c***keywords           ci 
c***author             saxe, paul 
c***source             %W%   %G%
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pm903
c
c
      implicit integer (a-z)
c
      character*4096 ops
      character*2 mcorci
      character*8 citype
      character*128 moden,tints
      logical logkey
      integer icore(*)
      real*8 rcore(*)
c
      common/io/inp,iout
c
c
      mcorci='ci'
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
      citype='m903'
      call iosys('write character "ci used" to rwf',0,0,0,citype)
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c
c     ----- call the ci routines -----
c
      maxcor=1
      call iosys('read character "transformed integral filename"'
     $          //' from rwf',0,0,0,tints)
      call iosys('open tints as old',0,0,0,tints)
c
      call iosys('read character "mo density filename"'
     $          //' from rwf',0,0,0,moden)
      call iosys('open moden as new',0,0,0,moden)
c
c     ----- ci section -----
c
      mcroot=0
      call mn903(icore,rcore,maxcor,'ci',0,0,'tints','moden',
     #           'ci',mcroot,mcorci)
c
c     --- routine to diagonalize 2pdm
      if(logkey(ops,'geminals',.false.,' ')) then
         call mn903(icore,rcore,maxcor,'density',0,0,'tints','moden',
     #              'ci',mcroot,mcorci)
      endif
c
c
c     ----- and exit with grace -----
c
      call iosys('close tints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end

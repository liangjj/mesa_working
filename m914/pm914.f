*deck  @(#)pm914.f	1.3 7/30/91
      subroutine pm914(rcore,icore)
c
c   3 december 1986   pws at lanl
c      changing 'namint' and iosys open to character.
c
      implicit integer (a-z)
c
      character*4096 ops
      character*2 mcorci
      character*8 citype
      character*128 gints,ciform
      integer icore(*)
      real*8 rcore(*)
c
c
      mcorci='ci'
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
      citype='m901'
      call iosys('write character "ci used" to rwf',0,0,0,citype)
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the ci routines -----
c
      call iosys('read character "guga integral filename" from rwf',
     $     0,0,0,gints)
c     call iosys('open gints as old',0,0,0,gints)
c
      call iosys('read character "ci formula filename" from rwf',
     $           0,0,0,ciform)
      call iosys('open ciform as new',0,0,0,ciform)
c
c
      call mn901(icore,rcore,maxcor,'ci',0,0,'guga integrals','gints')
c
c     ----- and exit with grace -----
c
      call iosys('close gints',0,0,0,' ')
      call iosys('close ciform',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end

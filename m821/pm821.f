*deck  @(#)pm821.f	5.1 11/6/94
      subroutine pm821
c
c   3 december 1986   pws at lanl
c      changing 'namgnt' and iosys open to character.
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 gints
c
c     ----- call the routines to make the supermatrices -----
c
      call iosys('read character "guga integral filename" from rwf',
     $            0,0,0,gints)
      call iosys('open gints as old',0,0,0,gints)
      call mn821('guga integrals','gints')
c
c     ----- and exit with grace -----
c
      call iosys('close gints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end

*deck @(#)pm830.f	5.1  11/6/94
      subroutine pm830(rcore,icore)
c***begin prologue     m830
c***date written       871027   (yymmdd)
c***revision date      910619   (yymmdd)
c
c   19 june    1991    rlm at lanl
c      adding an argment to mn830 and tocan so that the input density
c      matrix filename can be distinct from the output file name.
c***keywords           guga integral ordering
c***author             saxe, paul (lanl)
c***source             @(#)pm830.f	5.1   11/6/94
c
c***purpose            to sort mo ordered mo integrals to guga order.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       m830
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 gden,moden
      integer icore(*)
      real*8 rcore(*)
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the routines to make the supermatrices -----
c
c     the input density matrix file.
      call iosys('read character "guga density filename"'
     $          //' from rwf',0,0,0,gden)
      call iosys('open gden as old',0,0,0,gden)
c
c     create the output density matrix file.
      call iosys('read character "mo density filename"'
     $          //' from rwf',0,0,0,moden)
      call iosys('open moden as new',0,0,0,moden)
c
c     note that argument 6 is deliberately left blank.
      call mn830(rcore,icore,maxcor,'gden','moden',' ',
     $     'guga density matrix','mo 1pdm','mo 2pdm')
c
c
c     ----- and exit with grace -----
c
      call iosys('close gden',0,0,0,' ')
      call iosys('close moden',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end

*deck @(#)m910.f	1.4  8/3/91
      program m910
c
c***begin prologue     m910
c***date written       000811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, hamiltonian matrix
c***author             saxe, paul (lanl)
c***source             @(#)m910.f	1.4   8/3/91
c
c***purpose            to construct the hamiltonian matrix,
c                      or portions thereof.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       m910
c
c
      implicit integer (a-z)
c
c
      character*4096 ops
      common /io/ inp,iout
c
c
c     ----- open the read-write file -----
c
      call drum
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the ci routines -----
      call mn910(ops)
c
c
c     ----- and exit with grace -----
c
      call iosys('close gints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end

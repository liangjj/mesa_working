      program m1120
      implicit integer(a-z)
      character *4096 ops
      character *8 cpass
      common /io/ inp, iout
      call drum
c     ----- position input file -----
c
      call posinp ('$linalg',cpass)
c
c     ----- recover options string -----
c
      call iosys ('read character options from rwf',-1,0,0,ops)
c
c
      call onelin (ops)
      call chainx (0)
      stop
      end

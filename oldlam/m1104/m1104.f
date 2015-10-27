      program m1104
*
c----------------------------------------------------------------------c
c                                                                      c
c                    This is the original linear algebraic             c
c                    code to solve the equations from the              c
c                    input from other links                            c
c                                                                      c
c----------------------------------------------------------------------c
      common /io/ inp, iout
      character *4096 ops
      character *8 cpass
      logical logkey
      call drum
      call posinp ('$linalg',cpass)
      call iosys ('read character options from rwf',-1,0,0,ops)
*
*
*
*     ----- recover the basic data , set up the memory needed -----
*     ----- and call the routines needed for the scattering -----
*     -----                 calculation                       ----
*
      call driver (ops)
      call chainx (0)
      stop
      end

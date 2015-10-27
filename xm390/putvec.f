*deck @(#)putvec.f	2.1  10/10/91
      subroutine putvec(c,eigval,nbf)
c
c***purpose: to put the scf vector on rwf. currently handles only
c            no-sym vectors
c
c paul saxe                      20 august 1984                   lanl
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),eigval(nbf)
      character*16 xform
c
      common /io/ inp,iout
c
c     ----- store vector and eigenvalues on rwf -----
c
      call iosys('write real "scf vector" on rwf',nbf**2,c,0,' ')
      call iosys('write real "orbital energies" on rwf',nbf,eigval,
     $     0,' ')
      xform='"scf vector"'
      call iosys('write character "transformation vector" to rwf',
     $     0,0,0,xform)
c
c
      return
      end

*deck @(#)putvec.f	5.1  11/6/94
      subroutine putvec(c,eigval,nbf,noabort,usesym,orbsym,symlabl)
c
c***purpose: to put the scf vector on rwf. currently handles only
c            no-sym vectors
c
c paul saxe                      20 august 1984                   lanl
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),eigval(nbf)
      character*32 xform
      logical usesym
      integer orbsym(nbf)
      character*(*) symlabl(*)
c
c     ----- store vector and eigenvalues on rwf -----
c
      if(noabort.eq.0) then
         call iosys('write real "scf vector" on rwf',nbf**2,c,0,' ')
         call iosys('write real "orbital energies" on rwf',nbf,eigval,
     $              0,' ')
         xform='"scf vector"'
         call iosys('write character "transformation vector" to rwf',
     $               0,0,0,xform)
      else
         call iosys('write real "scf vector" on rwf',nbf**2,c,0,' ')
         call iosys('write real "orbital energies" on rwf',nbf,eigval,
     $               0,' ')
         call iosys('write real "unconverged scf vector" on rwf',
     $               nbf**2,c,0,' ')
         call iosys('write real "unconverged orbital energies" on rwf',
     $               nbf,eigval,0,' ')
         xform='"scf vector"'
         call iosys('write character "transformation vector" to rwf',
     $               0,0,0,xform)
      end if
      if(usesym) then
         call iosys('write integer "transformation vector symmetries"'
     $               //' to rwf',nbf,orbsym,0,' ')
         call iosys('write character "orbital symmetry labels" to rwf',
     $               len(symlabl(1))*nbf,0,0,symlabl)
      endif
c
c
      return
      end

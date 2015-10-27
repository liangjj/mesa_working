*deck @(#)putvec.f	5.1 11/6/94
      subroutine putvec(c,eigval,nbf,usesym,orbsym,symlabl)
c***begin prologue     putvec.f
c***date written       840820  (yymmdd)
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul (lanl) 
c***source             @(#)putvec.f	5.1   11/6/94
c***purpose            to put the scf vector on the rwf.
c                      currently handles nosym case only.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       putvec.f
      implicit none
c     --- input variables -----
      integer nbf
      logical usesym
c     --- input arrays (unmodified) ---
      real*8 c(nbf,nbf),eigval(nbf)
      integer orbsym(nbf)
      character*(*) symlabl(nbf)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      character*32 xform
      integer inp,iout
c
      common/io/inp,iout
c
c     --- store vector and eigenvalues on rwf ---
      call iosys('write real "scf vector" on rwf',nbf**2,c,0,' ')
      call iosys('write real "orbital energies" on rwf',nbf,eigval,
     $            0,' ')
      xform='"scf vector"'
      call iosys('write character "transformation vector" to rwf',
     $            0,0,0,xform)
      if(usesym) then
         call iosys('write integer "transformation vector symmetries"'
     $               //' to rwf',nbf,orbsym,0,' ')
         call iosys('write character "orbital symmetry labels" to rwf',
     $               nbf*len(symlabl(1)),0,0,symlabl) 
      endif
c
      return
      end

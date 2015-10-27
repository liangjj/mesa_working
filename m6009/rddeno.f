*deck @(#)rddeno.f	1.1 9/8/91
c***begin prologue     rddeno
c***date written       890528   (yymmdd)
c***revision date               (yymmdd)
c***keywords           read, bound
c***author             schneider, barry (lanl)
c***source             m6009
c***purpose            read in kohn energy dependent
c***                   denominator matrix
c*** 
c
c***references         none      
c
c***routines called    iosys
c***end prologue       rddeno
      subroutine rddeno(hambb,hambbc,energy,mxb,opt)
      implicit integer (a-z)
      real *8 energy, hambb
      complex *16 hambbc
      character *24 ftit
      character*16 fptoc
      character *3 opt
      common /io/ inp, iout
      dimension hambb(mxb,mxb), hambbc(mxb,mxb)
      ftit='bbdn-'//fptoc(energy)
      if (opt.eq.'yes') then
          call iosys ('read real '//ftit//' from kohnint',2*mxb*mxb,
     1        hambbc,0,' ')
      else
          call iosys ('read real '//ftit//' from kohnint',mxb*mxb,hambb,
     1                 0,' ')
      endif
      return
      end 


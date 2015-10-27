*deck @(#)wrdeno.f	1.1 9/8/91
c***begin prologue     wrdeno
c***date written       890528   (yymmdd)
c***revision date               (yymmdd)
c***keywords           write, bound
c***author             schneider, barry (lanl)
c***source             m6008
c***purpose            write out kohn energy dependent
c***                   denominator matrix
c*** 
c
c***references         none      
c
c***routines called    iosys
c***end prologue       wrdeno
      subroutine wrdeno(hambb,hambbc,energy,mxb,opt)
      implicit integer (a-z)
      real *8 energy, hambb
      complex *16 hambbc
      character *24 ftit
      character *3 opt
      character*16 fptoc
      dimension hambb(mxb,mxb), hambbc(mxb,mxb)
      ftit='bbdn-'//fptoc(energy)
      if (opt.ne.'yes') then
          call iosys ('write real '//ftit//' to kohnint',mxb*mxb,hambb,
     1                 0,' ')
      else
          call iosys ('write real '//ftit//' to kohnint',
     1                 2*mxb*mxb,hambbc, 0,' ')
      endif
      return
      end

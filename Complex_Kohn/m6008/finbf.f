*deck @(#)finbf.f	1.1 9/8/91
c***begin prologue     finbf
c***date written       890605   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           finbf, link 6008, kohn variational
c***author             schneider, barry (lanl)  
c***source             m6008
c***purpose            final bound-free integrals
c***                   
c***description        the bound-free integrals are transformed to a
c***                   free basis orthogonal to all bound molecular
c***                   orbitals.
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       finbf
      subroutine finbf (hpb,hmb,ovpb,ovmb,hambb,ntchn,matbb,prntfn)
      implicit integer(a-z)
      character *80 title
      logical prntfn
      real*8 hambb, rowv, hmb, ovmb
      complex*16 hpb, ovpb
      character *8 colt, rowt
      dimension hpb(ntchn,matbb), ovpb(ntchn,matbb), hambb(matbb,matbb)
      dimension hmb(ntchn,matbb), ovmb(ntchn,matbb)
      common /io/ inp,iout
c----------------------------------------------------------------------c
c             lets do bound-free integrals                             c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c                                                                      c
c            m(f,b) = m(f0,b) - sum <f0,b> m(b,b)                      c
c----------------------------------------------------------------------c
      call amcbc(hpb,ovpb,hambb,ntchn,matbb,matbb)
      call ambc(hmb,ovmb,hambb,ntchn,matbb,matbb)      
      if (prntfn) then
          rowv=-99.d0
          colv=-99
          title='hpb-o'
          call cmprir(hpb,rowv,colv,ntchn,matbb,ntchn,matbb,title,
     1                rowt,colt,iout)
          title='hmb-o'
          call prntrm(title,hmb,ntchn,matbb,ntchn,matbb,iout)
      endif
      return
      end




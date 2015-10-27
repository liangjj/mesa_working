*deck @(#)ci2dm.f	5.1  11/6/94
      subroutine ci2dm(twopdm,len,nv,lenblk,nij,prtpdm,index,npint,
     $     igrp,jgrp,kgrp,lgrp)
c
c***begin prologue     ci2dm
c***date written       871120   (yymmdd)
c***revision date      880104   (yymmdd)
c   january 4 1988     bhl at brl
c   changed to accumulate density in twopdm for mcscf
c
c***keywords           ci density matrix, mcscf density matrix
c***author             saxe, paul (lanl)
c***source             @(#)ci2dm.f	5.1   11/6/94
c
c***purpose            to shift this primitive block of density
c      matrix elements to 'twopdm'.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       ci2dm
c
      implicit integer (a-z)
c
      real*8 twopdm(nv,lenblk)
      real*8 prtpdm(npint,lenblk)
      real*8 fac
      integer index(len,6)
c
      common /io/ inp,iout
c
c     ----- set up a factor to account for repetitions of integrals ---
c
      if (igrp.ne.jgrp) then
         fac=2.0d+00
      else
         fac=1.0d+00
      end if
c
      if (kgrp.ne.lgrp) then
         fac=fac*2.0d+00
      end if
c
      if (igrp.ne.kgrp.or.jgrp.ne.lgrp) then
         fac=fac*2.0d+00
      end if
c
c
      do 5 i=1,nv
         ij=index(i,5)
         kl=index(i,6)
         ijkl=ij+(kl-1)*nij
         do 1 j=1,lenblk
            twopdm(i,j)=twopdm(i,j)+prtpdm(ijkl,j)*fac
 1       continue
 5    continue
c
c
c$$$      write (iout,500)
c$$$ 500  format (10x,'two-particle density matrix for ci')
c$$$      call matout(twopdm,nv,lenblk,nv,lenblk,iout)
c
c
      return
      end

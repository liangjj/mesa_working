*deck @(#)aoout.f	1.1 9/8/91
c***begin prologue     aoout
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           print
c***author             schneider, barry (lanl)
c***source             m6060
c***purpose            print routine for complex vlamda overlaps
c***
c***references         none
c
c***routines called    cebc
c***end prologue       aoout
 
      subroutine aoout(ovlmb,energy,nolam,nmo)
      implicit integer (a-z)
      real *8 energy
      complex *16 ovlmb
      character *3 itoc
      character *16 fptoc
      character *80 title
      common /io/ inp, iout 
      dimension ovlmb(nolam,nmo)
      write(iout,1)
      title='ovlmb energy-'//fptoc(energy)
      call prntcmn(title,ovlmb,nolam,nmo,nolam,nmo,iout,'e')
    1 format(//,5x,'integrals in atomic orbital representation')
      return
      end





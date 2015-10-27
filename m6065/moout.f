*deck @(#)moout.f	1.1 9/8/91
c***begin prologue     moout
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
c***end prologue       moout
 
      subroutine moout(ovlmbm,energy,nolam,nmo)
      implicit integer (a-z)
      real *8 energy
      complex *16 ovlmbm
      character *3 itoc
      character *16 fptoc
      character *80 title
      common /io/ inp, iout 
      dimension ovlmbm(nolam,nmo)
      write(iout,1)
      title='ovlmbm energy-'//fptoc(energy)
      call prntcmn(title,ovlmbm,nolam,nmo,nolam,nmo,
     1             iout,'e')
    1 format(//,5x,'integrals in molecular orbital representation')
      return
      end





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
 
      subroutine moout(ovlmbm,ovpbm,ovmbm,energy,nchan,maxlm,
     1                 nolam,nmo)
      implicit integer (a-z)
      real *8 energy
      complex *16 ovlmbm, ovpbm, ovmbm
      character *3 itoc
      character *16 fptoc
      character *80 title
      common /io/ inp, iout 
      dimension ovlmbm(nolam,nmo,nchan), ovpbm(1:maxlm,nmo,nchan)
      dimension ovmbm(1:maxlm,nmo,nchan)
      write(iout,1)
      do 10 ch=1,nchan
         title='ovlmbm-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovlmbm(1,1,ch),nolam,nmo,nolam,nmo,
     1                iout,'e')
         title='ovpbm-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovpbm(1,1,ch),maxlm,nmo,maxlm,nmo,
     1                iout,'e')
         title='ovmbm-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovmbm(1,1,ch),maxlm,nmo,maxlm,nmo,
     1                iout,'e')
   10 continue
    1 format(//,5x,'integrals in molecular orbital representation')
      return
      end





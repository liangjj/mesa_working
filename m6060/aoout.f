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
 
      subroutine aoout(ovplm,ovmlm,ovlmb,ovpb,ovmb,energy,nchan,
     1                  maxlm,nolam,nkept)
      implicit integer (a-z)
      real *8 energy
      complex *16 ovplm, ovmlm, ovlmb, ovpb, ovmb
      character *3 itoc
      character *16 fptoc
      character *80 title
      common /io/ inp, iout 
      dimension ovplm(1:maxlm,nolam,nchan), ovmlm(1:maxlm,nolam,nchan)
      dimension ovlmb(nolam,nkept,nchan), ovpb(1:maxlm,nkept,nchan)
      dimension ovmb(1:maxlm,nkept,nchan)
      write(iout,1)
      do 10 ch=1,nchan
         title='ovplm-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovplm(1,1,ch),maxlm,nolam,maxlm,nolam,
     1                iout,'e')
         title='ovmlm-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovmlm(1,1,ch),maxlm,nolam,maxlm,nolam,
     1                iout,'e')
         title='ovlmb-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovlmb(1,1,ch),nolam,nkept,nolam,nkept,
     1                iout,'e')
         title='ovpb-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovpb(1,1,ch),maxlm,nkept,maxlm,nkept,
     1                iout,'e')
         title='ovmb-chan-'//itoc(ch)//'  energy-'//fptoc(energy)
         call prntcmn(title,ovmb(1,1,ch),maxlm,nkept,maxlm,nkept,
     1                iout,'e')
   10 continue
    1 format(//,5x,'integrals in atomic orbital representation')
      return
      end





*deck @(#)prnvlm.f	1.1 9/7/91
c***begin prologue     prnvlm
c***date written       xxxxxx   (yymmdd)
c***revision date      890409   (yymmdd)
c***keywords           m6050, optical, potential
c***authors            schneider, barry (lanl)
c***source             m6050
c***purpose            print out optical potential
c***                   vlamda's
c***references       
c***routines called    iosys
c***end prologue       prnvlm
      subroutine prnvlm(vlamda,nolam,npass,pntbuf,nolst)
      implicit integer (a-z)
      complex *16 vlamda
      character *80 title
      dimension vlamda(*)
      common /io/ inp, iout
      call iosys ('rewind "complex vlamdas" on optint '//
     1            'read-and-write',0,0,0,' ')
      title='total complex w lambda'
      do 10 i=1,npass
         points=pntbuf
         wds=nolam*pntbuf
         if ( i.eq.npass) then
              points=nolst
              wds=nolam*nolst
         endif
         write(iout,1) npass
         call iosys ('read real "complex vlamdas" from optint '//
     1               'without rewinding',wds,vlamda,0,' ')
         call prntcmn(title,vlamda,points,nolam,points,nolam,iout,'e')
   10 continue
    1 format(a80)
      return
      end

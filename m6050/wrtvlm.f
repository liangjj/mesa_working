*deck @(#)wrtvlm.f	1.1 9/7/91
c***begin prologue     wrtvlm
c***date written       xxxxxx   (yymmdd)
c***revision date      890409   (yymmdd)
c***keywords           m6050, optical, potential
c***authors            schneider, barry (lanl)
c***source             m6050
c***purpose            write and print out optical potential
c***                   vlamda's
c***references       
c***routines called    iosys
c***end prologue       wrtvlm
      subroutine wrtvlm(vlamda,npnts,nolam,nwrite,nwds,ipass)
      implicit integer (a-z)
      real *8 vlamda
      dimension vlamda(*)
      common /io/ inp, iout
      wds=2*npnts*nolam
      call iosys ('write real "complex vlamdas" to optint without '//
     1               'rewinding',wds,vlamda,0,' ')
      nwds=nwds+wds
      nwrite=nwrite+1
      return
      end

*deck @(#)direct.f	5.1  11/6/94
      subroutine direct(tden,int,dirct,nnp,npass,ntriang)
      implicit integer(a-z)
      real*8 tden(nnp,npass),int(nnp,ntriang),dirct(nnp,npass)
      common /io/ inp,iout
c
c
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
c
c     ----- read through integrals if first iteration, or if
c
      maxkl=0
c
      do 2 i=1,nnp,ntriang
         minkl=maxkl+1
         maxkl=min(nnp,maxkl+ntriang)
         nread=maxkl-minkl+1
         lnread=(maxkl-minkl+1)*nnp
         call iosys('read real "sorted ao integrals" from ints '//
     $              'without rewinding',lnread,int,0,' ')
c
         call mxma(int,nnp,1,tden,1,nnp,dirct(minkl,1),1,nnp,
     $             nread,nnp,npass)
c
  2   continue
c
      return
      end

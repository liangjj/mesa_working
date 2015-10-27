*deck @(#)derint.f	5.1  11/6/94
      subroutine derint(g1,g2,buf,lbufso,nao,ldf)
      implicit integer(a-z)
c
c   add ta_der_ints to mo_der_ints
c
      real*8 g1(*),g2(*),buf(lbufso)
c
      common /io/ inp,iout
c
      nnao=nao*(nao+1)/2
      ntnao=nnao*nnao
c
      ipass=lbufso/nnao
      ipass=min(ldf,ipass)
c
      if(ipass.lt.1)call lnkerr(' m1003: derint increase lbufso')
c
      npass=(ldf-1)/ipass+1
c
      ix=0
c
      do 10 i=1,ldf,ipass
         kpass=min(ipass,ldf-i+1)
         nread=kpass*nnao
c
         call iosys('read real der_h from dints without rewinding',
     $               nread,buf,0,' ')
c..bhl
c      if(i.eq.1) then
c      write(iout,*)' derint: ndf=1 ta  h1 '
c      write(iout,1)(g1(jm),jm=1,nnao)
c      write(iout,*)' derint: ndf=1 der h1 '
c      write(iout,1)(buf(jm),jm=1,nnao)
c  1   format(5(1x,f12.8))
c      endif
c..bhl
         do 5 j=1,nread
            g1(ix+j)=g1(ix+j)+buf(j)
    5    continue
c
         ix=ix+nread
c
 10   continue
c
c
      ipass=lbufso/ntnao
      ipass=min(ldf,ipass)
      if(ipass.lt.1)call lnkerr(' m1003: derint increase lbufso')
      npass=(ldf-1)/ipass+1
c
      ix=0
c
      do 20 i=1,ldf,ipass
         kpass=min(ipass,ldf-i+1)
         nread=kpass*ntnao
c
         call iosys('read real der_g from dints without rewinding',
     $               nread,buf,0,' ')
c
c     write(iout,*)'  ta g2 '
c     write(iout,1)(g2(2*ntnao+j),j=1,ntnao)
c     write(iout,*)' der g '
c1    write(iout,1)(buf(2*ntnao+j),j=1,ntnao)
c
          do 15 j=1,nread
             g2(ix+j)=g2(ix+j)+buf(j)
   15     continue
c
         ix=ix+nread
c
 20   continue
c
      return
      end

*deck @(#)gttamo.f	1.2  7/30/91
      subroutine gttamo(tanco,tanao,scr,cm,buf,nco,nao,nob,lbufso,ldf)
c
      implicit integer(a-z)
      real*8 tanco(*),tanao(*),buf(*),scr(*),cm(*)
c
c  tanco(nob,nco,ldf)   tanao(nob,nao,ldf)
c
      nobnob=nob*nob
      nconob=nco*nob
      naonob=nao*nob
      ipass=lbufso/nobnob
      ipass=min(ldf,ipass)
      if(ipass.lt.1)call lnkerr(' m1001: gttamo increase lbufso')
c
      kx=1
      jx=1
c
      do 10 i=1,ldf,ipass
         kpass=min(ipass,ldf-i+1)
         nread=kpass*nobnob
c
         call iosys('read real mo_der_overlap from dints'//
     $              ' without rewinding',nread,buf,0,' ')
c
         ix=1
         do 11 j=1,kpass
c
            call maketa(tanco(jx),tanao(kx),buf(ix),nco,nao,nob)
c
            ix=ix+nobnob
            jx=jx+nconob
            kx=kx+naonob
  11     continue
c
 10   continue
c
      return
      end

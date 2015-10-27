*deck @(#)ytalag.f	1.3  8/8/91
      subroutine ytalag(talag,ta,ymat,cm,xlag,scr,nco,nao,noc,
     $                  nob,ndf,mxread,incorh)
c
c  talag(j,n,ndf) = talag(j,n,ndf) + ymay(j,n,i,m)*ta(j,n,ndf)
c
      implicit integer(a-z)
      real*8 ta(nob,nob),talag(nob,noc),ymat(*),cm(nob,nob)
      real*8 scr(*)
      real*8 xlag(nob,noc)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c..bhl
c     write(6,*)'  cm  '
c     call matout(cm,nob,nob,nob,nob,6)
c..bhl
c
c
      nocnob=noc*nob
      nlen=ndf*nocnob
c
      do 10 i=1,nlen
         talag(i,1)=0.0d+00
  10  continue
c
      nobnob=nob*nob
      ndfnob=ndf*nob
c
c      ix=1
c      if(nao.eq.0)go to 12
c      do 11 i=1,nao
c      call yfold(ymat(ix),scr,nob)
c      ix=ix+(i+1)*nobnob
c  11  continue
c  12  continue
c
c      ix=(nao*(nao+1)/2+nao*nco)*nobnob+1
c
c      if(nco.eq.0) go to 14
c      do 13 i=1,nco
c       call yfold(ymat(ix),scr,nob)
c       ix=ix+(i+1)*nobnob
c  13  continue
c  14  continue
c
      if(incorh.ne.0) then
         maxrd=(noc*(noc+1)/2)*nob*nob
      else
         maxrd=(nao*(nao+1)/2)*nob*nob
      endif
c
      ix=1
      lread=min(mxread,maxrd)
      call iosys('read real mc_ao_ymatrix from rwf',lread,ymat,0,' ')
      lred=lread
c
      if(nao.eq.0) go to 155
c
      naos=nco+1
      naoe=nco+nao
      do 100 i=naos,naoe
         do 200 j=naos,i
c
         if(ix.gt.lread)then
            lread=min(mxread,maxrd-lred)
            call iosys('read real mc_ao_ymatrix from rwf'//
     $                 ' without rewinding',lread,ymat,0,' ')
            lred=lred+lread
            ix=1
         endif
c
         if(i.eq.j) call yfold(ymat(ix),scr,nob)
c
c..cos      call mxmb(ymat(ix),1,nob,ta(1,j),1,nobnob,talag(1,i),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c..cos      call mxmb(ymat(ix),nob,1,ta(1,i),1,nobnob,talag(1,j),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c        call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,j),nobnob,
c    $              talag(1,i),nocnob,0,2)
c
c        call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,i),nobnob,
c    $              talag(1,j),nocnob,4,2)
         call sgemm('n','n',nob,ndf,nob,one,ymat(ix),nob,
     $              ta(1,j),nobnob,one,talag(1,i),nocnob)
c
         call sgemm('t','n',nob,ndf,nob,one,ymat(ix),nob,
     $              ta(1,i),nobnob,one,talag(1,j),nocnob)
c
         ix=ix+nobnob
c
 200     continue
100   continue
c
      if(nco.eq.0)goto 155
c
      if(incorh.eq.0) then
         ix=1
         maxrd=nco*nao*nob*nob
         lread=min(mxread,maxrd)
         call iosys('read real hess_ac from rwf',lread,ymat,0,' ')
         lred=lread
      endif
c
      do 140 i=1,nco
         do 150 j=naos,naoe
c
            if(ix.gt.lread)then
               ix=1
               lread=min(mxread,maxrd-lred)
               if(incorh.ne.0)then
                  call iosys('read real mc_ao_ymatrix from rwf'//
     $                       ' without rewinding',lread,ymat,0,' ')
               else
                  call iosys('read real hess_ac from rwf'//
     $                       ' without rewinding',lread,ymat,0,' ')
               endif
               lred=lred+lread
            endif
c
c..cos      call mxmb(ymat(ix),1,nob,ta(1,i),1,nobnob,talag(1,j),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c..cos      call mxmb(ymat(ix),nob,1,ta(1,j),1,nobnob,talag(1,i),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c        call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,j),nobnob,
c    $              talag(1,i),nocnob,0,2)
c
c        call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,i),nobnob,
c    $              talag(1,j),nocnob,4,2)
         call sgemm('n','n',nob,ndf,nob,one,ymat(ix),nob,
     $              ta(1,j),nobnob,one,talag(1,i),nocnob)
c
         call sgemm('t','n',nob,ndf,nob,one,ymat(ix),nob,
     $              ta(1,i),nobnob,one,talag(1,j),nocnob)
c
         ix=ix+nob*nob
c
  150    continue
  140 continue
c
  155 continue
c
      if(nco.eq.0) go to 250
c
      if(incorh.eq.0) then
         maxrd=(nco*(nco+1)/2)*nob*nob
         ix=1
         lread=min(mxread,maxrd)
         call iosys('read real hess_cc from rwf',lread,ymat,0,' ')
         lred=lread
      endif
c
      do 101 i=1,nco
         do 201 j=1,i
c
            if(ix.gt.lread)then
               ix=1
               lread=min(mxread,maxrd-lred)
               if(incorh.ne.0)then
                  call iosys('read real mc_ao_ymatrix from rwf'//
     $                       ' without rewinding',lread,ymat,0,' ')
               else
                  call iosys('read real hess_cc from rwf'//
     $                       ' without rewinding',lread,ymat,0,' ')
            endif
            lred=lred+lread
         endif
c
         if(i.eq.j) call yfold(ymat(ix),scr,nob)
c
c..cos      call mxmb(ymat(ix),1,nob,ta(1,j),1,nobnob,talag(1,i),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c..cos      call mxmb(ymat(ix),nob,1,ta(1,i),1,nobnob,talag(1,j),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c           call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,j),nobnob,
c    $                 talag(1,i),nocnob,0,2)
c
c           call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,i),nobnob,
c    $                 talag(1,j),nocnob,4,2)
            call sgemm('n','n',nob,ndf,nob,one,ymat(ix),nob,
     $                 ta(1,j),nobnob,one,talag(1,i),nocnob)
c
            call sgemm('t','n',nob,ndf,nob,one,ymat(ix),nob,
     $                 ta(1,i),nobnob,one,talag(1,j),nocnob)
c
            ix=ix+nob*nob
c
 201     continue
 101  continue
c
 250  continue
c
      jx=1
      do 300 i=1,ndf
         call ebtc(scr,cm,talag(jx,1),nob,nob,noc)
         call scopy(nocnob,scr,1,talag(jx,1),1)
         jx=jx+nocnob
  300 continue
c
      call addlag(talag,ta,xlag,nob,noc,ndf)
c
      return
      end

*deck @(#)xtalag.f	5.1  11/6/94
      subroutine xtalag(talag,ta,ymat,cm,xlag,scr,nco,nao,noc,
     1                  nob,ndf)
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
c
      common /io/ inp,iout
c..bhl
c      write(iout,*)'  xtalag  '
c      write(iout,*)'  mcscf ao_lagrangian   '
c      call matout(xlag,nob,noc,nob,noc,iout)
c      write(iout,*)'  ta ndf=1  '
c      call matout(ta,nob,nob,nob,nob,iout)
c..bhl
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
      ix=1
      if(nao.eq.0)go to 12
      do 11 i=1,nao
         call yfold(ymat(ix),scr,nob)
         ix=ix+(i+1)*nobnob
  11  continue
  12  continue
c
      ix=(nao*(nao+1)/2+nao*nco)*nobnob+1
c
      if(nco.eq.0) go to 14
      do 13 i=1,nco
         call yfold(ymat(ix),scr,nob)
         ix=ix+(i+1)*nobnob
  13  continue
  14  continue
c
      ix=1
c
      if(nao.eq.0) go to 155
c
      naos=nco+1
      naoe=nco+nao

      do 100 i=naos,naoe
         do 200 j=naos,i
c
c..cos      call mxmb(ymat(ix),1,nob,ta(1,j),1,nobnob,talag(1,i),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c..cos      call mxmb(ymat(ix),nob,1,ta(1,i),1,nobnob,talag(1,j),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c        call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,j),nobnob,
c    $              talag(1,i),nocnob,0,2)
         call sgemm('n','n',nob,ndf,nob,one,ymat(ix),nob,ta(1,j),nobnob,
     $              one,talag(1,i),nocnob)
c
c        call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,i),nobnob,
c    $              talag(1,j),nocnob,4,2)
         call sgemm('t','n',nob,ndf,nob,one,ymat(ix),nob,ta(1,i),nobnob,
     $              one,talag(1,j),nocnob)
c
         ix=ix+nob*nob

c
 200     continue
 100  continue
c..bhl
c      write(iout,*)'  talag ndf=1 after hess_aa  '
c      call matout(talag,nob,noc,nob,noc,iout)
c..bhl
      do 140 i=1,nco
         do 150 j=naos,naoe
c
c..cos      call mxmb(ymat(ix),1,nob,ta(1,i),1,nobnob,talag(1,j),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c..cos      call mxmb(ymat(ix),nob,1,ta(1,j),1,nobnob,talag(1,i),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c           call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,j),nobnob,
c    $                 talag(1,i),nocnob,0,2)
            call sgemm('n','n',nob,ndf,nob,one,ymat(ix),nob,
     $                  ta(1,j),nobnob,one,talag(1,i),nocnob)
c
c           call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,i),nobnob,
c    $                 talag(1,j),nocnob,4,2)
            call sgemm('t','n',nob,ndf,nob,one,ymat(ix),nob,
     $                  ta(1,i),nobnob,one,talag(1,j),nocnob)
c
            ix=ix+nob*nob
c
  150    continue
 140  continue
c..bhl
c      write(iout,*)'  talag ndf=1 after hess_ac  '
c      call matout(talag,nob,noc,nob,noc,iout)
c..bhl
  155 continue
c
      if(nco.eq.0) go to 250
c
      do 101 i=1,nco
          do 201 j=1,i
c
c..cos      call mxmb(ymat(ix),1,nob,ta(1,j),1,nobnob,talag(1,i),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c..cos      call mxmb(ymat(ix),nob,1,ta(1,i),1,nobnob,talag(1,j),1,nocnob,
c..cos     1 nob,nob,ndf)
c
c            call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,j),nobnob,
c    $                  talag(1,i),nocnob,0,2)
             call sgemm('n','n',nob,ndf,nob,one,ymat(ix),nob,
     $                  ta(1,j),nobnob,one,talag(1,i),nocnob)
c
c            call sgmm(nob,nob,ndf,ymat(ix),nob,ta(1,i),nobnob,
c    $                  talag(1,j),nocnob,4,2)
             call sgemm('t','n',nob,ndf,nob,one,ymat(ix),nob,
     $                  ta(1,i),nobnob,one,talag(1,j),nocnob)
c
             ix=ix+nob*nob
c
 201      continue
101   continue
c
 250  continue
c
      jx=1
      do 300 i=1,ndf
         call ebtc(scr,cm,talag(jx,1),nob,nob,noc)
         call scopy(nocnob,scr,1,talag(jx,1),1)
         jx=jx+nocnob
  300 continue
c..bhl
c      write(iout,*)'  talag ndf=1 after hess_cc  '
c      call matout(talag,nob,noc,nob,noc,iout)
c..bhl
      call addlag(talag,ta,xlag,nob,noc,ndf)
c
      return
      end

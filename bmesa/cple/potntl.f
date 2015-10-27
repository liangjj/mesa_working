*deck potntl
c***begin prologue     potntl
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential matrix elements
c***description        
c***references       
c
c***routines called
c***end prologue       potntl
      subroutine potntl(fci,foi,fcj,foj,vcij,pot,vij,pt,range,nci,noi,
     1                  ncj,noj,npts,n,m,type,chni,chnj,ntri,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 fci, foi, fcj, foj, vcij, vij, pot, pt, range
      character*(*) type
      dimension fci(npts,nci), foi(npts,noi)
      dimension fcj(npts,ncj), foj(npts,noj)
      dimension vcij(ntri), vij(n,m), pot(npts), pt(npts)
      logical prnt
      index=chni*(chni-1)/2+chnj
      call filpot(vcij(index),pot,pt,range,npts,chni,chnj,type)
      nbi=nci+1
      nfi=nci+noi
      nbj=ncj+1
      nfj=ncj+noj
      if (nci.ne.0.and.ncj.ne.0) then
          do 10 i=1,nci
             do 20 j=1,ncj
                vij(i,j)=0.d0
                do 30 k=1,npts
                   vij(i,j)=vij(i,j)+fci(k,i)*pot(k)*fcj(k,j)
   30           continue
   20         continue
   10     continue
      endif   
      if (nci.ne.0.and.noj.ne.0) then
          do 40 i=1,nci
             count=0
             do 50 j=nbj,nfj
                count=count+1
                vij(i,j)=0.d0
                do 60 k=1,npts
                   vij(i,j)=vij(i,j)+fci(k,i)*pot(k)*foj(k,count)
   60           continue
   50        continue
   40     continue
      endif
      if (noi.ne.0.and.ncj.ne.0) then
          counti=0
          do 70 i=nbi,nfi          
             counti=counti+1
             do 80 j=1,ncj
                vij(i,j)=0.d0
                do 90 k=1,npts
                   vij(i,j)=vij(i,j)+foi(k,counti)*pot(k)*fcj(k,j)
   90           continue
   80        continue
   70     continue                                     
      endif
      if (noi.ne.0.and.noj.ne.0) then
          counti=0
          do 100 i=nbi,nfi          
             counti=counti+1
             countj=0
             do 200 j=nbj,nfj
                countj=countj+1
                vij(i,j)=0.d0
                do 300 k=1,npts
                   vij(i,j)=vij(i,j)+
     1                      foi(k,counti)*pot(k)*foj(k,countj)
  300           continue
  200        continue
  100     continue
      endif                                  
      return
      end




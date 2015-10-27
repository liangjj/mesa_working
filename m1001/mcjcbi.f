*deck @(#)mcjcbi.f	1.2  7/30/91
      subroutine mcjcbi(grad,diag,b,t,c0,gl,olap,
     $     hess,thess,c,
     $     nmix,ncsf,mdim,energy,
     $     ndsv,i45,i46,ib2,noci,npassj,istopj,stopj,sqcdf,
     $     bufix,lbufso,rabcx,rabix,raibx,hss,incor,sg,tg,
     $     nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     buf,cr,icr,
     $     ncor2)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcjcbi.f	1.2   7/30/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cc
      integer nbf(2),nob(2),nfob(2),ncob(2),naob(2)
      integer locsym(2),len(2),lok(2),mix(2)
      integer icr(2)
      real*8 bufix(lbufso),sg(*),tg(*)
      real*8 rabcx(*),rabix(*),raibx(*),hss(*)
      real*8 c0(ncsf),hess(*),b(mdim,2),grad(*),diag(*),gl(*)
      real*8 t(mdim),olap(*),c(ncsf),thess(*)
      real*8 buf(2),cr(2),cv(2)
c
      logical debug
c
c
      common /io/ inp,iout
c
      common /number/zero,pt5,one,two,four,eight
      data debug/.false./
c
      if (debug) then
         write(iout,4011)
 4011    format(/,'   entry  jacobi ')
      end if
c
      npass=0
c
c
c       scale the problem
c
      if (debug) then
         write(iout,1071)(grad(i),i=1,mdim)
 1071    format(/,'  grad  mcjcbi '/4(2x,f14.8))
      end if
c
      zzzz=abs(grad(1))
      do 7 l=1,mdim
         if(abs(grad(l)).lt.zzzz)go to 7
         zzzz=abs(grad(l))
 7    continue
c
      aaaa=1.d0/zzzz
c temp
c      aaaa=1.d0
c      zzzz=1.d0
c temp
      do 8 l=1,mdim
         grad(l)=grad(l)*aaaa
 8    continue
 91   continue
c
      if (debug) then
         write(iout,1071)(grad(nmix+i),i=1,ncsf)
      end if
c
c     generate starting guess
c
      ltest=0
c
      if(noci.eq.1)go to 5020
      mx=nmix
      do 5015 k=1,ncsf
         diag(mx+k)=diag(mx+k)+c0(k)*c0(k)-2.d0*grad(mx+k)*c0(k)
 5015 continue
c
      if(debug) then
         write (iout,7654) (diag(k),k=1,mdim)
      end if
c
 7654 format(' *mcjcbi diag '//4(1x,f16.8))
 5020 continue
c
c
      do 6 k=1,mdim
         b(k,1)=grad(k)/diag(k)
 6    continue
c
c
 9    continue
c
      itot=0
      iter=0
c
c
 10   continue
c
c
      iter=iter+1
      iter1=iter+1
c
      if(noci.eq.1) go to 1011
c
c      project out  c0
c
      xx=0.d0
      do 1000 k=1,ncsf
         xx=xx+b(nmix+k,iter)*c0(k)
 1000 continue
      do 1010 k=1,ncsf
         b(nmix+k,iter)=b(nmix+k,iter)-xx*c0(k)
 1010 continue
 1011 continue
c
      xx=0.d0
      do 11 l=1,mdim
         xx=xx+b(l,iter)*b(l,iter)
 11   continue
      xnorm=1./sqrt(xx)
      do 12 l=1,mdim
         b(l,iter)=b(l,iter)*xnorm
 12   continue
c
c
c
c
      xx=0.d0
      do 20 k=1,mdim
         xx=xx+b(k,iter)*grad(k)
 20   continue
      gl(iter)=xx
c
c
      call mcmxvc(b(1,iter),nmix,mdim,energy,noci,iter,t,bufix,lbufso,
     $     rabcx,rabix,raibx,hss,incor,sg,tg,
     $     nsym,nbf,nob,nfob,ncob,naob,cv,c0,
     $     locsym,len,lok,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     buf,cr,icr,
     $     ncor2)
c
      iter1=iter+1
c
      if(noci.eq.1)go to 27
c
c        project out c0
c
      xx=0.d0
      do 25 k=1,ncsf
         xx=xx+c0(k)*t(nmix+k)
 25   continue
      do 26 k=1,ncsf
         t(nmix+k)=t(nmix+k)-xx*c0(k)
 26   continue
c
      if (debug) then
         write (iout,9123) (t(iii),iii=1,mdim)
      end if
c
 9123 format(' mcjcbi t projected '//4(1x,f16.8))
c
 27   continue
c
c       get new correction vector
c
 29   continue
      do 30 k=1,mdim
         b(k,iter1)=(t(k)-grad(k))/diag(k)
 30   continue
c
      if(noci.eq.1)go to 34
c
      yy=0.d0
      do 32 k=1,ncsf
         yy=yy+c0(k)*b(nmix+k,iter1)
 32   continue
      do 33 k=1,ncsf
         b(nmix+k,iter1)=b(nmix+k,iter1)-yy*c0(k)
 33   continue
 34   continue
c
      if(iter.eq.1) go to 50
c
      ipass=iter-1
      itot=ipass*(ipass+1)/2
c
      do 40 k=1,itot
         thess(k)=hess(k)
 40   continue
c
 50   continue
c
      do 70 k=1,iter
cc
         xx=0.0d0
         yy=0.0d0
         do 55 l=1,mdim
            xx=xx+b(l,k)*t(l)
            yy=yy+b(l,k)*b(l,iter1)
 55      continue
         thess(itot+k)=xx
         do 60 l=1,mdim
            b(l,iter1)=b(l,iter1)-yy*b(l,k)
 60      continue
 70   continue
c
      if(iter+1.gt.mdim)go to 86
cc
cc
      do 75 i=1,iter
cc
         xx=0.d0
         do 73 k=1,mdim
            xx=xx+b(k,i)*b(k,iter1)
 73      continue
         do 74 k=1,mdim
            b(k,iter1)=b(k,iter1)-xx*b(k,i)
 74      continue
 75   continue
      xx=0.0
      do 77 k=1,mdim
         xx=xx+b(k,iter1)*b(k,iter1)
 77   continue
      xx=1.d0/sqrt(xx)
      do 78 k=1,mdim
         b(k,iter1)=b(k,iter1)*xx
 78   continue
 86   continue
c
c
      itot=itot+iter
cc
      do 90 k=1,itot
         hess(k)=thess(k)
 90   continue
c
      if(iter.eq.1)go to 10
c
      icnvrg=0
c
      do 93 k=1,iter
         t(k)=gl(k)
 93   continue
c
      call linear(thess,t,c,iter)
c
      if (debug) then
         write (iout,9909) c(iter), olap(iter)
 9909    format(' *mcjcbi c,olap ',2(1x,f16.8))
      end if
c
      xx=0.0
      do 94 i=1,iter
         xx=xx+c(i)*c(i)
 94   continue
c
      testc=abs(c(iter))/sqrt(xx)
c
c     write(iout,1002) iter,testc
c1002 format(/,'  iter  ',i6,'  convergence ',f12.9)
c
c
      if(testc.lt.stopj) icnvrg=1
c
      if(icnvrg.eq.0.and.iter.lt.istopj)go to 10
c
      npass=npass+1
c
      do 115 l=1,mdim
         t(l)=0.d0
 115  continue
c
      do 130 k=1,iter
         do 120 l=1,mdim
            t(l)=t(l)+c(k)*b(l,k)
 120     continue
 130  continue
cc
c     write(iout,1005)
c1005 format('0 restart ')
cc
 135  continue
c
      call mcmxvc(t,nmix,mdim,energy,noci,iter,c,bufix,lbufso,
     $     rabcx,rabix,raibx,hss,incor,sg,tg,
     $     nsym,nbf,nob,nfob,ncob,naob,cv,c0,
     $     locsym,len,lok,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     buf,cr,icr,
     $     ncor2)
c
c
      if(noci.eq.1) go to 2011
c
      xx=0.d0
      do 2000 l=1,ncsf
         xx=xx+c0(l)*c(nmix+l)
 2000 continue
      do 2010 l=1,ncsf
         c(nmix+l)=c(nmix+l)-xx*c0(l)
 2010 continue
c
c
 2011 continue
c
      rmsd=0.d0
c
      do 2015 k=1,mdim
         rmsd=rmsd+(c(k)-grad(k))**2
 2015 continue
c
      rmsd=sqrt(rmsd)/float(mdim)
c
      if(rmsd.gt.stopj) then
         write(iout,181) rmsd,stopj
 181     format('0 rmsd = ',f15.10,' toler = ',f15.10)
         call lnkerr('linear equations failed to converge')
      endif
c
c        scale up the solution
c
      do 136 l=1,mdim
         t(l)=t(l)*zzzz
 136  continue
c
      sqcdf=0.d0
      do 137 l=1,mdim
         sqcdf=sqcdf+t(l)*t(l)
 137  continue
c
      if(noci.eq.1) return
c
c      calc. overlap of correction vector with  c(0)
c
c     xx=0.0d0
c     do 145 l=1,ncsf
c     xx=xx+c0(l)*t(nmix+l)
c145  continue
c
c     write(iout,146) xx
c146  format(//,'  <&cic> =  ',f12.9)
c
c     xx=0.d0
c     do 147 l=1,ncsf
c     xx=xx+t(nmix+l)*t(nmix+l)
c147  continue
c
c     xx=1.d0-xx/2.d0
c
c     do 150 l=1,ncsf
c     c(l)=xx*c0(l)-t(nmix+l)
c150  continue
c     zz=0.d0
c     do 151 l=1,ncsf
c     zz=zz+c(l)*c(l)
c151  continue
c
c     zz=1.d0/sqrt(zz)
c     do 155 l=1,ncsf
c     c(l)=c(l)*zz
c155  continue
c
      return
 200  continue
c     write(iout,210)
 210  format(/,'  stop  no solution to linear eqn.s')
      call lnkerr(' ')
      end

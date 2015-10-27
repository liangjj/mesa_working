*deck @(#)tocmplx.f
      subroutine tocmplx(eigc,vlc,vrc,eigr,eigi,vl,vr,ovlp,n,prn)
c***begin prologue     tocmplx
c***date written       
c***revision date      yymmdd  (yymmdd)
c***keywords           unsymmetric eigenvalue
c                      complex
c***author             schneider, barry (nsf)
c***source             
c***purpose            make real eigenvalue and eigenvector components 
c***description
c
c***references
c***routines called    
c***end prologue       tocmplx
      implicit integer(a-z)
c
      real*8 eigr, eigi, vl, vr
      complex*16 eigc, vlc, vrc, ovlp, eye, norm, cdotc
      character*80 title
      logical prn
      dimension eigc(n), vlc(n,n), vrc(n,n)
      dimension eigr(n), eigi(n), vl(n,n), vr(n,n)
      dimension ovlp(n,n)
      data eye/(0.d0,1.d0)/
      common/io/ inp, iout
      i=1 
      do while(i.le.n)
         if(eigi(i).eq.0.d0) then
            eigc(i)=eigr(i)
            do 10 k=1,n
               vlc(k,i)=vl(k,i)
               vrc(k,i)=vr(k,i)
 10         continue
            norm=sqrt(1.d0/cdotc(n,vlc(1,i),1,vrc(1,i),1))
            call cscal(n,norm,vrc(1,i),1)
            norm=conjg(norm)
            call cscal(n,norm,vlc(1,i),1)
            i=i+1
         else
            do 20 k=1,n
               vlc(k,i)=vl(k,i)+eye*vl(k,i+1)
               vlc(k,i+1)= vl(k,i)-eye*vl(k,i+1)
               vrc(k,i)=vr(k,i)+eye*vr(k,i+1)
               vrc(k,i+1)= vr(k,i)-eye*vr(k,i+1)
 20         continue
            eigc(i)=eigr(i)+eye*eigi(i)     
            eigc(i+1)=eigr(i+1)+eye*eigi(i+1)     
            norm=sqrt(1.d0/cdotc(n,vlc(1,i),1,vrc(1,i),1))
            call cscal(n,norm,vrc(1,i),1)
            norm=conjg(norm)
            call cscal(n,norm,vlc(1,i),1)
            norm=sqrt(1.d0/cdotc(n,vlc(1,i+1),1,vrc(1,i+1),1))
            call cscal(n,norm,vrc(1,i+1),1)
            norm=conjg(norm)
            call cscal(n,norm,vlc(1,i+1),1)
            i=i+2   
         endif
      enddo
      if(prn) then
         title='eigenvalues in complex form'
         call prntcm(title,eigc,n,1,n,1,iout)
         title='left eigenvectors in complex form'
         call prntcm(title,vlc,n,n,n,n,iout)
         title='right eigenvectors in complex form'
         call prntcm(title,vrc,n,n,n,n,iout)
         call cehbtc(ovlp,vlc,vrc,n,n,n)
         title='overlap matrix'
         call prntcm(title,ovlp,n,n,n,n,iout)
      endif
	 
c
      return
      end

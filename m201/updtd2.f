*deck @(#)updtd2.f	5.1  11/6/94
      subroutine updtd2(nvar,nvv,mxpt,xx,x,ff,f,fc,fs,ic,fcold,fcnew,
     $                  a,is,ipsav1,np,rmax,rmin,rlim,grderr)
c***begin prologue     updtd2.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)updtd2.f	5.1   11/6/94
c***purpose            
c***description
c     update the second derivative matrix using values of the first
c     derivative computed here and there on the surface.
c***references
c
c***routines called
c
c***end prologue       updtd2.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,mxpt,np
c     --- input arrays (unmodified) ---
      integer is(mxpt)
      real*8 x(nvar),f(nvar)
      real*8 fc(2*nvv),fs(mxpt)
c     --- input arrays (scratch) ---
      real*8 a(nvar,mxpt)
c     --- output arrays ---
      integer ic(mxpt)
      real*8 fcold(nvar),fcnew(nvar)
      real*8 xx(nvar,mxpt),ff(nvar,mxpt)
c     --- output variables ---
      integer ipsav1
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer lind,i,j,ij,ii,ip,nn,jp
      real*8 zero,dx,ddx,rx,deltak,sdot,error
      real*8 rmax,rmin,rlim,grderr
c
      parameter (zero=0.0d+00)
c
      common/io/inp,iout
c
 1000 format(5x,'update second derivatives using information from',
     $       ' cycles',5i3,(/5x,53x,5i3))
c
      lind(i,j)=(0.5*max(i,j)*(max(i,j)-1)) +min(i,j)
c
c
      ipsav1=0
      do 10 j=2, np
         call vsub(a(1,j-1),xx(1,j),xx(1,1),nvar)
   10 continue
      nn=np-1
      call gshmdt(a,is(2),nvar,nvar,nn,rmax,rmin,rlim)
c     gram-schmidt may have reduced the number of linearly independent
c     points.
      nn=nn+1
c
      do 70 ip=2,nn
         ii=is(ip)+1
         dx=zero
         rx=zero
         do 30 i=1,nvar
            ddx=(xx(i,ii)-x(i))
            dx=dx+ddx**2
            rx=rx+ddx*a(i,ip-1)
            fcold(i)=zero
            fcnew(i)=f(i)-ff(i,ii)
            do 20 j=1,nvar
               fcold(i)=fcold(i)+fc(lind(i,j))*(xx(j,ii)-x(j))
   20       continue
   30    continue
c
         dx=sqrt(dx)
         do 60 jp=ip,nn
            deltak=sdot(nvar,fcnew,1,a(1,jp-1),1)
     $             -sdot(nvar,fcold,1,a(1,jp-1),1)
            deltak=deltak/rx
c
c           test numerical accuracy and update the hessian.
            error=grderr*dx/rx
            if(abs(deltak).ge.error) then
               ij=0
               do 50 i=1,nvar
                  do 40 j=1,i
                     ij=ij+1
                     if(ip.ne.jp)
     $                     fc(ij)=fc(ij)+a(i,jp-1)*deltak*a(j,ip-1)
                     fc(ij)=fc(ij)+a(i,ip-1)*deltak*a(j,jp-1)
   40             continue
   50          continue
            endif
   60    continue
   70 continue
c
c     eliminate unsuitable points.
      do 90 ii=2,nn
         ip=is(ii)+1
         fs(ii)=fs(ip)
         ic(ii)=ic(ip)
         call vmove(xx(1,ii),xx(1,ip),nvar)
         call vmove(ff(1,ii),ff(1,ip),nvar)
   90 continue
c
      np=nn
      write(iout,1000) (ic(nn-i+1),i=1,nn)
c
c
      return
      end

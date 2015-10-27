*deck @(#)bsquad.f	5.1  11/6/94
      subroutine bsquad(nvar,nvv,ftemp,xnew,dxrms,dxmax,ok,nnege,
     $                  fc,vname,eigen,a,vec,xxx)
c***begin prologue     bsquad.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)bsquad.f	5.1   11/6/94
c***purpose            performs quadratic search
c***description
c     the current second derivative matrix and gradient are used to
c     compute the displacement from the current point needed to find
c     the minimum on a purely quadratic surface.   this is called the
c     quadratic search and uses the well known newton-raphson formula.
c
c***references
c
c***routines called
c
c***end prologue       bsquad.f
      implicit none
c     --- input variables -----
      integer nvar,nvv
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
      real*8 ftemp(nvar),vec(nvar,nvar),fc(nvv)
      real*8 eigen(nvar),xxx(5*nvar),a(nvv)
c     --- output arrays ---
      real*8 xnew(nvar)
c     --- output variables ---
      integer nnege
      logical ok
      real*8 dxrms,dxmax
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,bscycl,np,neg,ipsav1,ipsav2
      integer i,j,k,lind,ii,ij,ji
      logical prnt,tstcrv,updrwf,chkpt,clnup
      real*8 energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,eigmax,eigmin
      real*8 fswtch,fncerr,grderr,fnccnv
      real*8 dx,ddx,rx
      real*8 zero
c
      parameter (zero=0.0d+00)
c
      common/io/inp,iout
      common/bsinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fswtch,fncerr,grderr,fnccnv,
     $              mxcycl,bscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
c
 1010 format(5x,'eigenvectors of the second derivative matrix:')
 1020 format(5x,'eigenvectors required to have negative eigenvalues.')
 1030 format(5x,'eigenvalue',i3,' out of range, new value=',f12.6,'.',
     $      /5x,'eigenvector:')
 1040 format(5x,'eigenvalues -- ',1x,5f10.5)
c
      lind(i,j)=(0.5d00*max(i,j)*(max(i,j)-1)) + min(i,j)
c
c     --- diagonalize the second derivative matrix.
      ok=.true.
      call vmove(a,fc,nvv)
      call rzero(xxx,5*nvar)
      call rzero(eigen,nvar)
      call rzero(vec,nvar*nvar)
      call givens(nvar,nvar,nvar,a,xxx,eigen,vec)
c
      if(prnt) then
         write(iout,1010)
         call matprt(vec,nvar,nvar,nvar,nvar,2,0,vname,vname,
     $               0,eigen,.true.)
      else
         write(iout,1040) (eigen(k),k=1,nvar)
         if(neg.ne.0) then
            write(iout,1020)
            do 10 i=1,neg
               call matprt(vec(1,i),nvar,1,nvar,1,1,0,vname,vname,0,
     $                     eigen,.false.)
   10       continue
         endif
      endif
c
c     --- test eigenvalues of 2nd derivative matrix.
      nnege=0
      do 30 i=1,nvar
         ii=i
         rx=eigen(i)
         if(rx.lt.zero) nnege=nnege+1
         if(abs(eigen(i)).lt.eigmin) then
            eigen(i)=sign(eigmin,eigen(i))
         else if(abs(eigen(i)).gt.eigmax) then
            eigen(i)=sign(eigmax,eigen(i))
         endif
         if(((eigen(i).gt.zero).and.(i.le.neg)).or.((eigen(i).lt.zero)
     $        .and.(i.gt.neg))) eigen(i)=-eigen(i)
         if(eigen(i).ne.rx) then
            write(iout,1030) i,eigen(i)
            call matprt(vec(1,i),nvar,1,nvar,1,1,0,vname,vname,0,eigen,
     $                  .false.)
         endif
   30 continue
c
c     --- calculate displacements.
      do 70 i=1,nvar
         do 50 j=i,nvar
            ij=i+(j-1)*nvar
            ji=j+(i-1)*nvar
            rx=zero
            ddx=zero
            do 40 k=1,nvar
               rx=rx+vec(i,k)*vec(j,k)/eigen(k)
               ddx=ddx+vec(i,k)*vec(j,k)*eigen(k)
   40       continue
            fc(lind(i,j))=ddx
            xxx(j)=rx
   50    continue
         do 60 j=1,nvar
            if(j.lt.i) vec(i,j)=vec(j,i)
            if(j.ge.i) vec(i,j)=xxx(j)
   60    continue
   70 continue
c
c     --- generate new coordinates
      dxrms=zero
      dxmax=zero
      do 90 i=1,nvar
         dx=zero
         do 80 j=1,nvar
            dx=dx+vec(j,i)*ftemp(j)
   80    continue
         xnew(i)=xnew(i)+dx
         dxrms=dxrms+dx**2
         if(abs(dx).gt.dxmax) dxmax=abs(dx)
   90 continue
      dxrms=sqrt(dxrms/float(nvar))
c
c
      return
      end

*deck @(#)bigdav.f	1.2  2/27/91
      subroutine bigdav(c,b,smlham,diag,temp,eval,evec,buf,ipt,
     $     scr,thrv,thrn,thrs,
     $     n,nmax,nroot,nsml,nlook,msml,lenb)
c
c***begin prologue     bigdav
c***date written       871022   (yymmdd)
c***revision date      910225   (yymmdd)
c
c   25 february 1991   rlm at lanl
c      removing improper parameter nmx from subroutine statement. 
c   30 january 1988    bhl at brl
c      each iteration nroot new-vectors are created instead of nvec
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)bigdav.f	1.2   2/27/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       bigdav
c
c
      implicit real*8(a-h,o-z)
c
      dimension c(2),b(2),temp(2),diag(2),eval(2),evec(2),buf(2),
     $     ipt(2),scr(2),thrv(20),conv(20),smlham(2)
c
      common /io/ inp,iout
c
ccccc
c     setup programs
ccccc
c
      call bigues(buf,diag,smlham,ipt,n,nsml,nlook,lenb,n)
c
      call movham(smlham,temp,nsml)
c
      call givens(nsml,msml,nsml,temp,scr,eval,evec)
c
      call expand(evec,c,ipt,n,nsml,msml)
      call bighvc(buf,c,b,n,msml,lenb)
c.cc
c      call fasthv(buf,evec,b,n,nsml,msml,ipt)
c.cc
      iter=1
      lx=1
      nvec=msml
      mdim=nvec
c
      go to 110
ccccc
c     loop origin
ccccc
 100  continue
c
      iter=iter+1
      nvec=nroot
      mdim=mdim+nroot
cc
c     multiply matrix times vector
cc
      call bighvc(buf,c(lx),b(lx),n,nroot,lenb)
c
 110  continue
c-----
c     make the small matrix
c-----
      call maksml(smlham,c,b,n,mdim,nvec)
c-----
c     move the small matrix to scratch
c-----
      call movham(smlham,temp,mdim)
c-----
c     diagonalize the small matrix
c-----
      call givens(mdim,nroot,mdim,temp,scr,eval,evec)
c-----
c     form the next vector in the davidson sequence
c-----
      call newvec(diag,eval,b,c,evec,n,nroot,mdim,conv)
c
      konv=0
      do 145 ik=1,nroot
cc    write(iout,150) iter,ik,conv(ik),eval(ik)
         if(conv(ik).lt.thrv(ik))konv=konv+1
 145  continue
c
c 150 format(' iteration ',i5,'  root ',i5,
c    1'   conv ',d10.3,' eigenvalue ',f16.10)
c
      if(konv.eq.nroot)go to 200
c
      lx=lx+nvec*n
c-----
c     orthogonalize the vectors in the davidson sequence
c-----
      call schmdt(c,n,mdim,nroot,thrn,thrs,nroot)
c
      if(iter.lt.nmax)go to 100
ccccc
c     loop termination
ccccc
 200  continue
c-----
c     construct the converged eigenvectors
c-----
      call makvec(c,b,evec,n,mdim,nroot)
c
c     call eigout(b,eval,n,nroot)
c
      return
      end

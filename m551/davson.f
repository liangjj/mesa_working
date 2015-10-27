*deck @(#)davson.f	1.1  11/30/90
      subroutine davson(c,b,smlham,diag,temp,eval,evec,buf,ipt,
     $     scr,thrv,thrn,thrs,
     $     n,nmax,nroot,nsml,nlook,msml)
c
c***begin prologue     davson
c***date written       871022   (yymmdd)
c***revision date      880130   (yymmdd)
c   30 january 1988    bhl at brl
c      each iteration nroot new-vectors are created instead of nvec
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)davson.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       davson
c
      implicit real*8(a-h,o-z)
c
      dimension c(2),b(2),temp(2),diag(2),eval(2),evec(2),buf(2),
     $     ipt(2),scr(2),thrv(20),conv(20),smlham(2)
c
      common /io/ inp,iout
c
      call makdia(buf,diag,n)
c
      call guess(buf,diag,smlham,ipt,n,nsml,nlook)
c
      call movham(smlham,temp,nsml)
c
      call givens(nsml,msml,nsml,temp,scr,eval,evec)
c
      call expand(evec,c,ipt,n,nsml,msml)
c
      call hvec(buf,c,b,n,msml)
c
      iter=1
      lx=1
      nvec=msml
      mdim=nvec
c
      go to 110
c
c     loop origin
c
 100  continue
c
      iter=iter+1
c.
      nvec=nroot
      mdim=mdim+nroot
c.
c
c     multiply matrix times vector
c
      call hvec(buf,c(lx),b(lx),n,nroot)
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
         if(conv(ik).lt.thrv(ik))konv=konv+1
 145  continue
c
      if(konv.eq.nroot)go to 200
c
      lx=lx+nvec*n
c-----
c     orthogonalize the vectors in the davidson sequence
c-----
c.
      call schmdt(c,n,mdim,nroot,thrn,thrs,nroot)
c.
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
c
      return
      end

*deck @(#)sqzopn.f	5.2 4/18/95
      subroutine sqzopn(ng,ngtmp,mxgbsiz,nbf,phi,grad,
     $                  dengrida,dengrada,dengridb,dengradb,
     $                  grdwt,tmpgwt,nzptrs,dencut,becke,lyp)
c***begin prologue     sqzopn.f
c***date written       950229  
c***revision date      4/18/95      
c   february 29, 1995  rlm at lanl
c      segregating this operation in kmtxguts into a subroutine.
c
c***keywords           
c***author             russo, t.v. (lanl) 
c***source             @(#)sqzopn.f	5.2   4/18/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       sqzopn.f
      implicit none
c     --- input variables -----
      integer ng,mxgbsiz,nbf
      real*8 dencut
      logical becke,lyp
c     --- input arrays (unmodified) ---
      real*8 grdwt(mxgbsiz)
c     --- input arrays (scratch) ---
      integer nzptrs(mxgbsiz)
c     --- output arrays ---
      real*8 phi(mxgbsiz,nbf),grad(mxgbsiz,nbf,3)
      real*8 dengrida(mxgbsiz),dengrada(mxgbsiz,3)
      real*8 dengridb(mxgbsiz),dengradb(mxgbsiz,3)
      real*8 tmpgwt(mxgbsiz)
c     --- output variables ---
      integer ngtmp
c     --- scratch arrays ---
c     --- local variables ---
      integer mu
      integer inp,iout
c
      common/io/inp,iout
c
c     --- now we clean the density, and then use nzident 
c         to squeeze out all those values which are zero.
c         We also need to squeeze the grad.  Once the density's done,
c         squeeze out the grid weights corresponding to 
c         zero density points.
c
c     ---  clean based on total density
      call vclean(dengrida,dencut,ng)
      call vclean(dengridb,dencut,ng)
c     form total (alpha+beta) density using tmpgwt as a temp store area
      call vadd(tmpgwt,dengrida,dengridb,ng)
      ngtmp=ng 
      call nzident(tmpgwt,nzptrs,ngtmp)
      call gthr(dengrida,dengrida,nzptrs,ngtmp)
      call gthr(dengridb,dengridb,nzptrs,ngtmp)
      call gthr(tmpgwt,grdwt,nzptrs,ngtmp)
      if (becke.or.lyp)then
         call gthr(dengrada,dengrada,nzptrs,ngtmp)
         call gthr(dengrada(1,2),dengrada(1,2),nzptrs,ngtmp)
         call gthr(dengrada(1,3),dengrada(1,3),nzptrs,ngtmp)
         call gthr(dengradb,dengradb,nzptrs,ngtmp)
         call gthr(dengradb(1,2),dengradb(1,2),nzptrs,ngtmp)
         call gthr(dengradb(1,3),dengradb(1,3),nzptrs,ngtmp)
      endif
      do 10 mu=1,nbf
         call gthr(phi(1,mu),phi(1,mu),nzptrs,ngtmp)
         if (becke.or.lyp)then
            call gthr(grad(1,mu,1),grad(1,mu,1),nzptrs,ngtmp)
            call gthr(grad(1,mu,2),grad(1,mu,2),nzptrs,ngtmp)
            call gthr(grad(1,mu,3),grad(1,mu,3),nzptrs,ngtmp)
         endif
 10   continue 
c
c
      return
      end

*deck %W%  %G%
      subroutine md2ghi(g,h,i2,f00,b00,b10,bp01,c00,cp00,nmax,mmax,
     $                  imax,jmax,kmax,lmax,c,nat,symcen,nv,dg,dh,di,
     $                  dc00,dcp00,ndcen,nd2,d1exp,d2exp,npass,d2g,
     $                  d2h,d2i)
c***begin prologue     md2ghi.f
c***date written       851113  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl) 
c***source             %W%   %G%
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       md2ghi.f
      implicit none
c     --- input variables -----
      integer nmax,mmax,imax,jmax,kmax,lmax,nat,nv
      integer nd2,ndcen,npass
c     --- input arrays (unmodified) ---
      integer symcen(4)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 g(nv,3,0:nmax,0:mmax)
      real*8 h(nv,3,0:nmax,0:jmax,0:mmax)
      real*8 i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 f00(nv)
      real*8 b00(nv)
      real*8 b10(nv)
      real*8 bp01(nv)
      real*8 c00(nv,3)
      real*8 cp00(nv,3)
      real*8 dg(nv,3,0:nmax,0:mmax,ndcen)
      real*8 dh(nv,3,0:nmax,0:jmax,0:mmax,ndcen)
      real*8 di(nv,3,0:imax,0:jmax,0:mmax,0:lmax,ndcen)
      real*8 d2g(nv,3,0:nmax,0:mmax,nd2)
      real*8 d2h(nv,3,0:nmax,0:jmax,0:mmax,nd2)
      real*8 d2i(nv,3,0:imax,0:jmax,0:mmax,0:lmax,nd2)
      real*8 dc00(nv,ndcen)
      real*8 dcp00(nv,ndcen)
      real*8 d1exp(nv,3,ndcen)
      real*8 d2exp(nv,nd2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer cen(2,6,4)
      integer i,j
c
c     --- 'cen' gives the two derivatives for second-derivatives for each
c         integral type
      data ((cen(i,j,1),i=1,2),j=1,6) /1,1, 1,2, 2,2, 1,3, 2,3, 3,3/
      data ((cen(i,j,2),i=1,2),j=1,3) /1,1, 1,2, 2,2/
      data ((cen(i,j,3),i=1,2),j=1,3) /1,1, 1,2, 2,2/
      data  (cen(i,1,4),i=1,2)        /1,1/
      save cen
c
c
      call vmakg(g,f00,b00,b10,bp01,c00,cp00,nmax,mmax,nv,nv)
      call makdg(g,dg,f00,b00,b10,bp01,c00,cp00,dc00,dcp00,nmax,
     $           mmax,nv,nv,ndcen)
      call makd2g(g,dg,d2g,b00,b10,bp01,c00,cp00,dc00,dcp00,nmax,
     $            mmax,nv,nv,ndcen,nd2,cen(1,1,npass))
c
      call vmakh(h,g,c,nat,symcen,nmax,mmax,imax,jmax,nv)
      call makdh(h,dh,dg,c,nat,symcen,nmax,mmax,imax,jmax,nv,ndcen,
     $           npass)
      call makd2h(dh,d2h,d2g,c,nat,symcen,nmax,mmax,imax,jmax,nv,ndcen,
     $            nd2,npass)
c
      call vmaki2(i2,h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,nv)
      call makdi(i2,di,dh,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,
     $           nv,ndcen,npass)
      call makd2i(di,d2i,d2h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,
     $            lmax,nv,ndcen,nd2,npass)
c
      call d2expn(i2,di,d2i,imax,jmax,kmax,lmax,mmax,nv,ndcen,nd2,
     $            d1exp,d2exp,cen(1,1,npass))
      call derexp(i2,di,imax,jmax,kmax,lmax,mmax,nv,ndcen,d1exp)
c
c
      return
      end

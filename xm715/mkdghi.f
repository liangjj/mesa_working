*deck %W%  %G%
      subroutine mkdghi(g,h,i2,f00,b00,b10,bp01,c00,cp00,nmax,mmax,
     $                  imax,jmax,kmax,lmax,c,nat,symcen,nv,
     $                  dg,dh,di,dc00,dcp00,ndcen,d1exp,npass)
c
      real*8 g(nv,3,0:nmax,0:mmax),h(nv,3,0:nmax,0:jmax,0:mmax)
      real*8 i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 f00(nv),b00(nv),b10(nv),bp01(nv),c00(nv,3),cp00(nv,3)
      real*8 c(3,nat)
      real*8 dg(nv,3,0:nmax,0:mmax,ndcen)
      real*8 dh(nv,3,0:nmax,0:jmax,0:mmax,ndcen)
      real*8 di(nv,3,0:imax,0:jmax,0:mmax,0:lmax,ndcen)
      real*8 dc00(nv,ndcen),dcp00(nv,ndcen),d1exp(nv,ndcen)
      integer nmax,mmax,imax,jmax,kmax,lmax,nat,symcen(4),nv,ndcen
c
c     ----- timing -----
c
c
      call vmakg(g,f00,b00,b10,bp01,c00,cp00,nmax,mmax,nv,nv)
c
      call makdg(g,dg,f00,b00,b10,bp01,c00,cp00,dc00,dcp00,nmax,
     $           mmax,nv,nv,ndcen)
c
      call vmakh(h,g,c,nat,symcen,nmax,mmax,imax,jmax,nv)
c
      call makdh(h,dh,dg,c,nat,symcen,nmax,mmax,imax,jmax,nv,ndcen,
     $           npass)
c
      call vmaki2(i2,h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,nv)
c
      call makdi(i2,di,dh,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,
     $           nv,ndcen,npass)
c
      call derexp(i2,di,imax,jmax,kmax,lmax,mmax,nv,ndcen,d1exp)
c
c     ----- timing -----
c
c
c
      return
      end

*deck @(#)vmkghi.f	5.1  11/6/94
      subroutine vmkghi(g,h,i2,f00,b00,b10,bp01,c00,cp00,nmax,mmax,
     #                  imax,jmax,kmax,lmax,c,nat,symcen,nv)
c
      real*8 g(nv,0:nmax,0:mmax,3),h(nv,0:nmax,0:jmax,0:mmax,3)
      real*8 i2(nv,0:imax,0:jmax,0:mmax,0:lmax,3)
      real*8 f00(nv),b00(nv),b10(nv),bp01(nv),c00(nv,3),cp00(nv,3)
      real*8 c(3,nat)
      integer nmax,mmax,imax,jmax,kmax,lmax,nat,symcen(4),nv
c
c     ----- timing -----
c
c
      call vmakg(g,f00,b00,b10,bp01,c00,cp00,nmax,mmax,nv,nv)
      call vmakh(h,g,c,nat,symcen,nmax,mmax,imax,jmax,nv)
      call vmaki2(i2,h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,nv)
c
c     ----- timing -----
c
c
c
      return
      end

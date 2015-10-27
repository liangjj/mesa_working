*deck @(#)vmkghi.f	5.1  11/6/94
      subroutine vmkghi(g,h,i2,f00,b00,b10,bp01,c00,cp00,nmax,mmax,
     #                  imax,jmax,kmax,lmax,c,nat,symcen,nv,lenv)
c
      real*8 g(nv,3,0:nmax,0:mmax),h(nv,3,0:nmax,0:jmax,0:mmax)
      real*8 i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 f00(nv),b00(nv),b10(nv),bp01(nv),c00(lenv,3),cp00(lenv,3)
      real*8 c(3,nat)
      integer nmax,mmax,imax,jmax,kmax,lmax,nat,symcen(4),nv
c
c
      call vmakg(g,f00,b00,b10,bp01,c00,cp00,nmax,mmax,nv,lenv)
      call vmakh(h,g,c,nat,symcen,nmax,mmax,imax,jmax,nv)
      call vmaki2(i2,h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,lmax,nv)
c
c
      return
      end

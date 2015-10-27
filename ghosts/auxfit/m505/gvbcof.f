*deck @(#)gvbcof.f	3.1  11/20/92
      subroutine gvbcof(jmat,kmat,h,nnp,ncoul,nexch,y,nnv,shlmin,
     $                  shlmax,nshell,nv,first,c,root,eigval,
     $                  u,t1,t2)
c
c***begin prologue     gvbcof
c***date written       870524   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           gvb coefficients, ci coefficients
c***author             saxe, paul (lanl)
c***source             @(#)gvbcof.f	3.1   11/20/92
c
c***purpose
c
c***description        to calculate the ci coefficients for a gvb
c                       wavefunction.
c***references
c
c***routines called    (none)
c
c***end prologue       gvbcof
c
      implicit integer (a-z)
c
      real*8 jmat(nnp,ncoul),kmat(nnp,nexch),h(nnp),y(nnv),c(nv)
      real*8 eigval(nv),u(nv,nv),t1(nv),t2(nv)
      integer shlmin(nshell),shlmax(nshell)
      real*8 av,yv
      real*8 two,half
      real*8 f
c
      parameter (two=2.0d+00, half=0.5d+00)
c
c     ----- statement functions -----
c
      f(ij)=h(ij)+two*jmat(ij,1)-kmat(ij,1)
c
c     ----- form diagonals to y matrix -----
c
c       y(i,i)=hc(i,i)+1/2j(i,i)+ f(i)*(2j(k,i)-k(k,i))
c                where k not in this correlated pair
c
      do 5 i=1,nv
         iorb=first+i-1
         ia=iorb*(iorb-1)/2
         ii=ia+iorb
         iav=i*(i-1)/2
         y(iav+i)=f(ii)+half*jmat(ii,1+i)
c
c        ----- and pop in off-diagonals 1/2k(i,j) -----
c
         do 2 j=1,i-1
            jorb=first+j-1
            ij=ia+jorb
            y(iav+j)=half*kmat(ii,1+j)
    2    continue
    5 continue
c
c      call matout(y,nv,nv,nv,nv,iout)
c
c     ----- diagonalize this matrix -----
c
      call rsp(nv,nv,nnv,y,eigval,1,u,t1,t2,error)
c
      if (error.ne.0) call lnkerr('error diagonalizing gvb '//
     $                            'matrix')
c
c     ----- temporarily use two orbital formula -----
c
c      yv=y(2,2)-y(1,1)
c      av=yv/(two*y(1,2))-sqrt(1.0d+00+(yv/(two*y(1,2)))**2)
c
c      c1=1.0d+00/sqrt(1.0d+00+av**2)
c      c2=av/sqrt(1.0d+00+av**2)
c
c
c     ----- transfer the desired roots vector to coefficents array ----
c
      do 10 i=1,nv
         c(i)=u(i,root)
   10 continue
c
c      write (iout,100) c
c  100 format (' c: ',(10f10.6))
c
      return
      end

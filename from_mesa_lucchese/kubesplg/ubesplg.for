      program ubesplg
c
c bessel fuction spline code with new prescription for
c the outgoing wave continuum function.
c
      implicit real*8 (a-h,o-z)
      real*8 mu
      parameter (nbig=15000,lbig=#maxltop,nsmall=4015)
      parameter (nchmx=#maxchan,nemx=200)
      parameter (nsmall2=2*nsmall)
      dimension ipow(0:lbig),echan(nchmx),energy(nemx)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
      common/spl/rd,r(nsmall),w(nsmall,0:lbig),cs(nsmall),rd26,l
      complex*16 csc(nsmall,0:lbig),scc(nsmall)
      real*8 aj(0:(lbig+1)),ay(0:(lbig+1))
      real*8 aj1(0:(lbig+1)),ay1(0:(lbig+1))
      real*8 csl(nsmall,0:lbig),derl(nsmall,0:lbig)
      complex*16 x1,x2,f1,f2,f,fp,gp,sig(0:lbig),tv
      complex*16 yc1,yc2,zl
      real*8 css(nsmall)
      complex*16 c1,c2,ai,ww(nsmall),wwl(nsmall,0:lbig)
      external zero,gg
      dimension der(nsmall),xx(nbig),yy(nbig),zz(nbig)
      dimension scr(nsmall)
      dimension wwlr(nsmall2,0:lbig),cscr(nsmall2,0:lbig)
      equivalence (wwl(1,0),wwlr(1,0)),(csc(1,0),cscr(1,0))
      ai=cmplx(0.,1.)
c
c..unicos
c      call link ("unit5=infreeg,unit6=(outfree,hc,create)
c..unicos     1,unit8=bessplr//")
c
      open(5,file='infree')
      open(6,file='outfree')
      open(8,file='bessplr',form='unformatted')
c
      ipow(0)=0
      do 66 i=1,lbig
   66 ipow(i)=1
c     call keep80(4hwave)
c     call frid(10hxerox+film,1,4,1)
c     call dders(-1)
      eta=0.
      znuc=0.
      mu=1.
c
c lmax=max. l value desired
c nper=no. of points per unit a.u. for splined functions
c xstart = first r-point (.00001 is o.k.)
c xmax = last r-point(r=k*x!)
c nint = multiple of nper for integration(3 or 4 is reasonable)
c alpha = exponential parameter in test function
c
      read(5,*)lmax,nper,xstart,xmax,nint,alpha
      write(6,*)' Greens Function Code '
      write(6,101)lmax,nper,xstart,xmax,nint,alpha
 101  format(/,' L-Max  = ',i5,/,
     2         ' Nper   = ',i5,/,
     3         ' XStart = ',f20.12,/,
     4         ' XMax   = ',f20.12,/,
     5         ' Nint   = ',i5,/,
     6         ' Alpha  = ',f20.12)
      np=nper*xmax+1
      if(np.gt.nsmall)then
      write(6,100)
100   format(" np.gt.nsmall")
      stop
      endif
      rd=(xmax-xstart)/(np-1)
      rd26=rd*rd/6.
      do 1 i=1,np
1     r(i)=rd*(i-1.)+xstart
      aint=nint
      h=rd/aint
      nl=lmax+1
      do 668 i=1,np
      x=r(i)
      call sjymec(x,aj,ay,ncomp,lmax)
      do 668 j=0,lmax
668   w(i,j)=exp(-x*alpha)*aj(j)*x
c668   w(i,j)=exp(-x*x*alpha)*x
c668   w(i,j)=exp(-x*x*alpha)*x**(1+j)
667   continue
      do 33 l=0,lmax
      al=l
      ystart=xstart**(l+1)
      ypstart=(al+1.)*xstart**l
      yp1=(w(2,l)-w(1,l))/rd
      yp2=(w(np,l)-w(np-1,l))/rd
      call spline(r,w(1,l),np,yp1,yp2,scr,cs)
      call outward(xx,yy,index,zero)
      call outward(xx,zz,index,gg)
      i1=index
      i2=index-1./h
      x1=xx(i1)
      x2=xx(i2)
      det=yy(i1)*zz(i2)-yy(i2)*zz(i1)
      call sjymec(x1,aj,ay,ncomp,l)
      call sjymec(x2,aj1,ay1,ncomp,l)
      f1=x1*(-ay(l)+ai*aj(l))
      f2=x2*(-ay1(l)+ai*aj1(l))
      c1=f1*zz(i2)-f2*zz(i1)
      c1=c1/det
      c2=f2*yy(i1)-f1*yy(i2)
      c2=c2/det
      j=0
      do 3 i=1,index,nint
      j=j+1
      ww(j)=c1*yy(i)+c2*zz(i)
      ww(j)=ww(j)/xx(i)
      wwl(j,l)=ww(j)
3     continue
      yc1=(ww(2)-ww(1))/rd
      yc2=(ww(np)-ww(np-1))/rd
      call splinec(r,ww,np,yc1,yc2,scc,csc(1,l))
      c2r=c2
      do 5 i=1,np
      der(i)=w(i,l)*c2r
      der(i)=der(i)/r(i)**ipow(l)
      derl(i,l)=der(i)
5     continue
      yp1=(der(2)-der(1))/rd
      yp2=(der(np)-der(np-1))/rd
      call spline(r,der,np,yp1,yp2,scr,csl(1,l))
33    continue
      write(8)lmax,np,xstart,rd,alpha
      write(8)(r(i),i=1,np)
      np2=2*np
c      write(6,*)"wwl"
c      write(6,"(6f12.5)")((wwlr(i,k),i=1,np2),k=0,lmax)
c      write(6,*)"derl"
c      write(6,"(6f12.5)")((derl(i,k),i=1,np),k=0,lmax)
c      write(6,*)"csc"
c      write(6,"(6f12.5)")((cscr(i,k),i=1,np2),k=0,lmax)
c      write(6,*)"csl"
c      write(6,"(6f12.5)")((csl(i,k),i=1,np),k=0,lmax)
      write(8)((wwlr(i,k),i=1,np2),k=0,lmax)
      write(8)((derl(i,k),i=1,np),k=0,lmax)
      write(8)((cscr(i,k),i=1,np2),k=0,lmax)
      write(8)((csl(i,k),i=1,np),k=0,lmax)
      write(6,*)' BessplG Finished '
      stop
      end

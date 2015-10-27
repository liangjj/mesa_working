      subroutine radial(lmax,np,r,wt,ak,akp,aplus,aminus)
      implicit real*8(a-h,o-z)
      dimension aplus(1),aminus(1),r(1),wt(1)
      dimension aj(0:100),ajp(0:100),ay(0:100)
      do 1 i=1,np
      x1=ak*r(i)
      x2=akp*r(i)
      call sjymec(x2,ajp,ay,ncal,lmax)
      call sjymec(x1,aj,ay,ncal,lmax)
      do 1 l=1,lmax
      aplus(l)=aplus(l)+aj(l-1)*ajp(l)*wt(i)
      aminus(l)=aminus(l)+aj(l)*ajp(l-1)*wt(i)
1     continue
c     write(6,66)ak,akp
66    format(" ak,akp:",2e14.4)
c     do 12 i=1,lmax
c     write(6,100)i,aplus(i),aminus(i)
12    continue
100   format(i5,2e12.4)
      return
      end

*deck %W%  %G%
      subroutine quadr(lalo,lahi,lblo,lbhi,nlo,nhi,qq)
      implicit integer(a-z)
c
c     ----- arguments unchanged -----
      integer lalo,lahi,lblo,lbhi,nlo,nhi
      real*8 qq(0:6,0:6,0:6)
c     _____ local variables -----
      integer la,lb,n,i
c     ----- common -----
      real*8 ptpow,f,pt
      integer npts
      common/ptwtdat/ptpow(50,0:6),f(50,0:6,0:6),pt(50)
      common/ptwt/npts
c
c     ----- evaluate the integral by quadrature -----
      do 40 lb=lblo,lbhi
         do 30 la=lalo,lahi
            do 20 n=nlo,nhi
               do 10 i=1,npts
                  qq(n,la,lb)=qq(n,la,lb)+ptpow(i,n)*f(i,la,lb)
   10          continue
   20       continue
   30    continue
   40 continue
c
      return
      end

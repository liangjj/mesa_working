*deck @(#)quadr.f	5.1  11/6/94
      subroutine quadr(lmalo,lmahi,lmblo,lmbhi,nlo,nhi,qq)
      implicit real*8(a-h,o-z)
      common/ptwtdat/ptpow(50,11),f(50,9,9),pt(50)
      common/ptwt/npts
      dimension qq(11,9,9)
      common/io/inp,iout
c
c     this routine evaluates the integral by quadrature.
c
      do 40 lb=lmblo,lmbhi
         do 30 la=lmalo,lmahi
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

*deck rmamp
c***begin prologue     rmamp
c***date written       000623   (yymmdd)
c***revision date               (yymmdd)
c***keywords           r-matrix amplitudes for DVR basis
c***author             schneider, barry (nsf)
c***source             
c***purpose            construct r-matrix amplitudes
c***description        
c***references       
c
c***routines called
c***end prologue       rmamp
      subroutine rmamp(eigvec,p,rgamma,n,prn)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 eigvec, p, rgamma
      character*80 title
      logical prn
      dimension eigvec(n,n), p(2), rgamma(n,2)
      do 10 i=1,n
         rgamma(i,1) = eigvec(1,i)*p(1)
         rgamma(i,2) = eigvec(n,i)*p(2)
 10   continue
      if(prn) then
        title='r-matrix amplitudes'
        call prntrm(title,rgamma,n,2,n,2,iout)
      endif
      return
      end

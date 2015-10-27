*deck gam1.f
c***begin prologue     gam1
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            R-matrix amplitudes.
c***                   
c***references         
c
c***routines called    
c***end prologue       gam1
      subroutine gam1(pham,n,prn)
      implicit integer (a-z)
      integer*8 pham
      real*8 ham, eigc
      character*80 title
      logical prn
      common/io/inp, iout
      pointer (pham,ham(1))
      hmat=1
      v=hmat+n*n
      h0=v+n
      srf=h0+n*n      
      pt0=srf+2
      q1=pt0+1
      qwt=q1+n
      pq=qwt+n
      dpq=pq+n*n
      ddpq=dpq+n*n
      eigv=ddpq+n*n
      rgama=eigv+n*n+n
      one=1
      eigc=0.d0
      call iosys('write integer "number of channels" to ham',1,
     1            one,0,' ')
      call iosys('write real "channel energies" to ham',1,
     1            eigc,0,' ')
      call iosys('write real "h amplitudes" to ham',n,
     1            ham(rgama),0,' ')
      if(prn) then
         title='r-matrix amplitudes'
         call prntrm(title,ham(rgama),n,1,n,1,iout)
      endif
      return      
      end       







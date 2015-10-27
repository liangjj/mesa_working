*deck ov1d.f
c***begin prologue     ov1d
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            projection of perturbed onto unperturbed states.
c***                   
c***references         
c
c***routines called    
c***end prologue       ov1d
      subroutine ov1d(ov,vec,vec0,n,prn)
      implicit integer (a-z)
      real*8 ov, vec, vec0
      character*80 title
      logical prn
      dimension ov(n,n), vec(n,n), vec0(n,n)
      common/io/inp, iout
      fac=1.d0/sqrt(2.d0)
      call ebtc(ov,vec0,vec,n,n,n)
      if(prn) then
         title='< psi0(i) | psi(q) >'
         call prntfm(title,ov,n,n,n,n,iout)
      endif          
      return      
      end       







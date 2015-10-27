*deck fevodd.f
c***begin prologue     fevodd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            make polynomials of even and odd parity
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       fevodd
      subroutine fevodd(p,dp,ddp,pn,dpn,ddpn,neven,nodd,type,n)
      implicit integer (a-z)
      real*8 p, dp, ddp, pn, dpn, ddpn, fac
      character*80 title
      dimension p(n,*), dp(n,*), ddp(n,*)
      dimension pn(n,*), dpn(n,*), ddpn(n,*)
      common/io/inp, iout
      write(iout,1)
      fac=1.d0/sqrt(2.d0)
c
c     determine if n is even or odd.
c
      ntst=mod(n,2)
      if(ntst.eq.0) then
         neven=n/2
         nodd=neven
         nback=n
         do 10 i=1,neven
            do 20 j=1,n
               pn(j,i)     =     fac * ( p(j,i) + p(j,nback) )
               pn(j,nback) = fac * ( p(j,i) - p(j,nback) )
               dpn(j,i)     =     fac * ( dp(j,i) + dp(j,nback) )
               dpn(j,nback) = fac * ( dp(j,i) - dp(j,nback) )
               ddpn(j,i)     =     fac * ( ddp(j,i) + ddp(j,nback) )
               ddpn(j,nback) = fac * ( ddp(j,i) - ddp(j,nback) )
 20         continue 
            nback = nback - 1   
 10      continue
      else
         nodd=n/2
         neven=nodd+1
         nback = n
         do 30 i=1,nodd
            do 40 j=1,n
               pn(j,i)     =     fac * ( p(j,i) + p(j,nback) )
               pn(j,nback) = fac * ( p(j,i) - p(j,nback) )
               dpn(j,i)     =     fac * ( dp(j,i) + dp(j,nback) )
               dpn(j,nback) = fac * ( dp(j,i) - dp(j,nback) )
               ddpn(j,i)     =     fac * ( ddp(j,i) + ddp(j,nback) )
               ddpn(j,nback) = fac * ( ddp(j,i) - ddp(j,nback) )
 40         continue
            nback = nback - 1   
 30      continue
         call copy(p(1,nback),pn(1,nback),n)
      endif   
      write(iout,2) neven, nodd
c      title='good parity polynomials'
c      call prntfm(title,pn,n,n,n,n,iout)
c      title='first derivative of good parity polynomials'
c      call prntfm(title,dpn,n,n,n,n,iout)
c      title='second derivative of good parity polynomials'
c      call prntfm(title,ddpn,n,n,n,n,iout)
c      call lnkerr('quit')
      return
 1    format(/,5x,'forming solutions of good parity')
 2    format(/,5x,'number of even parity functions = ',i4,/,5x,
     1            'number of odd parity functions  = ',i4)
      end       

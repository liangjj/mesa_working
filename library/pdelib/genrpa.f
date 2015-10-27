*deck genrpa.f
c***begin prologue     genrpa
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c
c***routines called    
c***end prologue       genrpa
      subroutine genrpa(eps,kappa,q,n)
      implicit integer (a-z)
      real*8 fpkey, eps, kappa, q, sdot, norm
      real*8 gasdev
      character*80 cpass
      character*1600 card
      logical dollar
      dimension q(*) 
      common/io/inp, iout
      data iseed / 2777 /
      if ( dollar('$genrpa',card,cpass,inp) ) then
         eps=fpkey(card,'particle-hole-spacing',.1d0,' ')
         kappa=fpkey(card,'strength-of-coupling',10.d0,' ')
      endif
      write(iout,1) n, eps, kappa
c
      q(1)= (n - 1)*gasdev(iseed)     
      do 10 i=2,n
         q(i) = q(i-1) + n - i - i + 1
         q(i)=q(i)*gasdev(iseed)
 10   continue
      norm = 1.d0/sqrt(sdot(n,q,1,q,1))
      do 20 i=1,n
         q(i) = norm*q(i)
 20   continue                   
 1    format(/,1x,'rpa matrix of dimension = ',i3,/,1x,
     1            'particle-hole-spacing   = ',e15.8,/,1x,
     2            'strength-of-coupling    = ',e15.8)                 
      return
      end       

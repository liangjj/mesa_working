c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Moment}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck moment.f
***begin prologue     moment
***date written       020127   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           time development moment
***author             schneider, barry (nsf)
***source             
***purpose            calculate some moments and compare.
***references         
***routines called    
***end prologue       moment
      subroutine moment(chi,psi0,t,q1,q2,q3,wt1,wt2,wt3,n,nd,
     1                  dim,coord,typke,card)
      implicit integer (a-z)
      real*8 t, q1, q2, q3, wt1, wt2, wt3
      real*8 fpkey, sigma, beta, x_0 
      complex*16 chi, psi0
      real*8 avx, avxx, exavx, exavxx
      character*(*) coord, typke, card
      character*2 itoc
      logical prn, dollar, prnton
      character*80 title, chrkey
      dimension nd(3)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      dimension wt1(nd(1)), wt2(nd(2)), wt3(nd(3))
      dimension chi(n), psi0(n), sigma(3), beta(3), x_0(3)
      common/io/inp, iout
      if(dim.ne.1) then
         write(iout,1)
         call lnkerr('quit')
      endif
      if( dollar('$initial-state',card,title,inp) ) then  
          call fparr(card,'sigma',sigma,3,' ')  
          call fparr(card,'x_0',x_0,3,' ')
          call fparr(card,'beta',beta,3,' ')  
      endif  
      do 10 i=1,n
         chi(i) = chi(i) + psi0(i)
 10   continue   
      avx=0.d0
      avxx=0.d0
      if(typke.eq.'fd') then
         do 20 i=1,nd(1)
            avx = avx + wt1(i)*q1(i)*chi(i)*conjg(chi(i))
            avxx = avxx + wt1(i)*q1(i)*q1(i)*chi(i)*conjg(chi(i))
 20      continue
      else
         do 30 i=1,nd(1)
            avx = avx   + q1(i)*chi(i)*conjg(chi(i))
            avxx = avxx + q1(i)*q1(i)*chi(i)*conjg(chi(i))
 30      continue
      endif
      avxx = avxx - avx*avx    
      avxx=sqrt(avxx)
      exavx=beta(1)*t
      exavxx= .5d0*( sigma(1)*sigma(1) + t*t/( sigma(1)*sigma(1)) )
      exavxx=sqrt(exavxx)
      write(iout,2) exavx, avx, exavxx, avxx
       do 40 i=1,n
         chi(i) = chi(i) - psi0(i)
 40   continue   
      return
 1    format(/,1x,'This is a test for one-dimensional free wavepacket',
     1       /,1x,'Current calculation does not qualify')
 2    format(/,1x,'Exact <x>                 = ',e15.8,
     1       /,1x,'Calculated <x>            = ',e15.8,
     2       /,1x,'Exact <x*x> - <x><x>      = ',e15.8,
     3       /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
      end       

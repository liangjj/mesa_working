c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Subroutine for Radial Wavepacket}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle
c deck rad_paket.f
c***begin prologue     radpakt
c**date written       960723   (yymmdd)
c**revision date      yymmdd   (yymmdd)
c**keywords           time development
c**author             schneider, barry (nsf)
c**source
c**purpose            calculate zero time rad wavepacket
c**references
c**routines called
c**end prologue       radpakt
      subroutine radpakt(psi0,q,expnt,exppow,typpak,n)
      implicit integer(a-z)
      real*8 q, psi0, norm, expnt, sdot
      character*(*) typpak
      dimension q(n), psi0(n,2)
      common/io/inp, iout
      call rzero(psi0,2*n)
      write(iout,1)
      if(typpak.eq.'exponential') then
         write(iout,2)
         do i=1,n
            psi0(i,1) = q(i)**exppow * exp( -expnt * q(i) )
         end do
      elseif(typpak.eq.'gaussian') then
         write(iout,3)
         do i=1,n
            psi0(i,1) = q(i)**exppow * exp( -expnt * q(i)*q(i) )
         end do
      end if
      WRITE(iout,4)
      WRITE(iout,5) expnt
      WRITE(iout,6) exppow
      norm = 1.d0/sqrt(sdot(n,psi0(1,1),1,psi0(1,1),1))
      do i=1,n
         psi0(i,1) = psi0(i,1) * norm  
      end do
      return
 1    FORMAT(/,5X,'initial wavepacket at t=0')
 2    FORMAT(/,5X,'the form of the radial packet is:',///,5X,
     1            'psi = r**n * exp( - alpha * r)')
 3    FORMAT(/,5X,'the form of the radial packet is:',///,5X,  
     1            'psi = r**n * exp( - alpha * r *r)')
 4    FORMAT(/,1X,'radial wave packet parameters')
 5    FORMAT(/,1X,'alpha    = ',e15.8,1X)
 6    FORMAT(/,1X,'power    = ',i2)
 7    FORMAT(/,2x,'Normalization of Initial Wavepacket = ',
     1             e15.8,1x,e15.8)
      END 








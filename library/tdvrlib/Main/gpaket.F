c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Subroutine for Gaussian Wavepacket}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck gpaket.f
c***begin prologue     gpaket
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate zero time gaussian wavepacket 
c***                   
c***references         
c
c***routines called    
c***end prologue       gpaket
      subroutine gpaket(q1,q2,q3,norm,psi0,sigma,x_0,beta,n,
     1                  nd,dim,coord,prn)
      implicit integer (a-z)
      real*8 q1, q2, q3, norm, psi0
      real*8 sigma, x_0, beta
      real*8 z
      character*(*) coord
      character*80 title
      logical prn
      dimension nd(3)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      dimension sigma(3), x_0(3), beta(3)
      dimension psi0(n,2)
      common/io/inp, iout
      pointer (p,z(1))

c \center Gaussian Wave Packet
c \par
c\begin{eqnarray}
c \Psi({\bf r}) & = N exp( -\frac{ ( x -x_0 )^2 }{2 { {\sigma}_x }^2 }
c                                - i {\beta}_x (x-x_0) ) & \nonumber \\
c                         &exp( -\frac{ ( y -y_0 )^2 }{2 { {\sigma}_y }^2 }
c                               - i {\beta}_y (y-y_0) ) & \nonumber \\
c                         &exp( -\frac{ ( z -z_0 )^2 }{2 { {\sigma}_z }^2 }
c                               - i {\beta}_z (z-z_0) ) &
c \end{eqnarray}     

      nbig=0
      do 10 i=1,dim
         nbig=max(nbig,nd(i))
 10   continue   
      if(dim.ge.1) then
         i1=1
         need=i1+2*nbig
      endif
      if(dim.ge.2) then
         i2=need
         need= i2 + 2*nbig
      endif
      if(dim.eq.3) then
         i3=need
         need= i3 + 2*nbig
      endif  
      need=wpadti(need)
      call getmem(need,p,ngot,'gpaket',0)
      write(iout,1)
      call rzero(psi0,n*2)
      if(coord.eq.'cartesian') then
         write(iout,2)
         write(iout,3)
         write(iout,4) (sigma(i),i=1,dim)
         write(iout,5) (x_0(i),i=1,dim)
         write(iout,6) (beta(i),i=1,dim)
         if(dim.eq.1) then
            call zproj(q1,sigma(1),x_0(1),beta(1),z(i1),nd(1),coord,prn)
            call zfill1(psi0,z(i1),n)
         elseif(dim.eq.2) then
            call zproj(q1,sigma(1),x_0(1),beta(1),z(i1),nd(1),coord,prn)
            call zproj(q2,sigma(2),x_0(2),beta(2),z(i2),nd(2),coord,prn)
            call zfill2(psi0,z(i1),z(i2),nd,n)
         elseif(dim.eq.3) then
            call zproj(q1,sigma(1),x_0(1),beta(1),z(i1),nd(1),coord,prn)
            call zproj(q2,sigma(2),x_0(2),beta(2),z(i2),nd(2),coord,prn)
            call zproj(q3,sigma(3),x_0(3),beta(3),z(i3),nd(3),coord,prn)
            call zfill3(psi0,z(i1),z(i2),z(i3),nd,n)
         endif
      else
         write(iout,7)
         write(iout,8) sigma(1), x_0(1), beta(1)
         call zproj(q1,sigma(1),x_0(1),beta(1),z(i1),nd(1),coord,prn)
         call zfill1(psi0,z(i1),n)
      endif
      call getmem(-ngot,p,idum,'gpaket',idum) 
      return
 1    format(/,5x,'initial wavepacket')
 2    format(/,5x,'the form of the cartesian packet is:',///,5x,
     1            'psi = exp( - ( x - x_0 )**2 / (2*(a_x)**2)'
     2                      '     i*b_x*x )',/,5x,
     3            '      exp( - ( y - y_0 )**2 / (2*(a_y)**2)'
     4                      '     i*b_y*y )',/,5x,
     5            '      exp( - ( z - z_0 )**2 / (2*(a_z)**2)'
     6                      '     i*b_z*z )')
 3    format(/,1x,'gaussian wave packet parameters')
 4    format(/,1x,'sigma    = ',3(e15.8,1x))
 5    format(/,1x,'x_0      = ',3(e15.8,1x))
 6    format(/,1x,'beta     = ',3(e15.8,1x)) 
 7    format(/,5x,'the form of the radial packet is:',///,15x,
     1            'psi = r*exp(- ( r - r_0 )**2 / (2*(a_r)**2)'
     2                        '- i*b_r*r)')
 8    format(/,1x,'gaussian wave packet parameters'
     1        ,/,1x,'sigma = ',e15.8,/,1x,'r_0 = ',e15.8,
     2         /,1x,'beta  = ',e15.8)    
      end       








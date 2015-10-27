c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Reference Weight Function}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck genrwt.f
c***begin prologue     genrwt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           reference weight
c***author             schneider, barry (nsf)
c***source
c***purpose            calculate weight functions and derivatives.
c***                   
c***description             
c***                                                                       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       genrwt
      subroutine genrwt(rwt,drwt,ddrwt,pt,wt,alpha,beta,
     1                  deriv,edge,n)
      implicit integer (a-z)
      real*8 rwt, drwt, ddrwt, pt, alpha, beta, edge
      real*8 x1, x2, fac1, fac2
      character*(*) wt
      logical deriv
      dimension rwt(n), drwt(n), ddrwt(n), pt(n), edge(2)
      common/io/inp, iout 

c     Weight functions and their derivatives are computed for a variety of
c     cases

      if(wt.eq.'legendre') then
         call vfill(rwt,1.d0,n)
         if(deriv) then
            call rzero(drwt,n)
            call rzero(ddrwt,n)
         endif
      elseif(wt.eq.'r') then
         do 10 i=1,n
            rwt(i)=pt(i)
 10      continue   
         if(deriv) then
            call vfill(drwt,1.d0,n)
            call rzero(ddrwt,n)
         endif
      elseif(wt.eq.'rr') then
         do 20 i=1,n
            rwt(i)=pt(i)*pt(i)
 20      continue   
         if(deriv) then
            do 30 i=1,n
               drwt(i)=2.d0*pt(i)
               ddrwt(i)=2.d0
 30         continue   
         endif
      elseif(wt.eq.'hermite') then
         do 40 i=1,n
            rwt(i)=exp(-pt(i)*pt(i))
 40      continue   
         if(deriv) then
            do 50 i=1,n
               drwt(i) = -2.d0*rwt(i)
               ddrwt(i)= ( -2.d0 + 4.d0*pt(i)*pt(i) ) * rwt(i)
 50         continue   
         endif
      elseif(wt.eq.'chebyshev-1') then
         do 60 i=1,n
            rwt(i)= sqrt ( 1.d0/( 1.d0 - pt(i)*pt(i) ) )
 60      continue   
         if(deriv) then
            do 70 i=1,n
               ddrwt(i) = rwt(i)*rwt(i)*rwt(i)
               drwt(i) = pt(i)*ddrwt(i)
               ddrwt(i) = ddrwt(i) 
     1                     + 3.d0*pt(i)*pt(i)*ddrwt(i)*rwt(i)*rwt(i)
 70         continue
         endif   
      elseif(wt.eq.'chebyshev-2') then
         do 80 i=1,n
            rwt(i)= sqrt (  1.d0 -pt(i)*pt(i) ) 
 80      continue   
         if(deriv) then
            do 90 i=1,n
               drwt(i) = - pt(i)/rwt(i)
               ddrwt(i) = ( - 1.d0 + pt(i)*drwt(i)/rwt(i) )/rwt(i)
 90         continue   
         endif
      elseif(wt.eq.'laguerre') then
         do 100 i=1,n
            rwt(i)= pt(i)**alpha*exp(-pt(i)) 
 100     continue   
         if(deriv) then
            do 110 i=1,n
               drwt(i) = ( - pt(i) + alpha/pt(i) ) * rwt(i)
               ddrwt(i) = (- pt(i) + alpha/pt(i) ) * drwt(i) -
     1                    (  1.d0 + alpha/(pt(i)*pt(i))) * rwt(i)
 110        continue
         endif   
      elseif(wt.eq.'jacobi') then
         do 200 i=1,n
            rwt(i)= (1.d0-pt(i))**alpha * (1.d0+pt(i))**beta 
 200     continue   
         if(deriv) then
            do 210 i=1,n
               fac1 = -alpha/(1.d0-pt(i))
               fac2 = beta/(1.d0+pt(i))
               drwt(i) = ( fac1 + fac2 )*rwt(i)
               ddrwt(i) = ( fac1 + fac2 )*drwt(i)
               fac1=fac1/(1.d0-pt(i))
               fac2=fac2/(1.d0+pt(i))
               ddrwt(i) = ddrwt(i) + ( fac1 + fac2 )*rwt(i)
 210        continue
         endif        
      elseif(wt.eq.'rys') then
         do 300 i=1,n
            rwt(i)=exp(-alpha*pt(i)*pt(i))
 300     continue   
         if(deriv) then
            do 310 i=1,n
               fac1 = -2.d0*pt(i)*alpha
               drwt(i) = fac1*rwt(i)
               ddrwt(i) = fac1*drwt(i) - 2.d0*alpha*rwt(i)
 310        continue
         endif
      else   
         call vfill(rwt,1.d0,n)
      endif
      return
      end       











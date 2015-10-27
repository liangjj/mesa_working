*deck lgrnge
      subroutine lgrnge(fx,x,fi,xi,n)
c***begin prologue     lgrnge
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            lagrange interpolate a function using
c***                   (n+1) points and values
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       lgrnge
c
      implicit integer (a-z)
      real*8 fx, x, ln, ld, fi, xi
      dimension fi(0:n), xi(0:n)
      common /io/ inp, iout
      fx=0.d0
      do 10 i=0,n
         ln=1.d0
         ld=1.d0
         do 20 j=0,i-1
            ln=ln*(x-xi(j))
            ld=ld*(xi(i)-xi(j))
 20      continue    
         do 30 j=i+1,n           
            ln=ln*(x-xi(j))
            ld=ld*(xi(i)-xi(j))
 30      continue
         fx=fx+fi(i)*ln/ld
 10   continue    
      return
      end
















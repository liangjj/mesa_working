c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Band3}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck band3.f
***begin prologue     band3
***date written       951229   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           
***author             schneider, barry (nsf)
***source             
***purpose            3 point second derivative
***                   
***references         
***routines called    
***end prologue       band3
      subroutine band3(band,v,d,del,addv,n,prn)
      implicit integer (a-z)
      real*8 band, v, d, del
      logical addv, prn
      character*80 title
      dimension band(2,n), v(n), d(2)
      common/io/inp, iout
      d(1) = -2.d0/(del*del)
      d(2) = 1.d0/(del*del)

c     fill the diagonal and off-diagonal working bands.

      do 10 i=1,n
         band(1,i) = -.5d0*d(1)
         band(2,i) = -.5d0*d(2)
 10   continue   

c     add in the potential if desired

      if(addv) then
         do 20 i=1,n
            band(1,i) = band(1,i) + v(i)
 20      continue    
      endif

c     print

      if(prn) then
         title='hamiltonian bands'
         call prntrm(title,band,2,n,2,n,iout)
      endif
      return
      end       





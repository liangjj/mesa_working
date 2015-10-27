c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Points}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck points.f
***begin prologue     points
***date written       951229   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           
***author             schneider, barry (nsf)
***source             
***purpose            points
***references         
***routines called    
***end prologue       points
      subroutine points(q,wt,a,del,n,diff,prn)
      implicit integer (a-z)
      real*8 q, wt, a, del, tmp
      character*80 title
      logical prn
      dimension q(n), wt(n), tmp(10)
      common/io/inp, iout
      q(1)=a
      do 10 i=2,n
         q(i) = q(i-1) + del
 10   continue   

c     Trapezoidal rule

      call int_2(tmp,del)
      do 20 i=1,n
         wt(i)=2.d0*tmp(1)
 20   continue
      if(prn) then
         title='points'
         call prntrm(title,q,n,1,n,1,iout)
      endif
      return
      end       


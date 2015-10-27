c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Cartesian Bloch Operator Matrix Elements}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck blxyz.f
c***begin prologue     blxyz
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           bloch matrix elements
c***author             schneider, barry (nsf)
c***source
c***purpose            generate bloch matrix elements
c***                   in cartesian coordinates.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       blxyz
      subroutine blxyz(blmat,p,dp,xl,el,xr,er,n)
      implicit integer (a-z)
      real*8 blmat, p, dp, xl, el, xr, er
      dimension blmat(n,n), p(n,n), dp(n,n)
      common/io/inp, iout 
      if(xl.eq.el) then
         do 10 i=1,n
            blmat(1,i) = blmat(1,i) + p(1,1)*dp(1,i)
 10      continue
      endif   
      if(xr.eq.er) then
         do 20 i=1,n
            blmat(n,i) = blmat(n,i) - p(n,n)*dp(n,i)
 20      continue   
      endif
      return
      end       

c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Cylindrical Bloch Operator Matrix Elements}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck blcyl.f
c***begin prologue     blcyl
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           bloch operator matrix elements
c***author             schneider, barry (nsf)
c***source
c***purpose            generate bloch operator matrix elements
c***                   in cylindrical coordinates.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       blcyl
      subroutine blcyl(blmat,p,dp,xl,el,xr,er,n,parity)
      implicit integer (a-z)
      real*8 blmat, p, dp, xl, el, xr, er
      real*8 two
      character*(*) parity
      dimension blmat(n,n), p(n,n), dp(n,n)
      common/io/inp, iout 
      data two / 2.d0 /
      if(parity.eq.'none') then
         if(xl.eq.el) then
            do 10 i=1,n
               blmat(1,i) = blmat(1,i) + xl * p(1,1) * dp(1,i)
 10         continue
         endif   
         if(xr.eq.er) then
            do 20 i=1,n
               blmat(n,i) = blmat(n,i) - xr * p(n,n) * dp(n,i)
 20         continue   
         endif
      elseif(parity.eq.'even'.or.parity.eq.'odd') then
         if(xl.eq.el) then
            do 30 i=1,n
               blmat(1,i) = blmat(1,i) + 
     1                     two * xl * xl * p(1,1) * dp(1,i)
 30         continue
         endif   
         if(xr.eq.er) then
            do 40 i=1,n
               blmat(n,i) = blmat(n,i) - 
     1                     two * xr * xr * p(n,n) * dp(n,i)
 40         continue   
         endif
      else
         call lnkerr('error in argument type')
      endif
      return
      end       

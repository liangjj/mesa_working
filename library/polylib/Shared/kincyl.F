c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Cylindrical Kinetic Energy Matrix Elements}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck kincyl.f
c***begin prologue     kincyl
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           kinetic energy
c***author             schneider, barry (nsf)
c***source
c***purpose            generate kinetic energy matrix elements
c***                   in cartesian coordinates.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       kincyl
      subroutine kincyl(kmat,p,dp,ddp,pt,wt,n,parity)
      implicit integer (a-z)
      real*8 kmat, p, dp, ddp, pt, wt
      real*8 two, four, six
      character*(*) parity
      dimension kmat(n,n), p(n,n), dp(n,n), ddp(n,n)
      dimension pt(n), wt(n)
      common/io/inp, iout 
      data one, two, four, six / 1.d0, 2.d0, 4.d0, 6.d0 /
      if(parity.eq.'none') then
         do 10 i=1,n
            do 20 j=1,n
               kmat(i,j) = kmat(i,j)  + 
     1                     p(i,i) * wt(i) * ( ddp(i,j) + 
     2                                        (one/pt(i)) * dp(i,j) )
 20         continue
 10      continue 
      elseif(parity.eq.'even'.or.parity.eq.'odd') then
         do 30 i=1,n
            do 40 j=1,n
               kmat(i,j) = kmat(i,j)  + p(i,i) * wt(i) * four *
     2                   ( pt(i) * pt(i) * ddp(i,j) + dp(i,j) )
 40         continue
 30      continue 
      else
         call lnkerr('error in argument type')
      endif
      return
      end       

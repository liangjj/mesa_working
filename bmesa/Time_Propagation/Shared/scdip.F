c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Scdip}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck scdip.f
***begin prologue     scdip
***date written       960723   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           copy
***author             schneider, barry (nsf)
***source             
***purpose            compute spatial factors for dipole potentials.
***references         
***routines called    
***end prologue       scdip
      subroutine scdip(v,tfac,q1,q2,q3,dim,n,nd,nc)
      implicit integer (a-z)
      real*8 v, tfac, q1, q2, q3, fac
      dimension nd(3)
      dimension v(n,nc,nc)
      common/io/inp, iout
      if(dim.eq.1) then
         do 10 i=1,n
            fac=tfac*q1(i)
            do 20 j=1,nc
               do 30 k=1,j
                  v(i,j,k) = v(i,j,k) + fac
 30            continue
 20         continue
 10      continue   
         do 40 j=1,nc
            do 50 k=1,j
               call copy(v(1,j,k),v(1,k,j),n)
 50         continue
 40      continue   
      elseif(dim.eq.2) then
         count=0
         do 60 i=1,nd(1)
            do 70 j=1,nd(2)
               fac=tfac*sqrt(q1(i)*q1(i)+q2(j)*q2(j))
               count=count+1
               do 80 k=1,nc
                  do 90 l=1,k      
                     v(count,k,l) = v(count,k,l) + fac
 90               continue
 80            continue
 70         continue
 60      continue
         do 100 k=1,nc
            do 110 l=1,k
               call copy(v(1,k,l),v(1,l,k),n)
 110        continue
 100     continue   
      elseif(dim.eq.3) then
         count=0
         do 200 i=1,nd(1)
            do 210 j=1,nd(2)
               do 220 k=1,nd(3)
                  fac=tfac*sqrt(q1(i)*q1(i)+q2(j)*q2(j)+q3(k)*q3(k))
                  count=count+1
                  do 230 l=1,nc
                     do 240 m=1,l      
                        v(count,l,m) = v(count,l,m) + fac
 240                 continue
 230              continue   
 220           continue
 210        continue
 200     continue
         do 300 l=1,nc
            do 310 m=1,l
               call copy(v(1,l,m),v(1,m,l),n)
 310        continue
 300     continue   
      else
         call lnkerr('dimension error')
      endif
      return
      end       



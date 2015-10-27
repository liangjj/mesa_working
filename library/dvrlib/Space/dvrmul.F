c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{DVR Matrix Vector Mutiply}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle 

*deck dvrmul.f 
c***begin prologue     dvrmul
c***date written       020206   (yymmdd)
c***revision date               (yymmdd)
c***keywords           matrix vector multiply, dvr
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            multiply a packed DVR hamiltonian on a vector
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       dvrmul
      subroutine dvrmul(vecin,vecout,diag,hbuf,ihbuf,n,nvc,nonz)
c
      implicit integer (a-z)
      real*8 vecin, vecout, diag, hbuf, hij
      dimension vecin(n,nvc), vecout(n,nvc)
      dimension diag(n), hbuf(*), ihbuf(2,*)
      common/io/inp, iout  
      call rzero(vecout,n*nvc)

c     do the diagonal

      do 10 i=1,nvc
         do 20 j=1,n
            vecout(j,i) = vecout(j,i) + diag(j)*vecin(j,i)
 20      continue   
 10   continue

c     Now the off-diagonal in packed form.

      do 30 i=1,nonz
         ii=ihbuf(1,i)
         jj=ihbuf(2,i)
         hij=hbuf(i)
         do 40 j=1,nvc
            vecout(ii,j) = vecout(ii,j) + hij*vecin(jj,j)
            vecout(jj,j) = vecout(jj,j) + hij*vecin(ii,j)
 40      continue   
 30   continue   
      return
      end

















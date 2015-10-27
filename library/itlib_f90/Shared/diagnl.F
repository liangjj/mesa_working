c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{DIAGNL: Calculate DVR Diagonals}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck diagnl.f
c***begin prologue     diagnl
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonal
c***author             schneider, barry (nsf)
c***source
c***description        calculate and store the full diagonal
c***                   for the 1,2 or 3 dimensional matrix.
c***                   then zero the diagonal elements of the
c***                   one-dimensional matrices to avoid overcounting.
c***references         
c
c***routines called    
c***end prologue       diagnl
      subroutine diagnl(diag,hx,hy,hz,v,n,nx,ny,nz,dim)
      implicit integer (a-z)
      real*8 diag, hx, hy, hz, v
      dimension diag(n), hx(nx,nx), hy(ny,ny), hz(nz,nz), v(n)
      common/io/inp, iout
      call copy(v,diag,n)
      if(dim.eq.1) then
         do 10 i=1,nx
            diag(i) = diag(i) + hx(i,i)
            hx(i,i) = 0.d0
 10      continue
      elseif(dim.eq.2) then
         count=0
         do 20 i=1,nx
            do 30 j=1,ny
               count = count + 1   
               diag(count) = diag(count ) + hx(i,i)
     1                                    + hy(j,j)
 30         continue
 20      continue
         do 40 i=1,nx
            hx(i,i) = 0.d0
 40      continue
         do 50 i=1,ny
            hy(i,i) = 0.d0
 50      continue   
      elseif(dim.eq.3) then
         count=0
         do 60 i=1,nx
            do 70 j=1,ny
               do 80 k=1,nz
                  count = count +1
                  diag(count) = diag(count) + hx(i,i) 
     1                                      + hy(j,j) 
     2                                      + hz(k,k)   
 80            continue
 70         continue
 60      continue
         do 100 i=1,nx
            hx(i,i) = 0.d0
 100     continue
         do 110 i=1,ny
            hy(i,i) = 0.d0
 110     continue   
         do 120 i=1,nz
            hz(i,i) = 0.d0
 120     continue   
      else
         call lnkerr('error in dimension')
      endif   
      return
      end


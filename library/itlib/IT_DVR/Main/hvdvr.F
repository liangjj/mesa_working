c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{HVDVR: Multiplication of Hamiltonian in DVR Representation on 
c         a Vector}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck hvdvr.f
c***begin prologue     hvdvr
c***date written       010828   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           matrix vector
c***author             schneider, barry (nsf)
c***source             tproplib              
c***purpose            matrix-vector, multiply, dvr. 
c***description        multiply a hamiltonian matrix, which is the sum
c***                   of up to 3 one-body matrices, onto a vector.            
c**
c***                          definitions
c                      m = n3*n2*n1
c                      hi(ni,ni) = matrix representation of T(i) + V0(i)
c                                  with zero diagonals
c                      diag      = full diagonal
c                      vecin(m,nvc)  = vecin(n3,n2,n1,nvc)
c                      vecout(m,nvc) = vecout(n3,n2,n1,nvc)
c
c                      it is assumed that;
c                      H = h1 + h2 + h3 + diag
c***references         
c***routines called    
c***end prologue       hvdvr
      subroutine hvdvr(h1,h2,h3,diag,vecin,vecout,n,n1,n2,n3,
     1                 nvc,dim,prnt)
      implicit integer (a-z)
      real*8 h1, h2, h3, diag, vecin, vecout
      logical prnt
      character*80 title
      dimension h1(n1,n1), h2(n2,n2), h3(n3,n3)
      dimension diag(n), vecin(n,nvc), vecout(n,nvc)
      common/io/inp, iout
c
      call rzero(vecout,n*nvc)
c
c
      if(dim.eq.1) then
         call h1v(h1,vecout,vecin,n1,nvc)
      elseif(dim.eq.2) then
         call h2v(h1,h2,vecout,vecin,n1,n2,nvc)
      elseif(dim.eq.3) then
         call h3v(h1,h2,h3,vecout,vecin,n1,n2,n3,nvc)
      endif 
c
c     take care of the diagonal
c
      do 10 i=1,nvc
         do 20 j=1,n
            vecout (j,i) = vecout(j,i) + diag(j)*vecin(j,i)
 20      continue
 10   continue
      if(prnt) then
         title='hamiltonian on vectors'
         call prntrm(title,vecout,n,nvc,n,nvc,iout)
      endif	    
      return
      end       




c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Band5}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck band7.f
***begin prologue     band7
***date written       951229   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           
***author             schneider, barry (nsf)
***source             
***purpose            7 point second derivative
***                   
***references         
***routines called    
***end prologue       band7
      subroutine band7(band,v,d,del,addv,nrow,n,prn)
      implicit integer (a-z)
      real*8 band, v, d, del
      logical addv, prn
      character*80 title
      dimension band(nrow,n), v(n), d(nrow)
      common/io/inp, iout
      d(4) = 2.d0/(720.d0*del*del)
      d(3) = -54.d0*d(4)
      d(2) = 540.d0*d(4)
      d(1) = -980.d0*d(4)
      d(4) = 4.d0*d(4)

c The structure of the input for general band matrices is different
c than for the tridiagonal case.  They are stored as column vectors
c beginning at the diagonal element M(i,i) and proceeding down the column
c in the usual storage pattern.  Since there are only nrow non-zero 
c elements for any column, that is all the storage required for the
c matrix.  If eigenvectors are wanted, provide storage in eigvec for them.
c Do all but the last three columns

      do 10 j=1,n-nrow+1
         do 20 row=1,nrow
            band(row,j) = -.5d0*d(row)
 20      continue
 10   continue

c     finish off n-2, n-1 and n

      band(1,n-2) = -.5d0*d(1)
      band(2,n-2) = -.5d0*d(2)
      band(3,n-2) = -.5d0*d(3)  
      band(1,n-1) = -.5d0*d(1)
      band(2,n-1) = -.5d0*d(2)
      band(1,n)   = -.5d0*d(1)  

c     add in the potential if desired

      if(addv) then
         do 30 j=1,n
            band(1,j) = band(1,j) + v(j)
 30      continue   
      endif
      if(prn) then
         title='hamiltonian columns'
         call prntrm(title,band,nrow,n,nrow,n,iout)
      endif
      return
      end       


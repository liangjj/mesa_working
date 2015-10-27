*deck xmtrx.f
c***begin prologue     xmtrx
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       xmtrx
      subroutine xmtrx(p,pn,x,wt,mat,n,npt)
      implicit integer (a-z)
      real*8 p, pn, x, wt, mat
      character*80 title 
      dimension p(npt,0:n-1), pn(npt,0:n-1), x(npt), wt(npt)
      dimension mat(n,n)
      common/io/inp, iout 
      call rzero(mat,n*n)
      do 10 i=1,n
         do 20 j=1,i
            do 30 k=1,npt
               mat(i,j) = mat(i,j) + p(k,i-1) * p(k,j-1) * wt(k) *x(k)
   30       continue
            mat(j,i) = mat(i,j)
   20    continue
   10 continue
      title='x-matrix of polynomials'
      call prntrm(title,mat,n,n,n,n,iout)
      call rzero(mat,n*n)
      do 40 i=1,n
         do 50 j=1,i
            do 60 k=1,npt
               mat(i,j) = mat(i,j) + pn(k,i-1) * pn(k,j-1) * wt(k) *x(k)
   60       continue
            mat(j,i) = mat(i,j)
   50    continue
   40 continue
      title='x-matrix of dvr polynomials'
      call prntrm(title,mat,n,n,n,n,iout)         
      return
      end       

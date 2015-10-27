*deck fdf.f
c***begin prologue     fdf
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto matrix elements of function
c***                   and first derivative.
c***                   
c***references         
c
c***routines called    
c***end prologue       fdf
      subroutine fdf(fa,dfb,wt,mat,n,m)
      implicit integer (a-z)
      real*8 fa, dfb, wt, mat
      dimension fa(n,m), dfb(n,m), wt(m), mat(m,m)
      common/io/inp, iout
      call rzero(mat,n*n)
      do 10 i=1,m
         do 20 j=1,m
            mat(i,j) = fa(i,i)*wt(i)*dfb(i,j)
 20      continue   
 10   continue   
      return
      end       




*deck fafb.f
c***begin prologue     fafb
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto overlap matrix elements.
c***                   
c***references         
c
c***routines called    
c***end prologue       fafb
      subroutine fafb(fa,fb,wt,mat,n)
      implicit integer (a-z)
      real*8 fa, fb, wt, mat
      dimension fa(n,n), fb(n,n), wt(n), mat(n,n)
      common/io/inp, iout
      call rzero(mat,n*n)
      do 10 i=1,n
         mat(i,i) = fa(i,i)*wt(i)*fb(i,i)
 10   continue   
      return
      end       

*deck faddfb.f
c***begin prologue     faddfb
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto matrix elements of function with 
c***                   second derivative of function.
c***                   
c***references         
c
c***routines called    
c***end prologue       faddfb
      subroutine faddfb(fa,dfb,ddfb,pt,wt,mat,typarg,n)
      implicit integer (a-z)
      real*8 fa, dfb, ddfb, pt, wt, mat
      character*(*) typarg
      dimension fa(n,n), dfb(n,n), ddfb(n,n), pt(n), wt(n), mat(n,n)
      common/io/inp, iout
      call rzero(mat,n*n)
      if(typarg.eq.'linear') then
         do 10 i=1,n
            do 20 j=1,n
               mat(i,j) = fa(i,i)*wt(i)*ddfb(i,j)
 20         continue   
 10      continue   
      elseif(typarg.eq.'quadratic') then
         do 30 i=1,n
            do 40 j=1,n
               mat(i,j) = fa(i,i)*wt(i)* ( 2.d0*dfb(i,j) 
     1                                                   + 
     2                                     4.d0*pt(i)*pt(i)*ddfb(i,j) )
 40         continue
 30      continue
      else
         call lnkerr('error in argument type')
      endif   
      return
      end       




*deck maket.f
c***begin prologue     maket
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            symmetrized kinetic energy matrix elements.
c***                   
c***description        computes the matrix of second derivatives
c                      plus block operator.      
c***references         
c
c***routines called    
c***end prologue       maket
      subroutine maket(fa,dfb,ddfb,pt,wt,secder,bloch,mat,typarg,n)
      implicit integer (a-z)
      real*8 fa, dfb, ddfb, pt, wt, secder, bloch, mat
      character*(*) typarg
      dimension fa(n,n), dfb(n,n), ddfb(n,n), pt(n), wt(n)
      dimension secder(n,n), bloch(n,n), mat(n,n)
      common/io/inp, iout
      if(typarg.eq.'linear') then
         call copy(secder,mat,n*n)
         call vadd(mat,mat,bloch,n*n)
      elseif(typarg.eq.'quadratic') then
         do 10 i=1,n
            do 20 j=1,n
               mat(i,j) = fa(i,i) * wt(i) * 4.d0 * 
     1                             ( pt(i)*pt(i)*ddfb(i,j) + dfb(i,j) )
 20         continue
 10      continue
         call vadd(mat,mat,bloch,n*n)
      else
         call lnkerr('error in argument type')
      endif
      return
      end       




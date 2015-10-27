*deck dfadfb.f
c***begin prologue     dfadfb
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto matrix elements of first derivative
c***                   of function with first derivative of function.
c***                   
c***references         
c
c***routines called    
c***end prologue       dfadfb
      subroutine dfadfb(dfa,dfb,pt,wt,mat,typarg,n)
      implicit integer (a-z)
      real*8 dfa, dfb, pt, wt, mat
      character*(*) typarg
      dimension dfa(n,n), dfb(n,n), pt(n), wt(n), mat(n,n)
      common/io/inp, iout
      call rzero(mat,n*n)
      if(typarg.eq.'linear') then
         do 10 i=1,n
            do 20 j=1,n
	       do 30 k=1,n
                  mat(i,j) = mat(i,j) + dfa(k,i)*wt(k)*dfb(k,j)
 30            continue	       
 20         continue   
 10      continue
      elseif(typarg.eq.'quadratic') then
         do 40 i=1,n
            do 50 j=1,n
               do 60 k=1,n
                  mat(i,j) = mat(i,j) + 4.d0*pt(k)*pt(k)*wt(k)*
     1                                       dfa(k,i)*dfb(k,j)
 60            continue
 50         continue
 40      continue
       else   
         call lnkerr('error in argument type')
      endif
      return
      end       










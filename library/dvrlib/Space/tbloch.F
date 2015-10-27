*deck tbloch.f
c***begin prologue     tbloch
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto matrix elements of bloch 
c***                   operator.
c***                   
c***description        computes the matrix
c                            
c                      M = - fa(xright) * fb'(xright) 
c                          + fa(xleft)  * fb'(xleft) 
c***references         
c
c***routines called    
c***end prologue       tbloch
      subroutine tbloch(fa,dfb,pt,mat,typarg,n)
      implicit integer (a-z)
      real*8 fa, dfb, pt, mat
      character*(*) typarg
      dimension fa(n,n), dfb(n,n), pt(n), mat(n,n)
      common/io/inp, iout
      call rzero(mat,n*n)
      if(typarg.eq.'linear') then
         do 10 i=1,n
            mat(n,i) = mat(n,i) - fa(n,n)*dfb(n,i)
            mat(1,i) = mat(1,i) + fa(1,1)*dfb(1,i)
 10      continue   
      elseif(typarg.eq.'quadratic') then
         do 20 i=1,n
            mat(n,i) = mat(n,i) - 2.d0*pt(n)*pt(n)*fa(n,n)*dfb(n,i)
            mat(i,i) = mat(1,i) + 2.d0*pt(1)*pt(1)*fa(1,1)*dfb(1,i)
 20      continue
      else
         call lnkerr('error in argument type')
      endif   
      return
      end       




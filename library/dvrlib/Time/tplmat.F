*deck tplmat.f
c***begin prologue     tplmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto matrix elements.
c***                   
c***references         
c
c***routines called    
c***end prologue       tplmat
      subroutine tplmat(p,dp,wt,mat,n)
      implicit integer (a-z)
      real*8 wt, p, dp, mat
      dimension wt(n), p(n,n), dp(n,n), mat(n,n)
      character*80 title
      common/io/inp, iout
      call fdf(p,dp,wt,mat,n,n) 
      return
      end       







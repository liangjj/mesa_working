*deck modply.f
c***begin prologue     modply
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
c***end prologue       modply
      subroutine modply(p,pn,ni,nj,mi,mj)
      implicit integer (a-z)
      real*8 p, pn
      dimension p(nj,*), pn(nj,*)
      common/io/inp, iout
      do 10 i=1,mi
         call copy(p(1,i),pn(1,i),mj)
 10   continue   
      return
      end       

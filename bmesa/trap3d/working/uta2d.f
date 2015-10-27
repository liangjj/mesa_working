*deck uta2d.f
c***begin prologue     uta2d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           transformation of vectors for two dimensional problem 
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       uta2d
      subroutine uta2d(resid,vec,u01,u02,t,n1,n2,m)
      implicit integer (a-z)
      real*8 resid, vec, u01, u02, t
      character*80 title
      dimension u01(n1,n1), u02(n2,n2)
      dimension vec(n2,n1,m), resid(n2,n1,m), t(*)
      common/io/inp, iout
c
c         transform residual from dvr to h0 representation
c 
      do 10 i=1,m
c
c        perform transformation on n2 variable
c      
         call ebtc(t,u02,resid(1,1,i),n2,n2,n1)
c
c        perform transformation on n1 variable
c        
         call ebc(vec(1,1,i),t,u01,n2,n1,n1)
 10   continue   
      return
      end       




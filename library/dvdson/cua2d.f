*deck cua2d.f
c***begin prologue     cua2d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           transformation of vectors for two dimensional problem 
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       cua2d
      subroutine cua2d(vec,resid,u01,u02,t,n1,n2,m)
      implicit integer (a-z)
      complex*16 resid, vec, u01, u02, t
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
         call cebc(t,u02,vec(1,1,i),n2,n2,n1)
c
c        perform transformation on n1 variable
c        
         call cebct(resid(1,1,i),t,u01,n2,n1,n1)
 10   continue   
      return
      end       




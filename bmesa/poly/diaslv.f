*deck diaslv.f
c***begin prologue     diaslv
c***date written       960430   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           linear equations, diagonal approximation
c***author             schneider, barry (nsf)
c***source             
c***purpose            diagonal approximation to solution of linear
c***                   system.
c***                   
c***references         
c
c***routines called    
c***end prologue       diaslv
      subroutine diaslv(resid,root,d,n)
      implicit integer (a-z)
      real*8 resid, root, d, nrzero, test
      dimension resid(n), d(n)
      common/io/inp, iout 
      data nrzero / 1.0d-06 /
      do 10 i=1,n
         test=root-d(i)
         if (abs(test).ge.nrzero) then
             resid(i)=resid(i)/test
         else     
             resid(i)=1.0d+00
         endif                
 10   continue
      return
      end       


*deck fileig.f
c***begin prologue     fileig
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description          
c***references       
c
c***routines called
c***end prologue       fileig
      subroutine fileig(eigin,eig,point,n,nr)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 eigin, eig
      dimension eigin(n), eig(nr,n)
      do 10 i=1,n
         eig(point,i)=eigin(i)
 10   continue
      return
      end

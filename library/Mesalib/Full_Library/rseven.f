*deck rseven
      subroutine rseven(wt,stp)
c***begin prologue     rseven
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes seven point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rseven
c
      implicit integer (a-z)
      real*8 wt, stp, d
      dimension wt(7,6)
      data d / 60480.d0 /
      common /io/ inp, iout
      wt(1,1)=19087.d0/d
      wt(2,1)=65112.d0/d
      wt(3,1)=-46461.d0/d
      wt(4,1)=37504.d0/d
      wt(5,1)=-20211.d0/d
      wt(6,1)=6312.d0/d
      wt(7,1)=-863.d0/d
      wt(1,2)=-863.d0/d
      wt(2,2)=25128.d0/d
      wt(3,2)=46989.d0/d
      wt(4,2)=-16256.d0/d
      wt(5,2)=7299.d0/d
      wt(6,2)=-2088.d0/d
      wt(7,2)=271.d0/d
      wt(1,3)=271.d0/d
      wt(2,3)=-2760.d0/d
      wt(3,3)=30819.d0/d
      wt(4,3)=37504.d0/d
      wt(5,3)=-6771.d0/d
      wt(6,3)=1608.d0/d
      wt(7,3)=-191.d0/d
      wt(1,4)=-191.d0/d
      wt(2,4)=1608.d0/d
      wt(3,4)=-6771.d0/d
      wt(4,4)=37504.d0/d
      wt(5,4)=30819.d0/d
      wt(6,4)=-2760.d0/d
      wt(7,4)=271.d0/d
      wt(1,5)=271.d0/d
      wt(2,5)=-2088.d0/d
      wt(3,5)=7299.d0/d
      wt(4,5)=-16256.d0/d
      wt(5,5)=46989.d0/d
      wt(6,5)=25128.d0/d
      wt(7,5)=-863.d0/d
      wt(1,6)=-863.d0/d
      wt(2,6)=6312.d0/d
      wt(3,6)=-20211.d0/d
      wt(4,6)=37504.d0/d
      wt(5,6)=-46461.d0/d
      wt(6,6)=65112.d0/d
      wt(7,6)=19087.d0/d
      do 10 i=1,7
         do 20 j=1,6
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
      return
      end
















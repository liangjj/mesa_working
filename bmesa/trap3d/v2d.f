*deck v2d.f
c***begin prologue     v2d
c***date written       970609   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential, two-dimension
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       v2d
      subroutine v2d(v,parts,n,nd)
      implicit integer (a-z)
      real*8 v, parts
      dimension v(*), parts(n,2), nd(2)
      common/io/inp, iout
      count=0
      do 10 i=1,nd(1)
         do 20 j=1,nd(2)
            count=count+1
            v(count) = v(count) + parts(i,1) + parts(j,2)
  20     continue
  10  continue
      return
      end       

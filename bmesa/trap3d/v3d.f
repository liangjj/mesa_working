*deck v3d.f
c***begin prologue     v3d
c***date written       970609   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential, three-dimension
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       v3d
      subroutine v3d(v,parts,n,nd)
      implicit integer (a-z)
      real*8 v, parts
      dimension v(*), parts(n,3), nd(3)
      common/io/inp, iout
      count=0
      do 10 i=1,nd(1)
         do 20 j=1,nd(2)
            do 30 k=1,nd(3)
               count=count+1
               v(count) = v(count) + parts(i,1) 
     1                             + parts(j,2)
     2                             + parts(k,3)
  30        continue   
  20     continue
  10  continue
      return
      end       

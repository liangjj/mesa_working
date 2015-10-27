*deck fnlham.f 
c***begin prologue     fnlham
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             fnlham
c***purpose            finals hamiltonian
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       fnlham
      subroutine fnlham(h_0,v,n,first,last)
c
      implicit integer (a-z)
      real*8 h_0, v
      common/io/inp, iout
      dimension h_0(n,n), v(*)
c
      count=0
      do 10 i=first,last
         count=count+1
         h_0(count,count) = h_0(count,count) + v(i)
 10   continue   
      return
      end



















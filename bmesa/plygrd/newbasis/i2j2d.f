*deck i2j2d.f
c***begin prologue     i2j2d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            interpolate a 2d vector from grid i
c***                   to a grid j.  the storage in the vector
c***                   corresponds to a loop structure looking like
c***                         do i=1,n1
c***                            do j=1,n2
c***                               do k=1,n3 
c***                                  ........
c***                   
c***references         
c
c***routines called    
c***end prologue       i2j2d
      subroutine i2j2d(p1ji,p2ji,vi,vj,scr,n1j,n2j,n1i,n2i,nvc)
      implicit integer (a-z)
      real*8 p1ji, p2ji, vi, vj
      dimension p1ji(n1j,n1i), p2ji(n2j,n2i)
      dimension vi(n2i,n1i,nvc), vj(n2j,n1j,nvc), scr(n2j,n1i,nvc)
      common/io/inp, iout
      call ebc(scr,p2ji,vi,n2j,n2i,n1i*nvc)
      do 10 i=1,nvc
         call ebct(vj(1,1,i),scr(1,1,i),p1ji,n2j,n1i,n1j)
 10   continue
      return
      end       


*deck i2j3d.f
c***begin prologue     i2j3d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            interpolate a 3d vector from grid i
c***                   to grid j.  the storage in the vector
c***                   corresponds to a loop structure looking like
c***                         do i=1,n1
c***                            do j=1,n2
c***                               do k=1,n3 
c***                                  ........
c***                   the scratch array sc1 and vj can be
c***                   implicitly equivalenced in the call since
c***                   sc1 is only needed as an intermediate.       
c***references         
c
c***routines called    
c***end prologue       i2j3d
      subroutine i2j3d(p1ji,p2ji,p3ji,vi,vj,sc1,sc2,n1j,n2j,n3j,
     1                 n1i,n2i,n3i,nvc)
      implicit integer (a-z)
      real*8 p1ji, p2ji, p3ji
      real*8 sc1, sc2, vi, vj
      dimension p1ji(n1j,n1i), p2ji(n2j,n2i), p3ji(n3j,n3i)
      dimension vi(n3i,n2i,n1i,nvc), vj(n3j,n2j,n1j,nvc)
      dimension sc1(n3j,n2i,n1i,nvc), sc2(n3j,n2j,n1i,nvc)
      common/io/inp, iout
      call ebc(sc1,p3ji,vi,n3j,n3i,n2i*n1i*nvc)
      do 10 i=1,nvc
         do 20 j=1,n1i
            call ebct(sc2(1,1,j,i),sc1(1,1,j,i),p2ji,n3j,n2i,n2i)
 20      continue   
 10   continue
      do 30 i=1,nvc
         call ebct(vj(1,1,1,i),sc2(1,1,1,i),p1ji,n3j*n2j,n1i,n1j)
 30   continue
      return
      end       


*deck i2j4d.f
c***begin prologue     i2j4d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            interpolate a 4d vector from a grid i
c***                   to grid j.  the storage in the vector
c***                   corresponds to a loop structure looking like
c***                         do i=1,n1
c***                            do j=1,n2
c***                               do k=1,n3 
c***                                  do l=1,n4
c***                                  ........
c***                   the scratch array sc1 and sc3 can be
c***                   implicitly equivalenced in the call 
c***                   as well as sc2 and vj.       
c***references         
c
c***routines called    
c***end prologue       i2j4d
      subroutine i2j4d(p1ji,p2ji,p3ji,p4ji,vi,vj,sc1,sc2,sc3,
     1                 n1j,n2j,n3j,n4j,n1i,n2i,n3i,n4i,nvc)
      implicit integer (a-z)
      real*8 p1ji, p2ji, p3ji, p4ji
      real*8 sc1, sc2, vi, vj
      dimension p1ji(n1j,n1i), p2ji(n2j,n2i)
      dimension p3ji(n3j,n3i), p4ji(n4j,n4i)
      dimension vi(n4i,n3i,n2i,n1i,nvc), vj(n4j,n3j,n2j,n1j,nvc)
      dimension sc1(n4j,n3i,n2i,n1i,nvc), sc2(n4j,n3j,n2i,n1i,nvc)
      dimension sc3(n4j,n3j,n2j,n1i,nvc)
      common/io/inp, iout
      call ebc(sc1,p4ji,vi,n4j,n4i,n3i*n2i*n1i*nvc)
      do 10 i=1,nvc
         do 20 j=1,n1i
            do 30 k=1,n2i
               call ebct(sc2(1,1,k,j,i),sc1(1,1,k,j,i),p3ji,
     1                   n4j,n3i,n3j)
 30         continue   
 20      continue   
 10   continue
      do 40 i=1,nvc
         do 50 j=1,n1i
            call ebct(sc3(1,1,1,j,i),sc2(1,1,1,j,i),p2ji,
     1                n4j*n3j,n2i,n2j)
 50      continue
 40   continue   
      do 60 i=1,nvc
         call ebct(vj(1,1,1,1,i),sc3(1,1,1,1,i),p1ji,
     1             n4j*n3j*n2j,n1i,n1j)
 60   continue
      return
      end       


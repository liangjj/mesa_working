*deck tstmem
      subroutine tstmem(val,n)
      implicit integer (a-z)
      real*8 val, z, start       
      integer*8 pnt(100)
      common/io/inp, iout
      pointer(p,z(1)), (p,iz(1))
      do while(n.gt.0) 
         call filtst(pnt(n),val)
         p=pnt(n)
         write(iout,*) pnt(n)
         val=val*.5d0
         n=n-1
c         call tstmem(val,n)
         write(iout,1) (z(i),i=1,10)
c         call memory(-ngot,p,idum,'tstmem',idum)
      enddo
       return
 1    format(/,5x,5e15.8)
       end

*deck filtst
      subroutine filtst(p,val)
      implicit integer (a-z)
      real*8 val, z, start       
      common/io/inp, iout
      pointer(p,z(1))
      need=wptoin(10)
      call memory(need,p,ngot,'tstmem',0)
      write(iout,*) p
      start=val
      do 10 i=1,10
         z(i)=start
         start=start+1.d0
 10   continue
      write(iout,1) (z(i),i=1,10)
      return
 1    format(/,5x,5e15.8)
      end
      

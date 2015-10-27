      program testmem
      implicit integer(a-z)
      real*8 z
      pointer (p,z(1)), (p,a(1))
      common /io/ inp, iout
      call drum
      need=1000
      call memory(need,p,ngot,'test',0)
      do 100 i=1,100
         z(i)=i
 100  continue   
      write(iout,*) (z(i),i=1,100)
      call memory(-ngot,p,idum,'test',idum)
      call memory(need,p,ngot,'test',0)
      call chainx(0)
      stop 
      end

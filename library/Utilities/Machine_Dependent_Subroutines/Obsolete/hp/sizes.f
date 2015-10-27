*deck %W% %G%
      function sizes(mx)
      implicit integer(a-z)
c
c     this function map the maximum memory on a machine on the network
c     to its name.
c
      integer sizes
      character*(*) mx
c
c
      parameter (nmx=4)
      character*32 mch(nmx)
      integer mchsiz(4)
      data mchsiz/64,32,64,0/
      data mch/'lithium','sodium','potassium','unknown'/
      save mchsiz,mch
c
c     return the maximum memory available on various machines in megabytes.
      do 10 i=1,nmx
         if(mx.eq.mch(i)) then
            sizes=mchsiz(i)
         endif
   10 continue
c
c
      return
      end

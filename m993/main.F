*deck main.f 
      program main
c
      implicit integer (a-z)
      character*4096 ops
      character*80 chrkey, coord
      character*1 itoc
#ifdef DECPOINTER
      integer*8 phamil
#else
#ifdef SGIPOINTER
      integer*4 phamil
#endif SGIPOINTER
#endif DECPOINTER
      dimension pham(4)
      dimension coord(2)
      common/io/inp, iout      
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      do 10 i=1,2
         coord=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         call basis(pham(i),coord(i),nphy)
         write(iout,1) nphy
 10   continue   
 1    format(/,1x,'physical size of Hamiltonian = ',i4)
      call chainx(0)
      stop
      end

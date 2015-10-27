      subroutine geneng (energy,enbnd,en,rk,ncst,nsts,ncmax)
      implicit integer(a-z)
      common /io/ inp, iout
      real *8energy, enbnd, en, rk, ediff, engdif
      dimension en(ncmax,nsts), rk(ncmax,nsts), ncst(nsts)
      dimension enbnd(nsts)
      write (iout,50) energy
      write (iout,30)
      do 20 j=1,nsts
      ediff=enbnd(1)-enbnd(j)
      ediff=2.d+00*ediff
      engdif=energy-ediff
      write (iout,40) j,engdif
      nc=ncst(j)
      do 10 k=1,nc
      en(k,j)=engdif
      rk(k,j)=sqrt(abs(en(k,j)))
   10 continue
   20 continue
      return
c
   30 format (/,5x,'channel energies')
   40 format (/,5x,'state',2x,i3,5x,'energy',2x,f12.6)
   50 format (//,20x,'scattering calculation at energy',1x,f12.8)
      end

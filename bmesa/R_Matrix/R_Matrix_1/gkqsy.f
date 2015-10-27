*deck gkqsy.f
c***begin prologue     gkqsy
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            R-matrix sums for symmetric wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       gkqsy
      subroutine gkqsy(pham,pgam,ngot,n,dim,prn)
      implicit integer (a-z)
      integer*8 pham, pgam, ph, phamx0, pscr
      real*8 ham, gam, hamx0, scr
      character*80 title
      logical prn
      dimension pham(4,2), n(4,2)
      common/io/inp, iout
      pointer (ph,ham(1))
      pointer (pgam,gam(1))
      pointer (phamx0,hamx0(1))
      pointer (pscr,scr(1))
      phamx0=pham(1,2)
      ph=pham(dim+1,1)
      hmtx0=1
      vx0=hmtx0+n(1,2)*n(1,2)
      eigvcx0=vx0+n(1,2)
      gamx0=eigvcx0+n(1,2)*n(1,2)      
      gamx0=gamx0+n(1,2)
      eigvlx0=gamx0+n(1,2)
      srfx0=eigvlx0+n(1,2)
      srfx0=srfx0+1
      eigvc=1
      sum=1
      need=wpadti(sum+n(1,2)*n(dim+1,1)) 
      call getmem(need,pscr,junk,'scr',0) 
      gamma=1
      need=wpadti(gamma+n(1,2)*n(dim+1,1))
      call getmem(need,pgam,ngot,'gamma',0) 
      call gssums(ham(eigvc),hamx0(srfx0),scr(sum),n(dim+1,1),
     1                                             n(1,2),prn)
      call ebtc(gam(gamma),hamx0(eigvcx0),scr(sum),
     1          n(1,2),n(1,2),n(dim+1,1))
      if(prn) then
         title='r-matrix amplitudes'
         call prntrm(title,gam(gamma),n(1,2),n(dim+1,1),
     1                                n(1,2),n(dim+1,1),iout)
      endif
      call iosys('write integer "number of channels" to ham',
     1            1,n(1,2),0,' ')
      call iosys('write real "channel energies" to ham',
     1            n(1,2),hamx0(eigvlx0),0,' ')
      call iosys('write real "h amplitudes" to ham',gam(gamma),
     1            n(1,2)*n(dim+1,1),0,' ')
      call getmem(-junk,pscr,idum,'scr',idum) 
      return      
      end       



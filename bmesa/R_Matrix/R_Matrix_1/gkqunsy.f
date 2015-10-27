*deck gkqunsy.f
c***begin prologue     gkqunsy
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            R-matrix sums for non-symmetric wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       gkqunsy
      subroutine gkqunsy(pham,pgam,ngot,n,dim,prn)
      implicit integer (a-z)
      integer*8 pham, pgam, ph, phamx0, phamy0
      real*8 ham, gam, hamx0, hamy0, scr
      character*80 title
      logical prn
      dimension pham(4,2), n(4,2)
      common/io/inp, iout
      pointer (ph,ham(1))
      pointer (pgam,gam(1))
      pointer (phamx0,hamx0(1))
      pointer (phamy0,hamy0(1))
      pointer (pscr,scr(1))
c
c
      phamx0=pham(1,2)
      phamy0=pham(2,2)
      ph=pham(dim+1,1)
      hmtx0=1
      vx0=hmtx0+n(1,2)*n(1,2)
      eigvcx0=vx0+n(1,2)
      gmx0=eigvcx0+n(1,2)*n(1,2) 
      gmx0=gmx0+n(1,2)       
      eigvlx0=gmx0+n(1,2)
      srfx0=eigvlx0+n(1,2)
      srfx0=srfx0+1
      hmty0=1
      vy0=hmty0+n(2,2)*n(2,2)
      eigvcy0=vy0+n(2,2)
      gmy0=eigvcy0+n(2,2)*n(2,2)      
      gmy0=gmy0+n(2,2)
      eigvly0=gmy0+n(2,2)
      srfy0=eigvly0+n(2,2)
      srfy0=srfy0+1
      eigvc=1
c
c     the arrays sum and gam are filled as if there was one
c     matrix with total column length nc=n(1)+n(2)
c 
      nc=n(1,2)+n(2,2)
      call iosys('write integer "number of channels" to ham',
     1            1,nc,0,' ')
      sum=1
      sumx=1
      sumy=sumx+n(1,2)
      need=wpadti(sum+nc*n(dim+1,1))
      call getmem(need,pscr,junk,'scr',0) 
      gamma=1
      gamx=1
      gamy=gamx+n(1,2)
      need=wpadti(gamma+nc*n(dim+1,1))
      call getmem(need,pgam,ngot,'gamma',0) 
      call gusums(ham(eigvc),hamx0(srfx0),hamy0(srfy0),
     1            scr(sumx),scr(sumy),n(dim+1,1),n(1,2),n(2,2),nc,prn)
      call ebtcxx(gam(gamx),hamx0(eigvcx0),scr(sumx),
     1            n(1,2),n(1,2),n(dim+1,1),nc,n(1,2),nc)
      call ebtcxx(gam(gamy),hamy0(eigvcy0),scr(sumy),
     1            n(2,2),n(2,2),n(dim+1,1),nc,n(2,2),nc)   
      if(prn) then
         title='gammax r-matrix amplitudes'
         call prntrm(title,gam(gamx),n(1,2),n(dim+1,1),nc,
     1                                      n(dim+1,1),iout)
         title='gammay r-matrix amplitudes'
         call prntrm(title,gam(gamy),n(2,2),n(dim+1,1),nc,
     1                                      n(dim+1,1),iout)
      endif
      call iosys('write real "h amplitudes" to ham',
     1            nc*n(dim+1,1),gam(gamma),0,' ')
      call iosys('write integer "number of x channels" to ham',
     1            1,n(1,2),0,' ')
      call iosys('write real "x channel energies" to ham',
     1            n(1,2),hamx0(eigvlx0),0,' ')
      call iosys('write integer "number of y channels" to ham',
     1            1,n(2,2),0,' ')
      call iosys('write real "y channel energies" to ham',
     1            n(2,2),hamy0(eigvly0),0,' ')
      call getmem(-junk,pscr,idum,'scr',idum) 
      return      
      end       



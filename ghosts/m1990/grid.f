*deck @(#)grid.f	4.1  7/7/93
      subroutine grid(gridx,gridy,gridz,ngridx,ngridy,ngridz,
     $               maxpts,xstart,ystart,zstart,incrx,incry,incrz,
     $               c,nat,step)
c***begin prologue     grid
c***date written       891219   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           m1990, link 1990, grid
c***author             martin, richard (lanl)
c***source             @(#)grid.f	4.1   7/7/93
c***purpose            produces the grid for evaluation of the wavefunction.
c***description
c     grid produces a cubic grid of equally spaced points which should
c       encompass most of the amplitude of the wavefunction. the positions
c       of the atomic centers arrive in c, the grid points exit in cgrid.
c
c***references
c
c***routines called
c
c***end prologue       grid
      implicit integer(a-z)
      real*8 gridx(maxpts),gridy(maxpts),gridz(maxpts),c(3,nat)
      real*8 xmax,ymax,zmax,xmin,ymin,zmin
      real*8 xdim,ydim,zdim
      real*8 zero,two,step
      real*8 incrx,incry,incrz,xstart,ystart,zstart,min,max
c
      parameter (zero=0.0d+00,two=2.0d+00)
c
      common/io/inp,iout
c
 1000 format(1x,'grid parameters:')
 1010 format(5x,'xstart:',f10.5,2x,'stepx :',f10.5,2x,'ngridx:',i5,
     $      /5x,'ystart:',f10.5,2x,'stepy :',f10.5,2x,'ngridy:',i5,
     $      /5x,'zstart:',f10.5,2x,'stepz :',f10.5,2x,'ngridz:',i5)
c
c     determine the rough dimensions of the molecule.
      xmin=zero
      xmax=zero
      ymin=zero
      ymax=zero
      zmin=zero
      zmax=zero
      do 10 at=1,nat
         xmin=min(xmin,c(1,at))
         ymin=min(ymin,c(2,at))
         zmin=min(zmin,c(3,at))
         xmax=max(xmax,c(1,at))
         ymax=max(ymax,c(2,at))
         zmax=max(zmax,c(3,at))
   10 continue
c
c     xmin,xmax denote the range in x of atomic centers.
c     using a fixed stepsize, determine the number of points neede
c     to span each dimension.
      xdim=xmax-xmin
      ydim=ymax-ymin
      zdim=zmax-zmin
c     add 6.0 bohr and assume that encompasses the orbital.
      xdim=xdim+6.0d+00
      ydim=ydim+6.0d+00
      zdim=zdim+6.0d+00
c
c     determine the number of points. for right now, use a fixed step size
c     in all directions.
      ngridx=int(xdim/step) +1
      ngridy=int(ydim/step) +1
      ngridz=int(zdim/step) +1
c
c     load the grid.
      if(ngridx.ge.maxpts.or.ngridy.ge.maxpts.or.ngridz.ge.maxpts)
     $   call lnkerr('too many grid points..')
      pt=0
      incrx=step
      incry=step
      incrz=step
      xstart=xmin-3.0d+00
      ystart=ymin-3.0d+00
      zstart=zmin-3.0d+00
      do 40 pt=1,ngridx
         gridx(pt)=xstart+incrx*(pt-1)
   40 continue
      do 45 pt=1,ngridy
         gridy(pt)=ystart+incry*(pt-1)
   45 continue
      do 50 pt=1,ngridz
         gridz(pt)=zstart+incrz*(pt-1)
   50 continue
c
c
      write(iout,1000)
      write(iout,1010) xstart,incrx,ngridx,ystart,incry,ngridy,
     $                 zstart,incrz,ngridz
c
c
      return
      end

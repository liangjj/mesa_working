*deck dtacol.f
c***begin prologue     dtacol
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           collocation, scattering
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            data for collocation calculation 
c***references         
c
c***routines called 
c***                   
      subroutine dtacol(card,type,rbeg,rend,nreg,npts,nwts,ptot,wtot,
     1                  ptcl,wtcl,ncorfn,ng,nrmax,ngmax)
      implicit integer (a-z)
      real*8 fpkey, rbeg, rend
      character*(*) card
      character*3 itoc
      character*16 type, chrkey
      dimension npts(nrmax,ngmax), nwts(nrmax,ngmax), rbeg(nrmax)
      dimension rend(nrmax), ptot(ngmax), wtot(ngmax)
      common /io/ inp, iout
      call izero(ptot,ngmax)
      call izero(wtot,ngmax)
      ptcl=0
      wtcl=0
      write(iout,1)
      type=chrkey(card,'type-quadrature','legendre',' ')
      write(iout,2) type
      nreg=intkey(card,'number-of-coarse-grid-integration-regions',
     1            1,' ')
      if (type.eq.'newton-cotes') then
         do 10 i=1,nreg
            npts(i,1)=intkey(card,'number-of-coarse-grid-points-region-'
     1                     //itoc(i),3,' ')
            nwts(i,1)=npts(i,1)*(npts(i,1)-1)
            rbeg(i)=fpkey(card,'starting-r-region-'//itoc(i),0.d0,' ')
            rend(i)=fpkey(card,'final-r-region-'//itoc(i),10.d0,' ')
            ptcl=ptcl+npts(i,1)
            wtcl=wtcl+nwts(i,1)
            ptot(1)=ptot(1)+npts(i,1)
            wtot(1)=wtot(1)+nwts(i,1)
            do 20 j=2,ng
               npts(i,j)=npts(i,j-1)+npts(i,j-1)-1
               nwts(i,j)=npts(i,j)*(npts(i,j)-1)
               ptcl=ptcl+npts(i,j)
               wtcl=wtcl+nwts(i,j)
               ptot(j)=ptot(j)+npts(i,j)
               wtot(j)=wtot(j)+nwts(i,j)
 20         continue   
            write(iout,3) i, npts(i,1), rbeg(i), rend(i)
 10      continue   
      else
         do 30 i=1,nreg
            npts(i,1)=intkey(card,'number-of-coarse-grid-points-region-'
     1                       //itoc(i),4,' ')
            nwts(i,1)=npts(i,1)
            ptcl=ptcl+npts(i,1)
            wtcl=wtcl+nwts(i,1)
            ptot(1)=ptot(1)+npts(i,1)
            wtot(1)=wtot(1)+nwts(i,1)
            rbeg(i)=fpkey(card,'starting-r-region-'//itoc(i),0.d0,' ')
            rend(i)=fpkey(card,'final-r-region-'//itoc(i),10.d0,' ')
            do 40 j=2,ng
               npts(i,j)=npts(i,j-1)+npts(i,j-1)-1
               nwts(i,j)=npts(i,j)
               ptcl=ptcl+npts(i,j)
               wtcl=wtcl+nwts(i,j)
               ptot(j)=ptot(j)+npts(i,j)
               wtot(j)=wtot(j)+nwts(i,j)
 40         continue   
            write(iout,3) i, npts(i,1), rbeg(i), rend(i)
 30      continue
c     add an additional region with one point and no weight since gauss
c     quadratures do not contain the beginning or endpoint of the interval
         nreg=nreg+1
         rbeg(nreg)=rend(nreg-1)
         rend(nreg)=rbeg(nreg)
         do 50 i=1,ng
            npts(nreg,i)=1
            nwts(nreg,i)=1
            ptcl=ptcl+1
            wtcl=wtcl+1
            ptot(i)=ptot(i)+1
            wtot(i)=wtot(i)+1
 50      continue    
         i=nreg
         write(iout,3) i, npts(i,1), rbeg(i), rend(i)
      endif
      ncorfn=intkey(card,'number-of-coarse-grid-'//
     1                   'collocation-functions',ptcnt, ' ')
      write(iout,4)
      write(iout,5) ptcl, wtcl
      write(iout,6) (ptot(i),i=1,ng)
      write(iout,7) (wtot(i),i=1,ng)
      return
 1    format(/,1x,'collocation region coarse grid data')
 2    format(1x,'type quadrature =',a16)
 3    format(/,5x,'integration region = ',i3,1x,'number points = ',i4,
     1       /,5x,'starting value     = ',e15.8,1x,'ending value = ',
     2             e15.8)
 4    format(/,1x,'grid summary:')
 5    format(/,'total points all grids = ',i5,/,
     1         'total weights all grids = ',i5)
 6    format(/,'total points each subgrid =',(/,29x,10(i4,1x)))
 7    format(/,'total weights each subgrid =',(/,29x,10(i4,1x)))
      end

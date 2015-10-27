*deck dtaft.f
c***begin prologue     dtaft
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           collocation, scattering
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            data for fitting region
c***references         
c
c***routines called 
c***                   
      subroutine dtaft(card,type,rbeg,rend,nreg,npts,nwts,ptot,wtot,
     1                 ptft,wtft,ng,nrmax,ngmax)
      implicit integer (a-z)
      real*8 fpkey, rbeg, rend, delreg
      character*(*) card
      character*3 itoc
      character*16 type, chrkey
      dimension npts(nrmax,ngmax), nwts(nrmax,ngmax)
      dimension rbeg(nrmax), rend(nrmax), ptot(ngmax), wtot(ngmax)
      common /io/ inp, iout
      call izero(ptot,ngmax)
      call izero(wtot,ngmax)
      ptft=0
      wtft=0
      nreg=intkey(card,'number-of-coarse-grid-integration-regions',
     1                  1,' ')
      write(iout,1)
      type=chrkey(card,'type-step-size','variable',' ')
      write(iout,2) type
      if (type.eq.'default') then
          i=1
          rbeg(i)=fpkey(card,'starting-r',0.d0,' ')
          rend(i)=fpkey(card,'final-r',10.d0,' ')          
          npts(i,1)=intkey(card,'number-of-coarse-grid-points-per-'//
     1                         'region',5,' ')
          nwts(i,1)=npts(i,1)*(npts(i,1)-1)
          ptft=ptft+npts(i,1)
          wtft=wtft+nwts(i,1)
          ptot(1)=ptot(1)+npts(i,1)
          wtot(1)=wtot(1)+nwts(i,1)          
          delreg=(rend(i)-rbeg(i))/nreg
          rend(i)=rbeg(i) + delreg
          write(iout,3) i, npts(i,1), rbeg(i), rend(i)
          do 10 i=2,nreg
             npts(i,1)=npts(i-1,1)
             nwts(i,1)=nwts(i-1,1)
             rbeg(i)=rend(i-1)
             rend(i)=rend(i-1) + delreg
             write(iout,3) i, npts(i,1), rbeg(i), rend(i)
             ptft=ptft+npts(i,1)
             wtft=wtft+nwts(i,1)
             do 20 j=2,ng  
                npts(i,j)=npts(i,j-1)+npts(i,j-1)-1
                nwts(i,j)=npts(i,j)*(npts(i,j)-1)
                ptft=ptft+npts(i,j)
                wtft=wtft+nwts(i,j)
                ptot(j)=ptot(j)+npts(i,j)
                wtot(j)=wtot(j)+nwts(i,j)
 20          continue   
 10       continue   
      else
          nreg=intkey(card,'number-of-coarse-grid-integration-regions',
     1                1,' ')
          do 30 i=1,nreg
             npts(i,1)=intkey(card,'number-of-coarse-grid-points-'//
     1                        'region-'//itoc(i),3,' ')
             nwts(i,1)=npts(i,1)
             rbeg(i)=fpkey(card,'starting-r-region-'//itoc(i),0.d0,' ')
             rend(i)=fpkey(card,'final-r-region-'//itoc(i),10.d0,' ')
             write(iout,3) i, npts(i,1), rbeg(i), rend(i)
             ptft=ptft+npts(i,1)
             wtft=wtft+nwts(i,1)
             ptot(1)=ptot(1)+npts(i,1)
             wtot(1)=wtot(1)+nwts(i,1)
             do 40 j=2,ng
                npts(i,j)=npts(i,j-1)+npts(i,j-1)-1
                nwts(i,j)=npts(i,j)
                ptft=ptft+npts(i,j)
                wtft=wtft+nwts(i,j)
                ptot(j)=ptot(j)+npts(i,j)
                wtot(j)=wtot(j)+nwts(i,j)
 40          continue   
 30       continue
c     add an additional region with one point and no weight since gauss
c     quadratures do not contain the beginning or endpoint of the interval
          nreg=nreg+1
          rbeg(nreg)=rend(nreg-1)
          rend(nreg)=rbeg(nreg)
          do 50 i=1,ng
             npts(nreg,i)=1
             nwts(nreg,i)=1
             ptft=ptft+1
             wtft=wtft+1
             ptot(i)=ptot(i)+1
             wtot(i)=wtot(i)+1
 50       continue    
          i=nreg
          write(iout,3) i, npts(i,1), rbeg(i), rend(i)
      endif
      write(iout,4)
      write(iout,5) ptft, wtft
      write(iout,6) (ptot(i),i=1,ng)
      write(iout,7) (wtot(i),i=1,ng)
      return
 1    format(/,1x,'fitting region coarse grid data')
 2    format(1x,'type step size =',a16)
 3    format(/,5x,'integration region = ',i3,1x,'number points = ',i4,
     1       /,5x,'starting value     = ',e15.8,1x,'ending value = ',
     2             e15.8)
 4    format(/,1x,'grid summary:')
 5    format(/,'total points all grids = ',i5,/,
     1         'total weights all grids = ',i5)
 6    format(/,'total points each subgrid =',(/,29x,10(i4,1x)))
 7    format(/,'total weights each subgrid =',(/,29x,10(i4,1x)))
      end

*deck dtafd.f
c***begin prologue     dtafd
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           collocation, scattering
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            data for finite difference region calculation
c***references         
c
c***routines called 
c***                   
      subroutine dtafd(card,rbeg,rend,nreg,npts,nwts,ptfd,wtfd,
     1                 points,dim)
      implicit integer (a-z)
      real*8 fpkey, rbeg, rend, dele, rb, re
      character*(*) card
      character*3 itoc
      dimension npts(dim,2), nwts(dim,2), rbeg(dim), rend(dim)
      dimension ptfd(2), wtfd(2), points(2)
      logical type, logkey
      common /io/ inp, iout
      write(iout,1)
      ptfd(1)=0
      wtfd(1)=0
      ptfd(2)=0
      wtfd(2)=0
      points(1)=0
      points(2)=0
      type=logkey(card,'standard-mesh',.false.,' ')
      nreg=intkey(card,'number-of-integration-regions',1,' ')
      if (type) then
          rb=fpkey(card,'starting-r',0.d0,' ')
          re=fpkey(card,'final-r',10.d0,' ')         
          dele=(re-rb)/nreg
          npts(1,1)=5
          nwts(1,1)=20
          npts(1,2)=9
          nwts(1,2)=72
          rbeg(1)=rb
          rend(1)=rbeg(1)+dele
          ptfd(1)=npts(1,1)
          wtfd(1)=nwts(1,1)
          ptfd(2)=npts(1,2)
          wtfd(2)=nwts(1,2)
          points(1)=npts(1,1)-1
          points(2)=npts(1,2)-1
          i=1
          write(iout,2) i, npts(i,1), rbeg(i), rend(i)
          do 10 i=2,nreg
             npts(i,1)=npts(i-1,1)
             nwts(i,1)=nwts(i-1,1)
             npts(i,2)=npts(i-1,2)
             nwts(i,2)=nwts(i-1,2)
             rbeg(i)=rend(i-1)
             rend(i)=rbeg(i)+dele
             write(iout,2) i, npts(i,1), rbeg(i), rend(i)
             ptfd(1)=ptfd(1)+npts(i,1)
             wtfd(1)=wtfd(1)+nwts(i,1)
             ptfd(2)=ptfd(2)+npts(i,2)
             wtfd(2)=wtfd(2)+nwts(i,2)
             points(1)=points(1)+npts(i,1)-1
             points(2)=points(2)+npts(i,2)-1
 10       continue    
      else
          do 20 i=1,nreg
             npts(i,1)=intkey(card,'number-of-grid-points-'//
     1                            'region-'//itoc(i),3,' ')
             npts(i,2)=2*npts(i,1)-1
             if(npts(i,2).gt.9) then
                npts(i,1)=5
                npts(i,2)=9
             endif
             nwts(i,1)=npts(i,1)*(npts(i,1)-1)
             nwts(i,2)=npts(i,2)*(npts(i,2)-1)
             rbeg(i)=fpkey(card,'starting-r-region-'//itoc(i),0.d0,' ')
             rend(i)=fpkey(card,'final-r-region-'//itoc(i),10.d0,' ')
             write(iout,2) i, npts(i,1), rbeg(i), rend(i)
             ptfd(1)=ptfd(1)+npts(i,1)
             wtfd(1)=wtfd(1)+nwts(i,1)
             ptfd(2)=ptfd(2)+npts(i,2)
             wtfd(2)=wtfd(2)+nwts(i,2)
             points(1)=points(1)+npts(i,1)-1
             points(2)=points(2)+npts(i,2)-1
 20       continue
      endif
      points(1)=points(1)+1   
      points(2)=points(2)+1   
      write(iout,3)
      write(iout,4) ptfd(1), wtfd(1)
      write(iout,5) points(1)
      return
 1    format(/,1x,'finite difference region grid data')
 2    format(/,5x,'integration region = ',i3,1x,'number points = ',i4,
     1       /,5x,'starting value     = ',e15.8,1x,'ending value = ',
     2             e15.8)
 3    format(/,1x,'grid summary:')
 4    format(/,'total points = ',i5,/,
     1         'total weights = ',i5)
 5    format(/,'unique number of points = ',i5)
      end

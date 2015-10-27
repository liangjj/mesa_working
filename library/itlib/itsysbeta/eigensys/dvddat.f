*deck dvddat.f
c***begin prologue     dvddat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            data entry for davidson routine.
c***                   
c***references         
c
c***routines called    
c***end prologue       dvdddat
      subroutine dvddat(ops,nwks,nwksg,nroots,nattim,cnverg,thresh,
     1                  mxiter,left,prdvd,dvdall)
      implicit integer (a-z)
      real*8 cnverg, thresh, fpkey
      character*(*) ops
      logical logkey, prdvd, dvdall
      dimension prdvd(11)      
      common/io/inp, iout
      thresh=fpkey(ops,'ci=tolerance',1.0d-10,' ')
      cnverg=fpkey(ops,'ci=convergence',1.0d-08,' ')
      nattim=intkey(ops,'ci=nroots-at-a-time',nroots,' ')
      nwksg=intkey(ops,'ci=guess-size',0,' ')
c
c
c         ----- core allocation for guess routine -----
c
      if (nwksg.le.0) then
          do 10 nwksg=10,500,10
            if (wpadti(nwksg*nwksg+nwksg+max(5*nwksg,nwks)) + 
     1                 nwks.gt.left) go to 20
 10           continue
              nwksg=500
              go to 30
 20           continue
              nwksg=nwksg-10
 30           continue
      endif
c
      nwksg=min(nwks,nwksg)
      if(logkey(ops,'mcscf',.false.,' ')) then
         mxiter=intkey(ops,'mcscf=ci=iterations',20,' ')
      else
         mxiter=intkey(ops,'ci=iterations',15,' ')
      end if
      if(mxiter.lt.2*nattim) then
         mxiter=2*nattim
      endif   
      mxiter=min(mxiter,nwks)
      prdvd(1)=logkey(card,'print=ci=trials',.false.,' ')
      prdvd(2)=logkey(card,'print=ci=vectors',.false.,' ')
      prdvd(3)=logkey(card,'print=ci=h-on-vectors',.false.,' ')
      prdvd(4)=logkey(card,'print=ci=hamiltonian',.false.,' ')
      prdvd(5)=logkey(card,'print=ci=iteration-information',.false.,' ')
      prdvd(6)=logkey(card,'print=ci=residuals',.false.,' ')
      prdvd(7)=logkey(card,'print=ci=transformed-vectors',.false.,' ')
      prdvd(8)=logkey(card,'print=ci=transformed-h-on-vectors',
     1                               .false.,' ')               
      prdvd(9)=logkey(card,'print=ci=new-trial-vectors',.false.,' ')
      prdvd(10)=logkey(card,'print=ci=overlaps',.false.,' ')
      dvdall=logkey(card,'print=ci=all',.false.,' ')
      prdvd(11)=.false.
      if(dvdall) then
         do 50 i=1,11
            prdvd(i)=.true.
 50      continue
      endif
      write(iout,1) nroots, nattim, nwks, nwksg, mxiter, 
     1              thresh, cnverg
 1    format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots           = ',i4,/,5x,
     2             'number of roots at a time = ',i4,/,5x,
     3             'size of matrix            = ',i8,/,5x,
     4             'size of guess matrix      = ',i7,/,5x,
     5             'max. number of iterations = ',i6,/,5x,
     6             'overlap tolerance         = ',e15.8,/,5x,
     7             'convergence criterion     = ',e15.8)

      return
      end





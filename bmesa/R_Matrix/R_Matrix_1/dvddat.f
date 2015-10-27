*deck dvdddat.f
c***begin prologue     dvdddat
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
      subroutine dvddat(card,cpass,n,nroots,ntrials,nattim,cnverg,
     1                  thresh,niter,nvec,lenbuf,prdvd,cgrid,hamd,
     2                  n0,filham)
      implicit integer (a-z)
      real*8 cnverg, thresh, fpkey
      character*128 filham
      character*(*) card, cpass
      logical dollar, logkey, cgrid, hamd, prdvd
      dimension prdvd(*)      
      common/io/inp, iout
      if(dollar('$davidson',card,cpass,inp)) then
         nroots=intkey(card,'number-of-roots',1,' ')
         nroots=min(nroots,n)
         ntrials=intkey(card,'number-of-trial-vectors',nroots,' ')
         cnverg=fpkey(card,'convergence',1.d-08,' ')
         thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
         niter=intkey(card,'maximum-number-of-iterations',
     1                n,' ')
         nattim=intkey(card,'number-of-roots-at-a-time',
     1                 nroots,' ')
         nvec=intkey(card,'maximum-number-of-davidson-vectors',
     1               n,' ')
         if(nvec.lt.2*nattim) then
            nvec=2*nattim
         endif   
         nvec=min(nvec,n)
         lenbuf=intkey(card,'hamiltonian-buffer',min(1000000,n*n),' ')
         cgrid=logkey(card,'use-coarse-grid-vectors',.false.,' ')
         hamd=logkey(card,'h0-preconditioning',.false.,' ')
         n0=intkey(card,'to-h0=nkeep',n,' ')
         prdvd(1)=logkey(card,'print=davidson=all',.false.,' ')
         prdvd(2)=logkey(card,'print=davidson=trials',.false.,' ')
         prdvd(3)=logkey(card,'print=davidson=vectors',
     1                   .false.,' ')
         prdvd(4)=logkey(card,'print=davidson=h-on-vectors',
     1                   .false.,' ')
         prdvd(5)=logkey(card,'print=davidson=hamiltonian',
     1                   .false.,' ')
         prdvd(6)=logkey(card,'print=davidson=iteration-'//
     1                        'information',.false.,' ')
         prdvd(7)=logkey(card,'print=davidson=residuals',
     1                   .false.,' ')
         prdvd(8)=logkey(card,'print=davidson=transformed-vectors',
     1                  .false.,' ')
         prdvd(9)=logkey(card,'print=davidson=transformed-h-'//
     1                   'on-vectors',.false.,' ')               
         prdvd(10)=logkey(card,'print=davidson=new-trial-vectors',
     1                   .false.,' ')
         prdvd(11)=logkey(card,'print=davidson=overlaps',
     1                   .false.,' ')
         prdvd(12)=logkey(card,'print=davidson=buffered-hamiltonian',
     1                    .false.,' ') 
         if(prdvd(1)) then
            do 10 i=2,12
               prdvd(i)=.true.
 10         continue
         endif
         call iosys('write integer "buffer size" to ham',1,
     1               lenbuf,0,' ')           
         call iosys('write integer "davidson print options" to ham',
     1               12,prdvd,0,' ')         
      endif            
      return
      end


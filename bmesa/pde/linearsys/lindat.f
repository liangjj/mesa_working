*deck lindat.f
c***begin prologue     lindat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            data entry for iterative solver routine.
c***                   
c***references         
c
c***routines called    
c***end prologue       lindat
      subroutine lindat(card,cpass,n3d,nrhs,nattim,cnverg,thresh,
     1                  niter,nvec,lenbuf,cgrid,prbufh,hamd,
     2                  prlin,prital,filham)
      implicit integer (a-z)
      real*8 cnverg, thresh, fpkey
      character*128 filham
      character*(*) card, cpass
      logical dollar, logkey, cgrid, prbufh, hamd, prital, prlin
      dimension prlin(11)
      common/io/inp, iout
      if(dollar('$itsolve',card,cpass,inp)) then
         cnverg=fpkey(card,'convergence',1.d-08,' ')
         thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
         niter=intkey(card,'maximum-number-of-iterations',
     1                n3d,' ')
         nattim=intkey(card,'number-of-solutions-at-a-time',
     1                 nrhs,' ')
         nvec=intkey(card,'maximum-number-of-iteration-vectors',
     1               2*nrhs,' ')
         if(nvec.lt.2*nrhs) then
            nvec=2*nrhs
         endif   
c            nvec=min(nvec,n3d)
         lenbuf=intkey(card,'hamiltonian-buffer',
     1                 min(1000000,n3d*n3d),' ')
         cgrid=logkey(card,'use-coarse-grid-vectors',.false.,' ')
         prbufh=logkey(card,'print=buffered-hamiltonian',.false.,' ') 
         hamd=logkey(card,'h0-preconditioning',.false.,' ')
         prlin(1)=logkey(card,'print=iterative-solve=trials',
     1                   .false.,' ')
         prlin(2)=logkey(card,'print=iterative-solve=vectors',
     1                   .false.,' ')
         prlin(3)=logkey(card,'print=iterative-solve=h-on-vectors',
     1                   .false.,' ')
         prlin(4)=logkey(card,'print=iterative-solve=hamiltonian',
     1                   .false.,' ')
         prlin(5)=logkey(card,'print=iterative-solve=iteration-'//
     1                        'information',.false.,' ')
         prlin(6)=logkey(card,'print=iterative-solve=residuals',
     1                   .false.,' ')
         prlin(7)=logkey(card,'print=iterative-solve='//
     1                        'transformed-vectors',.false.,' ')
         prlin(8)=logkey(card,'print=iterative-solve='//
     1                        'transformed-h-on-vectors',
     2                         .false.,' ')               
         prlin(9)=logkey(card,'print=iterative-solve='//
     1                        'new-trial-vectors',.false.,' ')
         prlin(10)=logkey(card,'print=iterative-solve=overlaps',
     1                   .false.,' ')
         prital=logkey(card,'print=iterative-solve=all',.false.,' ')
         if(prital) then
            do 60 i=1,11
               prlin(i)=.true.
 60         continue
         endif
      endif             
      call iosys ('read character "hamiltonian filename" from rwf',
     1             -1,0,0,filham)
      call iosys('open ham as new',0,0,0,filham)
      call iosys('write integer "buffer size" to ham',1,
     1            lenbuf,0,' ')           
      call iosys('write integer "iterative solve print '//
     1           'options" to ham',11,prlin,0,' ')         
      return
      end







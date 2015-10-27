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
      subroutine lindat(card,cpass,n,cnverg,thresh,eps,precon,nblck,
     1                  lenbuf,prlin,filham)
      implicit integer (a-z)
      real*8 cnverg, thresh, fpkey, eps
      character*128 filham
      character*(*) card, cpass, precon
      character*80 chrkey
      logical dollar, logkey, prlin
      dimension prlin(*)
      common/io/inp, iout
      if( dollar('$gmres',card,cpass,inp) ) then 
         cnverg=fpkey(card,'convergence',1.d-08,' ')
         thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
         eps=fpkey(card,'restart-criterion',cnverg,' ')
         precon=chrkey(card,'preconditioner','none',' ')
         lenbuf=intkey(card,'hamiltonian-buffer',min(1000000,n*n),' ')
         nblck=0
         if(precon.ne.'none') then
            nblck=intkey(card,'maximum-size-of-preconditioning-block',
     1                         200,' ')                      
         endif
         prlin(1)=logkey(card,'print=iterative-solve=all',.false.,' ')
         prlin(2)=logkey(card,'print=iterative-solve=trials',
     1                  .false.,' ')
         prlin(3)=logkey(card,'print=iterative-solve=vectors',
     1                   .false.,' ')
         prlin(4)=logkey(card,'print=iterative-solve=h-on-vectors',
     1                   .false.,' ')
         prlin(5)=logkey(card,'print=iterative-solve=hamiltonian',
     1                   .false.,' ')
         prlin(6)=logkey(card,'print=iterative-solve=iteration-'//
     1                   'information',.false.,' ')
         prlin(7)=logkey(card,'print=iterative-solve=residuals',
     1                   .false.,' ')
         prlin(8)=logkey(card,'print=iterative-solve='//
     1                   'transformed-vectors',.false.,' ')
         prlin(9)=logkey(card,'print=iterative-solve='//
     1                   'transformed-h-on-vectors',
     2                   .false.,' ')               
         prlin(10)=logkey(card,'print=iterative-solve='//
     1                    'new-trial-vectors',.false.,' ')
         prlin(11)=logkey(card,'print=iterative-solve=overlaps',
     1                .false.,' ')
         if(prlin(1)) then
            do 10 i=2,11
               prlin(i)=.true.
 10         continue
         endif
      endif
      call iosys('write integer "iterative solve print '//
     1           'options" to ham',11,prlin,0,' ')         
      return
      end

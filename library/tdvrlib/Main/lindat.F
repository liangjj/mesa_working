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
      subroutine lindat(card,cpass,cnverg,thresh,eps,precon,nblck,
     1                  prn,filham,type)
      implicit integer (a-z)
      real*8 cnverg, thresh, fpkey, eps
      character*128 filham
      character*(*) card, cpass, precon, type
      character*80 chrkey, prnkey
      logical dollar, logkey, prn
      dimension prnkey(11)
      common/io/inp, iout
      data prnkey / 'trials', 'vectors', 'h-on-vectors',
     #              'hamiltonian', 'iteration-information',
     #              'residuals', 'transformed-vectors',
     #              'transformed-h-on-vectors',
     #              'new-trial-vectors', 'overlaps', 'all' /
     #
      cnverg=fpkey(card,'convergence',1.d-08,' ')
      thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
      eps=fpkey(card,'restart-criterion',cnverg,' ')
      precon=chrkey(card,'preconditioner','none',' ')
      nblck=0
      if(precon.ne.'none') then
         nblck=intkey(card,'maximum-size-of-preconditioning-block',
     1                      200,' ')                      
      endif
      do 10 i=1,10
         prnkey(i)='print=iterative-solve='//prnkey(i)
 10   continue
      prnkey(11)=chrkey(card,'print=iterative-solve=',prnkey(11),' ')
      if(prnkey(11).eq.'all') then
         call setprn(prn,10)
      else
         call setlog(prn,prnkey,card,11)
      endif    
      call iosys ('read character "hamiltonian filename" from rwf',
     1             -1,0,0,filham)
      call iosys('open ham as '//type,0,0,0,filham)
      call iosys('write integer "iterative solve print '//
     1           'options" to ham',11,prlin,0,' ')         
      return
      end














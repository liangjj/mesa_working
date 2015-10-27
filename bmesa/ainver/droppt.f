*deck droppt.f
c***begin prologue     droppt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       droppt
      subroutine droppt(q,wt,eigc,wtc,nfix,fl,fr,nin,nout)
      implicit integer (a-z)
      real*8 q, wt, eigc, wtc
      logical nfix
      dimension q(nin), wt(nin), eigc(*), wtc(*)
      dimension nfix(2)
      common/io/inp, iout
      nout=nin
      if(fl.ne.0.and.fr.ne.0) then
         write(iout,1)
         call copy(q,eigc,nin)
         call copy(wt,wtc,nin)
         return
      endif
      if(fl.eq.0.and.fr.ne.0) then
         if(nfix(1)) then
            write(iout,2)
            nout=nout-1
c
c           drop first function and first component
c
            call copy(q(2),eigc(1),nout)
            call copy(wt(2),wtc(1),nout)
            return
         endif
      endif
      if(fl.ne.0.and.fr.eq.0) then
         if(nfix(2)) then
            nout=nout-1
            write(iout,3)
c
c           drop last function and last component
c
            call copy(q(1),eigc(1),nout)
            call copy(wt(1),wtc(1),nout)
            return
         endif
      endif 
      if(fl.eq.0.and.fr.eq.0) then
         if(.not.nfix(1).and..not.nfix(2)) then
            write(iout,1)
            call copy(q,eigc,nin)
            call copy(wt,wtc,nin)
            return
         elseif(nfix(1).and..not.nfix(2)) then 
            write(iout,2)
            nout=nout-1
c
c           drop first function and first component
c
            call copy(q(2),eigc(1),nout)
            call copy(wt(2),wtc(1),nout)
            return
         elseif(.not.nfix(1).and.nfix(2)) then
            nout=nout-1
            write(iout,3)
c
c           drop last function and last component
c
            call copy(q(1),eigc(1),nout)
            call copy(wt(1),wtc(1),nout)
            return
c
c        drop first and last function and component
c 
         elseif(nfix(1).and.nfix(2)) then
            nout=nout-2
            call copy(q(2),eigc(1),nout)
            call copy(wt(2),wtc(1),nout)
         else
            call lnkerr('error in call to drop')
         endif
      endif
 1    format(/,5x,'neither boundary a nodal point',/,5x,
     1            '             or ',/,5x,
     2            'basis has nodal structure built in.',/,5x,
     3            '      no functions dropped')
 2    format(/,5x,'dropping function at left endpoint')
 3    format(/,5x,'dropping function at right endpoint')
 4    format(/,5x,'dropping function at left and right endpoint')

      end       

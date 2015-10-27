*deck dropfn.f
c***begin prologue     dropfn
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
c***end prologue       dropfn
      subroutine dropfn(pppji,fl,fr,nj,ni,mj,mi)
      implicit integer (a-z)
      real*8 p, dp, ddp, pn, dpn, ddpn
      logical nfix, nodrop
      dimension p(nj,ni), dp(nj,ni), ddp(nj,ni)
      dimension pn(mj,mi), dpn(mj,mi), ddpn(mj,mi)
      dimension nfix(2)
      common/io/inp, iout
      if(fl.ne.0.and.fr.ne.0) then
         write(iout,1)
         call copy(p,pn,nj*ni)
         call copy(dp,dpn,nj*ni)
         call copy(ddp,ddpn,nj*ni)
         return
      endif
      if(fl.eq.0.and.fr.ne.0) then
         if(nfix(1)) then
            write(iout,2)
c
c           drop first function and first component
c
            do 10 i=2,ni
               do 20 j=2,nj
                  pn(j-1,i-1)=p(j,i)
                  dpn(j-1,i-1)=dp(j,i)
                  ddpn(j-1,i-1)=ddp(j,i)
 20            continue
 10         continue   
            return
         endif
      endif
      if(fl.ne.0.and.fr.eq.0) then
         if(nfix(2)) then
            write(iout,3)
c
c           drop last function and last component
c
            do 30 i=1,ni-1
               do 40 j=1,nj-1
                  pn(j,i)=p(j,i)
                  dpn(j,i)=dp(j,i)
                  ddpn(j,i)=ddp(j,i)
 40            continue
 30         continue   
            return
         endif
      endif 
      if(fl.eq.0.and.fr.eq.0) then
         if(.not.nfix(1).and..not.nfix(2)) then
            write(iout,1)
            call copy(p,pn,ni*nj)
            call copy(dp,dpn,ni*nj)
            call copy(ddp,ddpn,ni*nj)
            return
         elseif(nfix(1).and..not.nfix(2)) then 
            write(iout,2)
c
c           drop first function and first component
c
            do 50 i=2,ni
               do 60 j=2,nj
                  pn(j-1,i-1)=p(j,i)
                  dpn(j-1,i-1)=dp(j,i)
                  ddpn(j-1,i-1)=ddp(j,i)
 60            continue
 50         continue   
            return
         elseif(.not.nfix(1).and.nfix(2)) then
            write(iout,3)
c
c           drop last function and last component
c
            do 70 i=1,ni-1
               do 80 j=1,nj-1
                  pn(j,i)=p(j,i)
                  dpn(j,i)=dp(j,i)
                  ddpn(j,i)=ddp(j,i)
 80            continue
 70         continue   
            return
c
c        drop first and last function and component
c 
         elseif(nfix(1).and.nfix(2)) then
            nout=nout-2
            do 100 i=2,ni-1
               do 200 j=2,nj-1
                  pn(j-1,i-1)=p(j,i)
                  dpn(j-1,i-1)=dp(j,i)
                  ddpn(j-1,i-1)=ddp(j,i)
 200           continue
 100        continue   
         else
            call lnkerr('error in call to dropfn')
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

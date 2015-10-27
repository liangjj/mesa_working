*deck vparm.f
c***begin prologue     vparm
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential parameters
c***                   
c***description        parameters to construct various standard potentials
c***                   
c***references         
c
c***routines called    
c***end prologue       vparm
      subroutine vparm(type,a,b,card,cpass,key,dim)
      implicit integer (a-z)
      real*8 a, b, pi, hbar, massau, lenau, timau, pmass, massn2p
      real*8 fpkey
      character*(*) type, card, cpass, key
      character*80 chrkey, units
      character*1 itoc
      character*24 phr
      logical dollar      
      dimension a(dim,dim), b(dim,dim), phr(4)
      common/io/inp, iout
      data pi  / 3.141592653589793238462643d+00 /
c     hbar in joule-sec
      data hbar/1.054571596d-34/                                  
      data massau, lenau, timau, pmass / 9.10938188d-31, 
     1                                   5.291772083d-11,
     2                                   2.418884326d-17, 
     3                                   1.67262158d-27 /
      data massn2p / 1.00137841887d0 / 
      if(dollar(key,card,cpass,inp) ) then
         type=chrkey(card,'potential','none',' ')
         write(iout,1) type
      endif
      if(type.eq.'none') then
         return
      elseif(type.eq.'well') then
         do 10 i=1,dim
            phr(1)='d'//itoc(i)
            phr(2)='l'//itoc(i)
            do 20 j=1,dim
               phr(3)=phr(1)(1:2)//itoc(j)
               phr(4)=phr(2)(1:2)//itoc(j)
               a(i,j)=fpkey(card,phr(3)(1:3),0.d0,' ')
               b(i,j)=fpkey(card,phr(4)(1:3),0.d0,' ')
 20         continue
 10      continue
         do 30 i=1,dim
            do 40 j=1,dim   
               if(a(i,j).ne.a(j,i).or.b(i,j).ne.b(j,i)) then
                  call lnkerr('error in off-diagonal well parameters')
               endif
 40         continue
 30      continue   
         cpass='well depth parameters'
         call prntrm(cpass,a,dim,dim,dim,dim,iout)
         cpass='well length parameters'
         call prntrm(cpass,b,dim,dim,dim,dim,iout)
      elseif(type.eq.'exponential') then
         do 50 i=1,dim
            phr(1)='a'//itoc(i)
            phr(2)='b'//itoc(i)
            do 60 j=1,dim
               phr(3)=phr(1)(1:2)//itoc(j)
               phr(4)=phr(2)(1:2)//itoc(j)
               a(i,j)=fpkey(card,phr(3)(1:3),0.d0,' ')
               b(i,j)=fpkey(card,phr(4)(1:3),1.d0,' ')
 60         continue
 50      continue   
         do 70 i=1,dim
            do 80 j=1,dim
               if(a(i,j).ne.a(j,i).or.b(i,j).ne.b(j,i)) then
                  call lnkerr('error in off-diagonal exponential '//
     1                        'parameters')
               endif
 80         continue
 70      continue   
         cpass='exponential strength parameters'
         call prntrm(cpass,a,dim,dim,dim,dim,iout)
         cpass='exponential range parameters'
         call prntrm(cpass,b,dim,dim,dim,dim,iout)
      elseif(type.eq.'coulomb') then
         do 90 i=1,dim
            phr(1)='z'//itoc(i)
            do 100 j=1,dim
               phr(2)=phr(1)(1:2)//itoc(j)
               a(i,j)=fpkey(card,phr(2)(1:3),-1.d0,' ')
 100        continue
 90      continue   
         do 110 i=1,dim
            do 120 j=1,dim
               if(a(i,j).ne.a(j,i)) then
                  call lnkerr('error in off-diagonal coulomb '//
     1                        'parameters')
               endif
 120        continue
 110     continue   
      elseif(type.eq.'harmonic-oscillator') then
         units=chrkey(card,'units','atomic-units',' ')
         ln=length(units)
         if(units(1:ln).eq.'atomic-units') then
            do 200 i=1,dim
               phr(1)='omega'//itoc(i)
               phr(2)='mass'//itoc(i)
               do 210 j=1,dim
                  phr(3)=phr(1)(1:6)//itoc(j)
                  phr(4)=phr(2)(1:5)//itoc(j)
                  b(i,j)=fpkey(card,phr(3)(1:7),1.d0,' ')
                  a(i,j)=fpkey(card,phr(4)(1:6),1.d0,' ')
 210           continue
 200        continue
         else
            do 300 i=1,dim
               phr(1)='omega'//itoc(i)
               phr(2)='number-of-protons'//itoc(i)
               phr(3)='number-of-neutrons'//itoc(i)
               do 310 j=1,dim            
                  phr(4)=phr(1)(1:6)//itoc(j)
                  phr(5)=phr(2)(1:18)//itoc(j)
                  phr(6)=phr(2)(1:19)//itoc(j)
                  b(i,j) = fpkey(card,phr(5)(1:7),1.d0,' ')
                  b(i,j)=b(i,j)*2.d0*pi
                  a(i,j) = pmass * 
     1                     fpkey(card,phr(5)(1:19),1,' ') +
     4                     pmass * fpkey(ops,phr(6)(1:20),0,' ') *
     7                     massn2p

 310           continue
 300        continue   
         endif
      else
         call lnkerr('error in potential type')
      endif
      return
 1    format(/,5x,'interaction potential = ',a32)
      end       






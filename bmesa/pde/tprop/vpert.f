*deck vpert.f
c***begin prologue     vpert
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            time-dependent potential
c***                   
c***description        calculate the time and space dependent potential
c***                   matrix elements in the dvr representation.
c***                   
c***references         
c
c***routines called    
c***end prologue       vpert
      subroutine vpert(v,eig1,eig2,eig3,t,efield,omega,width,shift,
     1                 nd,dim,type)
      implicit integer (a-z)
      real*8 v, eig1, eig2, eig3, t, efield, omega, width, shift
      real*8 dt, pi
      character*(*) type
      character*8 cpass
      logical posinp
      dimension nd(3)      
      dimension v(*), eig1(*), eig2(*), eig3(*)
      common/io/inp, iout
      data pi/3.1415926535897932384d0/
      count = 0
      if(type.eq.'t') then
         if(dim.eq.1) then
            do 10 i=1,nd(1)
               v(i) = t                  
 10         continue
         endif
         if(dim.eq.2) then
            do 20 i=1,nd(1)
               do 30 j=1,nd(2)
                  count = count + 1
                  v(count) = t                     
 30            continue
 20         continue
         endif
         if(dim.eq.3) then
            do 40 i=1,nd(1)
               do 50 j=1,nd(2)
                  do 60 k=1,nd(3)
                     count = count + 1
                     v(count) = t                     
 60               continue
 50            continue
 40         continue   
         endif   
      elseif(type.eq.'cosine') then
         dt = efield*cos(omega*t)
         if(dim.eq.1) then
            do 100 i=1,nd(1)
               v(i) = dt*eig1(i)
 100        continue
         endif   
         if(dim.eq.2) then
            do 110 i=1,nd(1)
               do 120 j=1,nd(2)
                  count = count + 1
                  v(count) = dt*( eig1(i) + eig2(j) )
 120              continue
 110           continue
         endif
         if(dim.eq.3) then
            do 130 i=1,nd(1)
               do 140 j=1,nd(2)
                  do 150 k=1,nd(3)
                     count = count + 1
                     v(count) = dt*( eig1(i) + eig2(j) +eig3(k) )
 150                 continue
 140              continue
 130           continue   
         endif   
      elseif(type.eq.'gaussian-pulse') then
         dt=efield * cos( omega * t ) * exp( - width *
     1                  ( t - shift ) * ( t - shift ) )
         if(dim.eq.1) then
            do 200 i=1,nd(1)
               v(i) = dt*eig1(i)
 200        continue
         endif
         if(dim.eq.2) then
            do 210 i=1,nd(1)
               do 220 j=1,nd(2)
                  count = count + 1
                  v(count) = dt*( eig1(i) + eig2(j) )
 220           continue
 210        continue
         endif   
         if(dim.eq.3) then
            do 230 i=1,nd(1)
               do 240 j=1,nd(2)
                  do 250 k=1,nd(3)
                     count = count + 1
                     v(count) = dt*( eig1(i) + eig2(j) + eig3(k) )
 250              continue
 240           continue
 230        continue
         endif
      endif   
      return
      end       



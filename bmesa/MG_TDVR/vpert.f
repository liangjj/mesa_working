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
      subroutine vpert(v,vt,eig1,eig2,eig3,eigt,n,nd,nt,dim,type,type0)
      implicit integer (a-z)
      real*8 v, vt, eig1, eig2, eig3, eigt, omega, scale, width, shift
      real*8 dt, pi, fpkey, scr
      character*(*) type, type0
      character*1600 card
      character*80 cpass, units, chrkey
      logical dollar
      dimension nd(3)      
      dimension v(n*nt), eig1(nd(1)), eig2(nd(2)), eig3(nd(3)), eigt(nt)
      dimension vt(nt)
      common/io/inp, iout
      data pi/3.1415926535897932384d0/
      pointer(pscr,scr(1))
      write(iout,1) 
      call rzero(v,n*nt)
      if(type.ne.'none') then
         call memory(nt,pscr,ngot,'scr',0)
      endif  
      if(dollar('$vt',card,cpass,inp) ) then
         units=chrkey(card,'units','atomic-units',' ')
         omega=fpkey(card,'electric-field-frequency',10.d0,' ')
         scale=fpkey(card,'electric-field-strength',1.d0,' ')
         if(type.eq.'cosine-gaussian-pulse'.or.
     1      type.eq.'sine-gaussian-pulse') then
            width=fpkey(card,'width-of-gaussian-pulse',0.d0,' ')
            shift=fpkey(card,'shift-of-gaussian-pulse',0.d0,' ')
            scale=fpkey(card,'scale-of-gaussian-pulse',1.d0,' ')
         endif
         write(iout,2) units
         if(units.eq.'hertz') then
            omega=2.d0*pi*omega   
         endif
         if(type.eq.'cosine'.or.
     1      type.eq.'cosine-gaussian-pulse'.or.
     2      type.eq.'cosine-dipole-field') then
            do 10 i=1,nt
               scr(i)=scale*cos(omega*eigt(i))
 10         continue  
         endif
         if(type.eq.'sine'.or.
     1      type.eq.'sine-gaussian-pulse'.or.
     2      type.eq.'sine-dipole-field') then
            do 20 i=1,nt
               scr(i)=scale*sin(omega*eigt(i))
 20         continue  
         endif
      endif
      if(type.eq.'cosine'.or.type.eq.'sine') then
         if(type.eq.'cosine') then
            write(iout,3) omega, scale
         else
            write(iout,4) omega, scale
         endif
         if(dim.eq.1) then
            do 100 i=1,nt
               do 110 j=1,nd(1)
                  count=count + 1     
                  v(count) = scr(i)
 110           continue   
 100        continue
         elseif(dim.eq.2) then
            do 120 i=1,nt
               do 130 j=1,nd(1)
                  do 140 k=1,nd(2)
                     count = count + 1
                     v(count) = scr(i)
 140              continue   
 130           continue
 120        continue
         elseif(dim.eq.3) then
            do 150 i=1,nt
               do 160 j=1,nd(1)
                  do 170 k=1,nd(2)
                     do 180 l=1,nd(3)
                        count = count + 1
                        v(count) = scr(i)
 180                 continue
 170              continue   
 160           continue
 150        continue
         endif   
      elseif(type.eq.'cosine-dipole-field'.or.
     1       type.eq.'sine-dipole-field') then
         if(dim.ne.1) then
            call lnkerr('cannot handle larger than 1-D case')
         endif
         if(type.eq.'cosine-dipole-field') then
            write(iout,5) omega, scale
         else
            write(iout,6) omega, scale
         endif
         count=0
         do 190 i=1,nt
            do 200 j=1,nd(1)
               count=count+1
               v(count)=scr(i)*eig1(j)
 200        continue
 190     continue
      elseif(type.eq.'cosine-gaussian-pulse'.or.
     1       type.eq.'sine-gaussian-pulse' ) then
         if(type.eq.'cosine-gaussian-pulse') then
            write(iout,7) scale, width, shift, omega     
         else
            write(iout,8) scale, width, shift, omega     
         endif
         if(dim.eq.1) then
c
c              assume the polarization is along the only coordinate axis
c
            do 300 i=1,nt
               dt = scr(i) * exp( - width *
     1                       ( eigt(i) - shift ) * 
     2                       ( eigt(i) - shift ) )
               do 310 j=1,nd(1)
                  count=count+1
                  v(count) = dt*eig1(j)
 310           continue
 300        continue   
         elseif(dim.eq.2) then
c
c              assume polarization in (x,y) plane, i.e. along rho
c
            do 320 i=1,nt
               dt = scr(i) * exp( - width *
     1                          ( eigt(i) - shift ) * 
     2                          ( eigt(i) - shift ) )
               do 330 j=1,nd(1)
                  do 340 k=1,nd(2)
                     count = count + 1
                     v(count) = dt*sqrt( eig1(j)*eig1(j) + 
     1                                   eig2(k)*eig2(k) )
 340              continue    
 330           continue   
 320        continue
         elseif(dim.eq.3) then
c
c              assume polarization is radial
c
            do 350 i=1,nt
               dt = scr(i) * exp( - width *
     1                          ( eigt(i) - shift ) * 
     2                          ( eigt(i) - shift ) )
               do 360 j=1,nd(1)
                  do 370 k=1,nd(2)
                     do 380 l=1,nd(3)
                        count = count + 1
                        v(count) = dt*sqrt( eig1(j)*eig1(j) + 
     1                                      eig2(k)*eig2(k) + 
     2                                      eig3(l)*eig3(l) )
 380                 continue
 370              continue
 360           continue
 350        continue
         endif

      endif
c
c     add in any pure time potential
c
      if(type0.ne.'none') then
         write(iout,9) type0
         count = 0
         do 1000 i=1,nt
            do 2000 j=1,n
               count = count + 1
               v(count) = v(count) + vt(i)
 2000       continue
 1000    continue
      endif
c      write(iout,*) v   
      if(type.ne.'none') then
         call memory(-ngot,pscr,idum,'scr',idum)
      endif 
      return
 1    format(/,5x,'constructing factors for a space and time-dependent '
     1            'perturbation')
 2    format(/,5x,'frequency units = ',a16)
 3    format(/,5x,'perturbation = E0 * cosine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 4    format(/,5x,'perturbation = E0 * sine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 5    format(/,5x,'perturbation = E0 * x * cosine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 6    format(/,5x,'perturbation = E0 * x * sine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 7    format(/,5x,'perturbation = A * x * cos(omega * t) * ',
     1       /,5x,'               exp(- W * (t-S) * (t-S))',
     2       /,5x,'    A        = ',e15.8,
     3       /,5x,'    omega    = ',e15.8,
     4       /,5x,'    W        = ',e15.8,
     5       /,5x,'    S        = ',e15.8)
 8    format(/,5x,'perturbation = A * x * sin(omega * t) * ',
     1       /,5x,'               exp(- W * (t-S) * (t-S))',
     2       /,5x,'    A        = ',e15.8,
     3       /,5x,'    omega    = ',e15.8,
     4       /,5x,'    W        = ',e15.8,
     5       /,5x,'    S        = ',e15.8)
 9    format(/,5x,'adding in a pure time dependent potential = ',a24)
      end       



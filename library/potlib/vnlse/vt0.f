*deck vt0.f
c***begin prologue     vt0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            pure time dependent potential
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vt0
      subroutine vt0(card,t,vt,hbar,units,n,typ)
      implicit integer (a-z)
      real*8 t, vt
      real*8 omega, shift
      real*8 width, scale, pi, hbar, massau, lenau, timau, pmass
      real*8 fpkey
      character*(*) card, typ, units
      dimension t(n), vt(n)
      common/io/inp, iout
      data pi  / 3.141592653589793238462643d+00 /
      data massau, lenau, timau, pmass / 9.10938188d-31, 
     1                                   5.291772083d-11,
     2                                   2.418884326d-17, 
     3                                   1.67262158d-27 /
      if(units.eq.'atomic-units') then
         hbar=1.d0
      endif
      if(typ.eq.'none') then
         return
      elseif(typ.eq.'t') then
         write(iout,1)
         do 10 i=1,n
            vt(i) = t(i)
 10      continue   
      elseif(typ.eq.'cosine') then
         omega=fpkey(card,'electric-field-frequency',10.d0,' ')
         scale=fpkey(card,'scale-of-cosine-pulse',1.d0,' ')
         omega=2.d0*pi*omega   
         if(units.eq.'atomic-units') then
            omega=omega*timau
         endif                      
         write(iout,2) omega, scale
         do 20 i=1,n   
            vt(i) = cos(omega*t(i))
 20      continue
         call vscale(vt,vt,scale,n)
      elseif(typ.eq.'gaussian-pulse') then
         omega=fpkey(card,'electric-field-frequency',10.d0,' ')
         omega=2.d0*pi*omega   
         if(units.eq.'atomic-units') then
            omega=omega*timau
         endif                      
         width=fpkey(card,'width-of-gaussian-pulse',0.d0,' ')
         shift=fpkey(card,'shift-of-gaussian-pulse',0.d0,' ')
         scale=fpkey(card,'scale-of-gaussian-pulse',1.d0,' ')
         write(iout,3) 
         write(iout,4) scale, width, shift, omega     
         do 30 i=1,n
            vt(i)=scale * cos( omega * t(i) ) * exp( - width *
     1                      ( t(i) - shift ) * 
     2                      ( t(i) - shift ) )
 30     continue         
      endif
      return
 1    format(/,5x,'constructing factors for a time-dependent '
     1            'perturbation = t')
 2    format(/,5x,'constructing factors for a time-dependent '
     1            'perturbation = A * cosine(omega*t)',/,5x,
     2 '           omega = ',e15.8,/,5x,
     3 '           A     = ',e15.8)
 3    format(/,5x,'the form of the potential is:',///,10x,
     1            'v = scale * cos(omega*t) * '
     2            'exp(-width*(t-shift)*(t-shift))')
 4    format(/,5x,'constructing factors for a time-dependent '
     1            'perturbation = A * cos(omega*t) * '
     2            '               exp(-W*(t-S)*(t-S))',/,5x,
     3 '                          A        = ',e15.8,/,5x,
     4 '                          omega    = ',e15.8,/,5x,
     5 '                          W        = ',e15.8,/,5x,
     6 '                          S        = ',e15.8)
      end       



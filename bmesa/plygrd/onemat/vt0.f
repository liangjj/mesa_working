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
      subroutine vt0(eig,v,n,type,typ,dim)
      implicit integer (a-z)
      real*8 eig, omega, shift
      real*8 width, scale, v, pi, fpkey
      character*240 card
      character*(*) type, typ
      character*1 itoc
      character*80 cpass 
      logical dollar     
      dimension eig(n), v(n)
      common/io/inp, iout
      data pi/3.1415926535897932384d0/
      lst=cskipb(typ,' ') 
      call rzero(v,n)
      if(type.eq.'t') then
         write(iout,1)
         do 10 i=1,n
            v(i) = eig(i)
 10      continue   
      elseif(type.eq.'cosine') then
         if(dollar('$vt0-'//typ(1:lst)//'-'//itoc(dim),card,
     1              cpass,inp) ) then
            omega=fpkey(card,'electric-field-frequency',10.d0,' ')
            scale=fpkey(card,'scale-of-cosine-pulse',1.d0,' ')
            omega=2.d0*pi*omega   
            write(iout,2) omega, scale
            do 20 i=1,n   
               v(i) = cos(omega*eig(i))
 20         continue
            call vscale(v,v,scale,n)
         endif
      elseif(type.eq.'gaussian-pulse') then
         if(dollar('$vt0-'//typ(1:lst)//'-'//itoc(dim),card,
     1              cpass,inp) ) then
            omega=fpkey(card,'electric-field-frequency',10.d0,' ')
            omega=2.d0*pi*omega   
            width=fpkey(card,'width-of-gaussian-pulse',0.d0,' ')
            shift=fpkey(card,'shift-of-gaussian-pulse',0.d0,' ')
            scale=fpkey(card,'scale-of-gaussian-pulse',1.d0,' ')
            write(iout,3) 
            write(iout,4) scale, width, shift, omega     
            do 30 i=1,n
               v(i)=scale * cos( omega * eig(i) ) * exp( - width *
     1                         ( eig(i) - shift ) * 
     2                         ( eig(i) - shift ) )
 30        continue         
         endif
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



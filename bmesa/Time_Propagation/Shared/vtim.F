c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Vtim}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck vtim.f
c***begin prologue     vtim
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           1-d hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate one body time potential.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vtim
      subroutine vtim(vt,time,key,typ)
      implicit integer (a-z)
      real*8 vt, time
      real*8 omega, shift
      real*8 width, scale, pi, hbar, massau, lenau, timau, pmass
      real*8 fpkey
      character*(*) key, typ 
      character*80 loctyp
      character*80 title, cpass, units, chrkey
      character*240 card
      logical prn, dollar, reuse
      common/io/inp, iout
      save loctyp, units
      data first / 0 /
      data hbar/1.054571596d-34/
      data pi  / 3.141592653589793238462643d+00 /
      data massau, lenau, timau, pmass / 9.10938188d-31, 
     1                                   5.291772083d-11,
     2                                   2.418884326d-17, 
     3                                   1.67262158d-27 /
      if(first.eq.0) then
         if(dollar(key,card,cpass,inp)) then
            typ=chrkey(card,'potential','none',' ')
            loctyp=typ
            units=chrkey(card,'units','atomic-units',' ') 
         endif
         first=1
      endif
      if(units.eq.'atomic-units') then
         hbar=1.d0
      endif
      if(typ.eq.'none') then
         vt=0.d0
         return
      elseif(typ.eq.'t') then
         write(iout,1)
         vt = time
      elseif(typ.eq.'cosine') then
         omega=fpkey(card,'electric-field-frequency',10.d0,' ')
         scale=fpkey(card,'scale-of-cosine-pulse',1.d0,' ')
         omega=2.d0*pi*omega   
         if(units.eq.'atomic-units') then
            omega=omega*timau
         endif                      
         write(iout,2) omega, scale
         vt = scale*cos(omega*time)
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
         vt = scale * cos( omega * time ) * exp( - width *
     1                      ( time - shift ) * 
     2                      ( time - shift ) )
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





c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Vcple}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck vcple.f
***begin prologue     vcple
***date written       960723   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           time development
***author             schneider, barry (nsf)
***source             
***purpose            time-dependent potential
***                   
***description        input data for coupling potential, including the
***                   mean field interaction for a contact potential.
***references         
***routines called    
***end prologue       vcple 
      subroutine vcple(units,typ,hbar,scale,omega,width,shift,gamma,prn)
      implicit integer (a-z)
      real*8 omega, scale, width, shift, gamma
      real*8 fpkey, mass, sctlen
      real*8 pi, hbar, massau, lenau, timau, pmass, massn2p
      character*80 typ
      character*1600 card
      character*80 cpass, units, chrkey
      character*2 atom
      logical dollar, prn
      dimension typ(2)      
      common/io/inp, iout
      data pi  / 3.141592653589793238462643d+00 /
      data massau, lenau, timau, pmass / 9.10938188d-31, 
     1                                   5.291772083d-11,
     2                                   2.418884326d-17, 
     3                                   1.67262158d-27 /
      data massn2p / 1.00137841887d0 / 
      if(dollar('$vnlse',card,cpass,inp) ) then
         units=chrkey(card,'units','atomic-units',' ')
         typ(1)=chrkey(card,'v-space-time','none',' ')
         typ(2)=chrkey(card,'v-non-linear','none',' ')
         if(typ(2).ne.'none') then
            atom=chrkey(card,'atom','cs',' ')
            sctlen=fpkey(card,'scattering-length',52.d0,' ')
            natom=intkey(card,'number-of-atoms',0,' ')
         endif
         write(iout,1) (typ(i),i=1,2)
         write(iout,2) units
         scale=fpkey(card,'electric-field-strength',1.d0,' ')
         omega=fpkey(card,'electric-field-frequency',10.d0,' ')
         if(typ(1).eq.'cosine-gaussian-pulse'.or.
     1      typ(1).eq.'sine-gaussian-pulse') then
            width=fpkey(card,'width-of-gaussian-pulse',0.d0,' ')
            shift=fpkey(card,'shift-of-gaussian-pulse',0.d0,' ')
            scale=fpkey(card,'scale-of-gaussian-pulse',1.d0,' ')
         endif
         omega=2.d0*pi*omega   
         if(units.eq.'atomic-units') then
            omega=omega*timau
         endif                      
         if(typ(1).eq.'cosine-pulse') then
            write(iout,3) omega, scale
         elseif(typ(1).eq.'sine-pulse') then
            write(iout,4) omega, scale
         elseif(typ(1).eq.'t') then
            write(iout,5)
         elseif(typ(1).eq.'cosine-dipole-field') then
            write(iout,6) omega, scale
         elseif(typ(1).eq.'sine-dipole-field') then
            write(iout,7) omega, scale
         elseif(typ(1).eq.'cosine-gaussian-pulse') then
            write(iout,8) scale, width, shift, omega
         elseif(typ(1).eq.'sine-gaussian-pulse') then
            write(iout,9) scale, width, shift, omega
         endif
         if(natom.ne.0) then
            write(iout,11) atom, natom, sctlen
            if(atom.eq.'cs') then
               mass=2.2d-25
            elseif(atom.eq.'na') then
               mass=3.8176d-26
            endif
            if(units.eq.'atomic-units') then
               hbar=1.d0
               mass=mass/massau
            endif
            gamma = natom*natom*4.d0*pi*hbar*hbar *
     1              sctlen/(mass*mass)
         endif
      endif
      return
 1    format(/,1x,'space and time    = ',a32,/,1x,
     1            'non-linear        = ',a32)
 2    format(/,5x,'units = ',a16)
 3    format(/,5x,'perturbation = E0 * cosine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 4    format(/,5x,'perturbation = E0 * sine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 5    format(/,5x,'perturbation =  t')
 6    format(/,5x,'perturbation = E0 * x * cosine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 7    format(/,5x,'perturbation = E0 * x * sine(omega * t)',
     1       /,5x,'       omega = ',e15.8,
     2       /,5x,'       E0    = ',e15.8)
 8    format(/,5x,'perturbation = A * x * cos(omega * t) * ',
     1       /,5x,'               exp(- W * (t-S) * (t-S))',
     2       /,5x,'    A        = ',e15.8,
     3       /,5x,'    omega    = ',e15.8,
     4       /,5x,'    W        = ',e15.8,
     5       /,5x,'    S        = ',e15.8)
 9    format(/,5x,'perturbation = A * x * sin(omega * t) * ',
     1       /,5x,'               exp(- W * (t-S) * (t-S))',
     2       /,5x,'    A        = ',e15.8,
     3       /,5x,'    omega    = ',e15.8,
     4       /,5x,'    W        = ',e15.8,
     5       /,5x,'    S        = ',e15.8)
 11   format(/,15x,'atomic parameters',/,/,5x,
     1             'atom                    = ',a2,/,5x,
     2             'number of atoms in trap = ',i8,/,5x,
     3             'scattering length       = ',e15.8)
      end       



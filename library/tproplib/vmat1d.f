*deck vmat1d.f
c***begin prologue     vmat1d
c***date written       970797   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            dvr representation of one dimension potential
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vmat1d
      subroutine vmat1d(eigc,v,hbar,mass,scale,n,dim,prnt,ops)
      implicit integer (a-z)
      real*8 eigc, v
      real*8 hbar, mass, scale, zii
      real*8 pi, zero, half, one, two, three, ten
      real*8 timau
      real*8 aii, sii, omegii, vii
      real*8 facii, fpkey
      character*32 pot, chrkey
      character*80 cpass
      character*320 card
      character*(*) ops
      character*1 itoc
      logical prnt, dollar, toau, useau, logkey
      dimension eigc(n), v(n)
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
      data zero, half, one, two, three, ten / 0.d0, .5d0, 1.d0, 2.d0, 
     1                                        3.d0, 10.d0 /
      data timau / 2.418884d-17 /
c
c
c
c     lets assume pairwise interactions of some sort
c
c        if the model potential is a harmonic oscillator the form is;
c
c                    m*omega*omega*r*r
c              v =   _________________
c                           2
c        if the model potential is an exponential the form is:
c
c        here the potential is
c             v = - a*exp(-b*r)
c
c        if the model potential is a coulomb potential the form is:
c        
c             v = z/r
c
c
      if( dollar('$v0-'//itoc(dim),card,cpass,inp) ) then
          pot=chrkey(card,'potential-type','none',' ')
          scale=one
          zii=fpkey(card,'z',-1.d0,' ')
          vii=fpkey(card,'v',zero,' ')
          omegii=fpkey(card,'omega',zero,' ')
          omegii=omegii*two*pi
          if(pot.eq.'harmonic-oscillator') then
             scale=one/(omegii*hbar)
          endif   
          aii=fpkey(card,'a',zero,' ')
          sii=fpkey(card,'s',zero,' ')
      endif
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      if(toau) then
         omegii=omegii*timau
      endif
      if(useau) then
         omegii=one
      endif
      call rzero(v,n)
      if(pot.eq.'harmonic-oscillator') then
         facii=mass*omegii*omegii*half
         write(iout,1) dim, omegii
         do 10 i=1,n
            v(i)=facii*eigc(i)*eigc(i)
 10      continue
      elseif(pot.eq.'exponential') then
         write(iout,2) dim, aii, sii
         do 20 i=1,n
            v(i)=aii*exp(-sii*abs(eigc(i))) 
 20      continue   
      elseif(pot.eq.'well') then
         write(iout,3) dim, vii
         do 30 i=1,n
            v(i) = vii
 30        continue                  
      elseif(pot.eq.'coulomb') then
         write(iout,4) dim, zii
         do 40 i=1,n
            v(i)=zii/eigc(i) 
 40      continue   
      elseif(pot.eq.'inverse-r4') then
         write(iout,5)
         do 50 i=1,n
            v(i)=one/( 1.d0 + eigc(i) )**4
 50      continue   
      endif
      return
 1    format(/,5x,'oscillator frequency:',/,5x,
     1            'dimension = ',i1,/,5x,
     2            'omegaii   = ',e15.8)
 2    format(/,5x,'exponential constant:',/,5x,
     1            'dimension = ',i1,/,5x,
     2            'scale         = ',e15.8,/,5x,
     3            'exponent      = ',e15.8)
 3    format(/,5x,'well constant:',/,5x,
     1            'dimension = ',i1,/,5x,
     2            'constant  = ',e15.8)
 4    format(/,5x,'coulomb charge:',/,5x,
     1            'dimension = ',i1,/,5x,
     2            'z         = ',e15.8)
 5    format(/,5x,'1/r4')
      end       

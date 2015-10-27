*deck vmat.f
c***begin prologue     vmat
c***date written       970797   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            interaction potential in dvr representation
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vmat
      subroutine vmat(eigc,v,n,prnt,ops)
      implicit integer (a-z)
      real*8 eigc
      complex*16 v
      real*8 vr, vi, omegar, omegai, ar, ai, sr, si, zr, zi
      real*8 pi, zero, half, one, two, three, ten
      complex*16 a, s, omega, vfac, z, scale
      real*8 fpkey
      character*32 pot, chrkey
      character*80 cpass
      character*320 card
      character*(*) ops
      logical prnt, dollar
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
c                     omega*omega*r*r
c              v =   _________________
c                           2
c        if the model potential is an exponential the form is:
c
c        here the potential is
c             v = - a*exp(-b*r)
c
c        if the model potential is an constant the form is:
c        
c             v = - constant
c
c
c     everything can be complex
c
      if( dollar('$v',card,cpass,inp) ) then
          pot=chrkey(card,'potential-type','none',' ')
          vr=fpkey(card,'real-v',zero,' ')
          vi=fpkey(card,'imag-v',zero,' ')
          vfac=dcmplx(vr,vi)
          omegar=fpkey(card,'real-omega',zero,' ')
          omegai=fpkey(card,'imag-omega',zero,' ')
          omega=dcmplx(omegar,omegai)
          omega=omega*two*pi
          if(pot.eq.'harmonic-oscillator') then
             scale=one/(omega)
          endif   
          ar=fpkey(card,'real-a',zero,' ')
          ai=fpkey(card,'imag-a',zero,' ')
          a=dcmplx(ar,ai)
          sr=fpkey(card,'real-s',zero,' ')
          si=fpkey(card,'imag-s',zero,' ')
          s=dcmplx(sr,si)
          zr=fpkey(card,'real-z',one,' ')
          zi=fpkey(card,'imag-z',zero,' ')
          z=dcmplx(zr,zi)
      endif          
      call czero(v,n)
      if(pot.eq.'harmonic-oscillator') then
         write(iout,1) omega
         do 10 i=1,n
            v(i)=half*omega*omega**eigc(i)*eigc(i)
 10      continue
      elseif(pot.eq.'exponential') then
         write(iout,2) a, s 
         do 20 i=1,n
            v(i)=a*exp(-s*abs(eigc(i))) 
 20      continue   
      elseif(pot.eq.'well') then
         write(iout,3) vfac 
         do 30 i=1,n
            v(i) = vfac
 30      continue
      elseif(pot.eq.'coulomb') then
         write(iout,4) z 
         do 40 i=1,n
            v(i) = -z/eigc(i)
 40      continue             
      endif
      return
 1    format(/,5x,'frequency:',/,5x,
     1            'omega = ('e15.8,','e15.8')')
 2    format(/,5x,'exponential constants:',/,5x,
     2            'a = ('e15.8,',',e15.8')',/,5x,
     3            's = ('e15.8,',',e15.8')')
 3    format(/,5x,'well constant:',/,5x,
     1            'v = ('e15.8,',',e15.8')')
 4    format(/,5x,'charge:',/,5x,
     1            'z = ('e15.8,',',e15.8')')
      end       

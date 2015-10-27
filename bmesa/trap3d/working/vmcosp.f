*deck vmcosp.f
c***begin prologue     vmcosp
c***date written       970797   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            interaction potential fbr representation
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vmcosp
      subroutine vmcosp(p,q,wt,v,w,t,hbar,mass,scale,n,npt,prnt,ops)
      implicit integer (a-z)
      real*8 p, q, wt, v, w, t
      real*8 hbar, mass, scale
      real*8 pi, zero, half, one, two, three, ten
      real*8 timau
      real*8 aij, sij, omegij, vij
      real*8 facij, fpkey
      character*24 pot1, chrkey
      character*80 title
      character*800 card
      character*(*) ops
      logical prnt, posinp, toau, useau, logkey
      dimension p(npt,0:n-1), q(npt), wt(npt), v(n,n)
      dimension w(npt), t(npt,n)
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
c        if the model potential is an constant the form is:
c        
c             v = - constant
c
c
      if( posinp('$v0',title) ) then
          call cardin(card)      
          pot1=chrkey(card,'unperturbed-potential-type','none',' ')
          scale=one
          vij=fpkey(card,'v11',zero,' ')
          omegij=fpkey(card,'omega-11',zero,' ')
          omegij=omegij*two*pi
          if(pot1.eq.'harmonic-oscillator') then
             scale=one/(omegij*hbar)
          endif   
          aij=fpkey(card,'a11',zero,' ')
          sij=fpkey(card,'s11',zero,' ')
      endif
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      if(toau) then
         omegij=omegij*timau
      endif
      if(useau) then
         omegij=one
      endif
      call rzero(v,n*n)
      if(pot1.eq.'harmonic-oscillator') then
         facij=mass*omegij*omegij*half
         write(iout,1) omegij
         do 10 i=1,npt
            w(i)=facij*q(i)*q(i)*wt(i)
 10      continue
      elseif(pot1.eq.'exponential') then
         write(iout,2) aij, sij
         do 20 i=1,npt
            w(i)=aij*wt(i)*exp(-sij*abs(q(i))) 
 20      continue   
      elseif(pot1.eq.'well') then
         write(iout,3) vij
         do 30 i=1,npt
            w(i)=vij*wt(i)
 30      continue   
      endif
      do 40 i=1,n
         do 50 j=1,npt
            t(j,i)=w(j)*p(j,i-1)
 50      continue
 40   continue   
      call ebtc(v,p,t,n,npt,n)
      return
 1    format(/,5x,'diagonal frequency:',/,5x,
     1               '                  omegaii = ',e15.8)
 2    format(/,5x,'diagonal exponential constant:',/,5x,
     1               '                             aii = ',e15.8,/,5x,
     2               '                             sii = ',e15.8)
 3    format(/,5x,'diagonal well constant:',/,5x,
     1               '                      vii = ',e15.8)
      end       

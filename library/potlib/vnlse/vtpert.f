*deck vtpert.f
c***begin prologue     vtpert
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
c***end prologue       vtpert 
      subroutine vtpert(v,vt,q1,q2,q3,qt,p1,p2,p3,psi,n,nd,nt,nc,
     1                  dim,coord,tim,prn)
      implicit integer (a-z)
      real*8 psi, v, vt, q1, q2, q3, qt
      real*8 omega, scale, width, shift, gamma
      real*8 p1, p2, p3
      real*8 fpkey, scr, mass, sctlen
      real*8 pi, hbar, massau, lenau, timau, pmass, massn2p
      character*80 typ
      character*1600 card
      character*80 cpass, units, chrkey
      character*2 atom
      character* 8 itoc, ii, jj, it
      logical dollar, prn
      dimension nd(3), typ(2)      
      dimension psi(n,nt,nc,2), v(n,nt,nc,nc), vt(*)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3)), qt(nt)
      dimension p1(nd(1),nd(1)), p2(nd(2),nd(2)), p3(nd(3),nd(3))
      dimension atom(2), natom(2), mass(2)
      common/io/inp, iout
      data pi  / 3.141592653589793238462643d+00 /
c     hbar in joule-sec
      data hbar/1.054571596d-34/                                  
      data massau, lenau, timau, pmass / 9.10938188d-31, 
     1                                   5.291772083d-11,
     2                                   2.418884326d-17, 
     3                                   1.67262158d-27 /
      data massn2p / 1.00137841887d0 / 
      pointer (pscr,scr(1))
c
c     deal with pure time potentials first
c 
      call rzero(v,nc*nc*n*nt)
      do 10 i=1,nc
         do 20 j=1,i
            call vcprow(v,vt,i,j,n,nt,nc)
 20      continue
 10   continue   
c	    
      it='t'//itoc(tim)
      len=length(it)
      call getmem(wptoin(nt),pscr,ngot,'scr',0)
      units='atomic-units'
      do 30 i=1,nc
         ii=itoc(i)
         leni=length(ii)
         do 40 j=1,i
            jj=itoc(j)
            lenj=length(jj)
            if(dollar('$vtnlse('//ii(1:leni)//','//jj(1:lenj)//','
     1                         //it(1:len)//')',card,cpass,inp) ) then
               units=chrkey(card,'units',units,' ')
               typ(1)=chrkey(card,'v-space-time','none',' ')
               typ(2)=chrkey(card,'v-non-linear','none',' ')
               write(iout,1) i,j, (typ(k),k=1,2)
               write(iout,2) units
               scale=fpkey(card,'electric-field-strength',1.d0,' ')
               omega=fpkey(card,'electric-field-frequency',10.d0,' ')
               if(typ(1).eq.'cosine-gaussian-pulse'.or.
     1            typ(1).eq.'sine-gaussian-pulse') then
                  width=fpkey(card,'width-of-gaussian-pulse',0.d0,' ')
                  shift=fpkey(card,'shift-of-gaussian-pulse',0.d0,' ')
                  scale=fpkey(card,'scale-of-gaussian-pulse',1.d0,' ')
               endif
               omega=2.d0*pi*omega   
               if(units.eq.'atomic-units') then
                  omega=omega*timau
               endif                      
               if(typ(1).eq.'cosine') then
                  write(iout,3) omega, scale
	          call trgfac(scr,scale,omega,qt,nt,'cosine')
                  call vcprow(v,scr,i,j,n,nt,nc)                 
               elseif(typ(1).eq.'sine') then
                  write(iout,4) omega, scale
	          call trgfac(scr,scale,omega,qt,nt,'sine')
                  call vcprow(v,scr,i,j,n,nt,nc)                 
               elseif(typ(1).eq.'t') then
                  write(iout,5) 
                  call copy(qt,scr,nt)
                  call vcprow(v,scr,i,j,n,nt,nc)                 
               elseif(typ(1).eq.'cosine-dipole-field') then
                  write(iout,6) omega, scale
	          call trgfac(scr,scale,omega,qt,nt,'cosine')
                  call vecprd(v,scr,q1,i,j,n,nt,nc)
               elseif(typ(1).eq.'sine-dipole-field') then
                  write(iout,7) omega, scale
	          call trgfac(scr,scale,omega,qt,nt,'sine')
                  call vecprd(v,scr,q1,i,j,n,nt,nc)
               elseif(typ(1).eq.'cosine-gaussian-pulse') then
                  write(iout,8) scale, width, shift, omega
	          call trgfac(scr,scale,omega,qt,nt,'cosine')
                  call gfac(scr,qt,width,shift,nt)
                  call vfac(v,scr,q1,q2,q3,i,j,n,nd,nt,nc,dim)
               elseif(typ(1).eq.'sine-gaussian-pulse') then
                  write(iout,9) scale, width, shift, omega
	          call trgfac(scr,scale,omega,qt,nt,'sine')
                  call gfac(scr,qt,width,shift,nt)
                  call vfac(v,scr,q1,q2,q3,i,j,n,nd,nt,nc,dim)
               endif
               if(typ(2).ne.'none') then
                  if(tim.eq.1) then
                     call rzero(psi,n*nt*2*nc)
                     return
                  endif
                  atom(1)=chrkey(card,'atom-'//ii(1:leni),'cs',' ')
                  atom(2)=chrkey(card,'atom-'//jj(1:lenj),'cs',' ')
                  sctlen=fpkey(card,'scattering-length',52.d0,' ')
                  natom(1)=intkey(card,'number-of-atoms-'//ii(1:leni),
     1                            0,' ')
                  natom(2)=intkey(card,'number-of-atoms-'//jj(1:lenj),
     1                            0,' ')
                  if(natom(1).ne.0.or.natom(2).ne.0) then
                     write(iout,11) atom(1), atom(2), natom(1), 
     1                                                natom(2),sctlen
                     if(atom(1).eq.'cs') then
                        mass(1)=2.2d-25
                     elseif(atom(1).eq.'na') then
                        mass(1)=3.8176d-26
                     endif
                     if(atom(2).eq.'cs') then
                        mass(2)=2.2d-25
                     elseif(atom(2).eq.'na') then
                        mass(2)=3.8176d-26
                     endif
                     if(units.eq.'atomic-units') then
                        hbar=1.d0
                        mass(1)=mass(1)/massau
                        mass(2)=mass(2)/massau
                     endif
                     gamma = natom(1)*natom(2)*4.d0*pi*hbar*hbar *
     1                       sctlen/(mass(1)*mass(2))
                     call vnlse(v,psi,gamma,p1,p2,p3,i,j,n,nd,nt,nc,dim)
                  endif
               endif
            endif
 40      continue
 30   continue   
      if(prn) then
	 do 50 ic=1,nc
	    do 60 jc=1,ic
               cpass='time-dependent perturbation ch(i) = '//itoc(ic)
     1                //' ch(j) = '//itoc(jc)
	       call prntrm(cpass,v(1,1,ic,jc),n,nt,n,nt,iout)
 60         continue
 50      continue
      endif   
      call getmem(-ngot,pscr,idum,'scr',idum)
      return
 1    format(/,1x,'channel interaction for ',/,1x,
     1            'channel i = ',i3,1x,'channel j = ',i3,/,1x,
     2            'space and time    = ',a32,/,1x,
     3            'non-linear        = ',a32)
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



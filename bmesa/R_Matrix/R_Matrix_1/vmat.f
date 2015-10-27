deck vmat.f
c***begin prologue     vmat
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential matrix elements of DVR/FEM basis.
c***                   
c***references         
c
c***routines called    
c***end prologue       vmat
      subroutine vmat(key,pt,v,dscale,typ,n,zeroit,prn)
      implicit integer (a-z)
      real*8 pt, v, dscale
      real*8 hbar, massau, lenau, timau, mass, pmass, lenscl, escal
      real*8 massn2p, pi, zero, half, one, two, three, ten
      real*8 charg, depth, len, amp, expnt, shift, omega
      real*8 awell, fac, fpkey, scale
      character*(*) key, typ
      character*800 card
      character*80 chrkey, atom, units
      logical dollar, prn, logkey, axis
      logical toau, useau, totrap, zeroit
      dimension pt(n), v(n), amp(2), expnt(2)
      common/io/inp, iout
      data zero, half, one, two, three, ten / 0.d0, .5d0, 1.d0, 2.d0, 
     1                                        3.d0, 10.d0 /
      data pi  / 3.141592653589793238462643d+00 /
c     hbar in joule-sec
      data hbar/1.054571596d-34/                                  
      data massau, lenau, timau, pmass / 9.10938188d-31, 
     1                                   5.291772083d-11,
     2                                   2.418884326d-17, 
     3                                   1.67262158d-27 /
      data massn2p / 1.00137841887d0 /   
      if(zeroit) then
         call rzero(v,n)
      endif            
      if( dollar(key,card,atom,inp) ) then
         units=chrkey(card,'units','atomic-units',' ')
         typ=chrkey(card,'potential','none',' ')
         write(iout,1) typ
         scale=one
         charg=fpkey(card,'charge',-one,' ')
         angmom=intkey(card,'angular-momentum',0,' ')
         depth=fpkey(card,'well-depth',zero,' ')
         len=fpkey(card,'well-size',pt(n),' ')
         number=1
         if(typ.eq.'expres') then
            number=2
         endif
         call fparr(card,'amplitude',amp,number,' ')
         call fparr(card,'exponent',expnt,number,' ')
         shift=fpkey(card,'exponent-shift',zero,' ')
         nwell=intkey(card,'n-well',ten,' ')
         awell=fpkey(card,'a-well',14.d0,' ')
         if(typ.eq.'harmonic-oscillator') then
c
c           take the Cesium atom as model.
c        
            atom=chrkey(card,'atom','generic',' ')
            if(atom.eq.'cs') then
               mass=fpkey(card,'mass',2.2d-25,' ')
               omega=fpkey(card,'omega',ten,' ')
            elseif(atom.eq.'na') then
               mass=3.8176d-26
               omega=13.846d0
            elseif(atom.eq.'na-3d') then
               mass=3.8176d-26
               omega=177.0d0
               axis=logkey(card,'short-axis',.false.,' ')
               if(axis) then
                  omega=sqrt(2.d0)*omega
               endif
            elseif(atom.eq.'generic') then
               mass=fpkey(card,'mass',massau,' ')
               omega=fpkey(card,'omega',1.d0/(two*pi*timau),' ')
            endif
            omega=omega*two*pi
            fac=mass*omega*omega*half
            if(units.eq.'atomic-units') then
               write(iout,*) '   converting to atomic units'
               omega=omega*timau
               mass=mass/massau
               fac=mass*omega*omega*half
               hbar=one
            endif
            write(iout,2) mass, omega
         else
            if(units.eq.'atomic-units') then
               hbar=one
               mass=one
            endif
         endif
      endif
      dscale = - half*hbar*hbar/mass
c
c     calculate the potential
c
      if(typ.eq.'none') then
         call none(v,pt,n,prn)
      elseif(typ.eq.'well') then
         call vwell(v,pt,len,depth,n,prn)
      elseif(typ.eq.'exponential') then
         call vexp(v,pt,amp,expnt,n,prn)
      elseif(typ.eq.'coulomb') then
         call vcoul(v,pt,charg,n,prn)
      elseif(typ.eq.'inverse-r4') then
         call vir4(v,pt,n,prn)
      elseif(typ.eq.'rounded-well') then
         call vrwell(v,pt,awell,nwell,n,prn)
      elseif(typ.eq.'harmonic-oscillator') then
         call vhmo(v,pt,fac,n,prn)
      elseif(typ.eq.'anharmonic-oscillator') then
         call vanhmo(v,pt,n,prn)
      elseif(typ.eq.'expres') then
         call vres(v,pt,amp,expnt,shift,n,prn)
      else            
         call lnkerr('error in potential')
      endif
      if(angmom.ne.0) then
         call addang(v,pt,angmom,key,n,prn)
      endif
      return
 1    format(/,1x,'potential type = ',a32)
 2    format(/,1x,'oscillator mass      = ',e15.8,
     1       /,1x,'oscillator-frequency = ',e15.8)
 3    format(/,1x,'length scale = ',e15.8,1x,'energy scale = ',e15.8)
      end       




!deck vone.f
!***begin prologue     vone
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, barry (nsf)
!***source             
!***purpose            potential matrix elements of DVR/FEM basis.
!***                   
!***references         
!
!***routines called    
!***end prologue       vone
      subroutine vone(coord,pt,v,dscale,typ,nglobal,zeroit,prnt)
      USE input_output
      USE grid_global
      implicit none
      integer                     :: nglobal
      real*8, dimension(nglobal)  :: pt, v
      real*8                      :: dscale, escal, charg, depth
      real*8                      :: len, shift, omega
      real*8                      :: awell, fac, fpkey, scale
      character*(*)               :: coord, typ
      character*80                :: chrkey, atom
      logical                     :: dollar, logkey, axis
      logical                     :: toau, useau, totrap, zeroit, prnt
      real*8, dimension(3)        :: amp, expnt
      integer                     :: number, nwell, n_p, intkey
      integer                     :: lenth, lenkey
      if(zeroit) then
         v=0.d0
      endif            
      lenkey=lenth(coord)
      if( dollar('$v0('//coord(1:lenkey)//')',card,atom,inp) ) &
          then
         units=chrkey(card,'units','atomic-units',' ')
         typ=chrkey(card,'potential','none',' ')
         write(iout,1) typ
         scale=one
         angmom=intkey(card,'angular-momentum',0,' ')
         if(typ.eq.'harmonic-oscillator') then
!
!           take the Cesium atom as model.
!        
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
      else
         call lnkerr('quit')
      end if
      dscale = - half*hbar*hbar/mass
!
!     calculate the potential
!
      if(typ.eq.'none') then
         call none(v,pt,nglobal,prnt)
      elseif(typ.eq.'well') then
         depth=fpkey(card,'well-depth',zero,' ')
         len=fpkey(card,'well-size',pt(nglobal),' ')
         call vwell(v,pt,len,depth,nglobal,prnt)
      elseif(typ.eq.'exponential') then
         amp(1)=fpkey(card,'amplitude',amp,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         call vexp(v,pt,amp,expnt,nglobal,prnt)
      elseif(typ.eq.'power-exponential') then
         amp(1)=fpkey(card,'amplitude',amp,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         n_p=intkey(card,'power',0,' ')
         call v_pow_exp(v,pt,amp,expnt,n_p,nglobal,prnt)
      elseif(typ.eq.'sum-exponential') then
         call fparr(card,'amplitudes',amp,2,' ')
         call fparr(card,'exponents',expnt,2,' ')
         call vexp_sum(v,pt,amp,expnt,nglobal,prnt)
      elseif(typ.eq.'coulomb') then
         charg=fpkey(card,'charge',-one,' ')
         call vcoul(v,pt,charg,nglobal,prnt)
      elseif(typ.eq.'inverse-r4') then
         call vir4(v,pt,nglobal,prnt)
      elseif(typ.eq.'rounded-well') then
         nwell=intkey(card,'n-well',ten,' ')
         awell=fpkey(card,'a-well',14.d0,' ')
         call vrwell(v,pt,awell,nwell,nglobal,prnt)
      elseif(typ.eq.'harmonic-oscillator') then
         call vhmo(v,pt,fac,nglobal,prnt)
      elseif(typ.eq.'anharmonic-oscillator') then
         call vanhmo(v,pt,nglobal,prnt)
      elseif(typ.eq.'expres') then
         call fparr(card,'amplitude',amp,2,' ')
         call fparr(card,'exponent',expnt,2,' ')
         shift=fpkey(card,'exponent-shift',zero,' ')
         call vres(v,pt,amp,expnt,shift,nglobal,prnt)
      elseif(typ.eq.'eberly') then
         charg=fpkey(card,'charge',-one,' ')     
         n_p=intkey(card,'power',0,' ')
         amp(1)=fpkey(card,'amplitude',one,' ')
         shift=fpkey(card,'shift',one,' ')
         call v_eberlonium(v,pt,charg,amp(1),shift,n_p, &
                           nglobal,prnt)
      elseif(typ.eq.'morse') then
         amp(1)=fpkey(card,'amplitude',one,' ')     
         amp(2)=fpkey(card,'constant',amp(1),' ')     
         amp(3)=fpkey(card,'static-field',0.d0,' ')
         expnt(1)=fpkey(card,'first-exponent',2.d0/sqrt(8.d0*amp(1)),' ')
         expnt(2)=fpkey(card,'second-exponent',.5d0*expnt(1),' ')
         call vmorse(v,pt,amp(1),amp(2),expnt(1),expnt(2), &
                     amp(3),nglobal,prnt)
      else            
         call lnkerr('error in potential')
      endif
      return
 1    format(/,1x,'potential type = ',a32)
 2    format(/,1x,'oscillator mass      = ',e15.8, &
             /,1x,'oscillator-frequency = ',e15.8)
 3    format(/,1x,'length scale = ',e15.8,1x,'energy scale = ',e15.8)
      end       




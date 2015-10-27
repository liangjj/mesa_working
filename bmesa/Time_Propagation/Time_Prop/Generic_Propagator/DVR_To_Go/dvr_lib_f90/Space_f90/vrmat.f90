! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Generalized Potential Energy Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck vrmat.f
!***begin prologue     vrmat
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate potential energy matrix elements
!***references
!***routines called
!***end prologue       vrmat
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE vrmat(v,pt,coord,dscale,nr,region)
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL*8, DIMENSION(nr)                  :: v, pt
  CHARACTER (LEN=*)                      :: coord
  CHARACTER (LEN=3)                      :: itoc
  CHARACTER (LEN=16)                     :: datkey
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: region
  INTEGER                                :: i
  REAL*8                                 :: dscale, escal, charg, depth
  REAL*8                                 :: len, shift, omega
  REAL*8                                 :: awell, fac, fpkey, scale
  CHARACTER (LEN=80)                     :: typ, chrkey, atom
  LOGICAL                                :: dollar, logkey, axis
  LOGICAL                                :: toau, useau, totrap, prnt
  REAL*8, DIMENSION(2)                   :: amp, expnt
  REAL*8                                 :: e_c
  INTEGER                                :: number, nwell, n_p, intkey, n_scale
  INTEGER                                :: lendat, lencord
  INTEGER                                :: lenth
  v=0.d0
  if(reuse) then
     datkey = '$v_reg_'//itoc(1)
  else
     datkey = '$v_reg_'//itoc(region)
  END IF
  lendat=lenth(datkey)
  lencord=lenth(coord)
  if( dollar( datkey(1:lendat)//'('// coord(1:lencord)//')', &
              card,atom,inp) )  then
      typ=chrkey(card,'potential','none',' ')
      scale=one
      angmom=intkey(card,'angular-momentum',angmom,' ')
      m_val=angmom
      write(iout,1) typ, angmom
!
!     calculate the potential
!
      if(typ == 'none') then
         call none(v,pt,nr,prnt)
      elseif(typ == 'well') then
         depth=fpkey(card,'well-depth',zero,' ')
         call vwell(v,depth,nr,prnt)
      elseif(typ == 'exponential') then
         amp(1)=fpkey(card,'amplitude',-1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         call vexp(v,pt,amp,expnt,nr,prnt)
      elseif(typ == 'yukawa') then
         amp(1)=fpkey(card,'amplitude',-1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         call vyukawa(v,pt,amp,expnt,nr,prnt)
      elseif(typ == 'power-exponential') then
         amp(1)=fpkey(card,'amplitude',1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         n_p=intkey(card,'power',0,' ')
         call v_pow_exp(v,pt,amp,expnt,n_p,nr,prnt)
      elseif(typ == 'sum-exponential') then
         call fparr(card,'amplitudes',amp,2,' ')
         call fparr(card,'exponents',expnt,2,' ')
         call vexp_sum(v,pt,amp,expnt,nr,prnt)
      elseif(typ == 'coulomb') then
         charg=fpkey(card,'charge',-one,' ')
         call vcoul(v,pt,charg,nr,prnt)
      elseif(typ == 'eberlonium') then
         charg=fpkey(card,'charge',-one,' ')         
         n_p=intkey(card,'power',0,' ')
         amp(1)=fpkey(card,'a',1.d0,' ')
         amp(2)=fpkey(card,'b',1.d0,' ')
         call v_eberlonium(v,pt,charg,amp(1),amp(2),n_p,nr,prnt)
      elseif(typ == 'inverse-r4') then
         call vir4(v,pt,nr,prnt)
      elseif(typ == 'rounded-well') then
         nwell=intkey(card,'n-well',ten,' ')
         awell=fpkey(card,'a-well',14.d0,' ')
         call vrwell(v,pt,awell,nwell,nr,prnt)
      elseif(typ == 'harmonic-oscillator') then
!
!        take the Cesium atom as model.
!         
         atom=chrkey(card,'atom','generic',' ')
         if(atom == 'cs') then
            mass=fpkey(card,'mass',2.2d-25,' ')
            omega=fpkey(card,'omega',ten,' ')
         elseif(atom == 'na') then
            mass=3.8176d-26
            omega=13.846d0
         elseif(atom == 'na-3d') then
            mass=3.8176d-26
            omega=177.0d0
            axis=logkey(card,'short-axis',.false.,' ')
            if(axis) then
              omega=sqrt(2.d0)*omega
            endif
         elseif(atom == 'generic') then
           mass=fpkey(card,'mass',massau,' ')
           omega=fpkey(card,'omega',1.d0/(two*pi*timau),' ')
         endif
         omega=omega*two*pi
         fac=mass*omega*omega*half
         if(units == 'atomic-units') then
            write(iout,*) '   converting to atomic units'
            omega=omega*timau
            mass=mass/massau
            fac=mass*omega*omega*half
            hbar=one
         endif
         write(iout,2) mass, omega
         call vhmo(v,pt,fac,nr,prnt)
      elseif(typ == 'anharmonic-oscillator') then
         call vanhmo(v,pt,nr,prnt)

      elseif(typ == 'expres') then
         call fparr(card,'amplitude',amp,2,' ')
         call fparr(card,'exponent',expnt,2,' ')
         shift=fpkey(card,'exponent-shift',zero,' ')
         call vres(v,pt,amp,expnt,shift,nr,prnt)
      elseif(typ == 'periodic') then
         n_scale=intkey(card,'n_i',10,' ')         
         e_c=fpkey(card,'e_c',.001d0,' ')         
         awell=n_scale/e_c
         call vperiod(v,pt,awell,nr,prnt)
      else
         call lnkerr('error in potential')
      endif
  END IF
  dscale = -.5d0*hbar*hbar/mass
  if(units == 'atomic-units') then
     dscale = -.5d0
  END IF     
  if(angmom /= 0) then
     call addang(v,pt,dscale,angmom,coord,nr,.false.)
  end if
  if(prn(3)) then
     title='potential matrix elements for region = '//itoc(region)
     call prntrm(title,v,nr,1,nr,1,iout)
  endif
 1    format(/,1x,'potential type = ',a32, &
             /,1x,'angular momentum = ',i3)
 2    format(/,1x,'oscillator mass      = ',e15.8, &
             /,1x,'oscillator-frequency = ',e15.8)
 3    format(/,1x,'length scale = ',e15.8,1x,'energy scale = ',e15.8)
END SUBROUTINE vrmat

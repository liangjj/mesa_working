!deck v_mat.f
!***begin prologue     v_mat
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate potential energy matrix elements
!***references
!***routines called
!***end prologue       v_mat
  SUBROUTINE v_mat(v,pt,case,coord,n)
  USE input_output
  IMPLICIT NONE
  REAL*8, DIMENSION(n)                  :: v
  REAL*8, DIMENSION(n)                  :: pt
  INTEGER                               :: case
  INTEGER                               :: n
  CHARACTER (LEN=*)                     :: coord
  CHARACTER (LEN=3)                     :: itoc
  LOGICAL                               :: dollar
  LOGICAL                               :: logkey
  CHARACTER (LEN=80)                    :: chrkey
  CHARACTER (LEN=16)                    :: datkey
  CHARACTER (LEN=80)                    :: title
  INTEGER                               :: i
  REAL*8                                :: zero =0.d0
  REAL*8                                :: one = 1.d0
  REAL*8                                :: two = 1.d0
  REAL*8                                :: ten = 10.d0
  REAL*8                                :: half = .5d0
  REAL*8                                :: pi=3.141592653589793238462643D+00
  REAL*8                                :: mass
  REAL*8                                :: hbar
  REAL*8                                :: dscale
  REAL*8                                :: escal
  REAL*8                                :: charge
  REAL*8                                :: depth
  REAL*8                                :: len
  REAL*8                                :: shift
  REAL*8                                :: omega
  REAL*8                                :: awell
  REAL*8                                :: fac
  REAL*8                                :: fpkey
  REAL*8                                :: scale
  CHARACTER (LEN=400)                   :: card
  CHARACTER (LEN=80)                    :: typ
  CHARACTER (LEN=80)                    :: atom
  CHARACTER (LEN=8)                     :: chrcase
  LOGICAL                               :: axis
  LOGICAL                               :: toau
  LOGICAL                               :: useau
  LOGICAL                               :: totrap
  LOGICAL                               :: prnt
  REAL*8, DIMENSION(2)                  :: amp
  REAL*8, DIMENSION(2)                  :: expnt
  REAL*8                                :: e_c
  INTEGER                               :: number
  INTEGER                               :: nwell
  INTEGER                               :: n_p
  INTEGER                               :: intkey
  INTEGER                               :: n_scale
  INTEGER                               :: lencord
  INTEGER                               :: lencase
  INTEGER                               :: lenth
  v(:)=0.d0
  chrcase=itoc(case)
  lencase=lenth(chrcase)
  lencord=lenth(coord)
  datkey = '$v_'//chrcase(1:lencase)//'('// coord(1:lencord)//')'
  lencord=lenth(datkey)
  if( dollar( datkey,card,atom,inp) )  then
      typ=chrkey(card,'potential','none',' ')
      scale=one
      write(iout,1) typ
!
!     calculate the potential
!
      if(typ == 'none') then
         call none(v,pt,n,prnt)
      elseif(typ == 'well') then
         depth=fpkey(card,'well_depth',zero,' ')
         call vwell(v,depth,n,prnt)
      elseif(typ == 'exponential') then
         amp(1)=fpkey(card,'amplitude',-1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         call vexp(v,pt,amp,expnt,n,prnt)
      elseif(typ == 'yukawa') then
         amp(1)=fpkey(card,'amplitude',-1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         call vyukawa(v,pt,amp,expnt,n,prnt)
      elseif(typ == 'power_exponential') then
         amp(1)=fpkey(card,'amplitude',1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         n_p=intkey(card,'power',0,' ')
         call v_pow_exp(v,pt,amp,expnt,n_p,n,prnt)
      elseif(typ == 'sum_exponential') then
         call fparr(card,'amplitudes',amp,2,' ')
         call fparr(card,'exponents',expnt,2,' ')
         call vexp_sum(v,pt,amp,expnt,n,prnt)
      elseif(typ == 'coulomb') then
         charge=fpkey(card,'charge',-one,' ')
         call vcoul(v,pt,charge,n,prnt)
      elseif(typ == 'eberlonium') then
         charge=fpkey(card,'charge',-one,' ')         
         n_p=intkey(card,'power',0,' ')
         amp(1)=fpkey(card,'a',1.d0,' ')
         amp(2)=fpkey(card,'b',1.d0,' ')
         call v_eberlonium(v,pt,charge,amp(1),amp(2),n_p,n,prnt)
      elseif(typ == 'inverse_r4') then
         call vir4(v,pt,n,prnt)
      elseif(typ == 'rounded_well') then
         nwell=intkey(card,'n_well',ten,' ')
         awell=fpkey(card,'a_well',14.d0,' ')
         call vrwell(v,pt,awell,nwell,n,prnt)
      elseif(typ == 'harmonic_oscillator') then
         mass=one
         omega=fpkey(card,'omega',1.d0,' ')
         fac=mass*omega*omega*half
         hbar=one
         write(iout,2) mass, omega
         call vhmo(v,pt,fac,n,prnt)
      elseif(typ == 'anharmonic_oscillator') then
         call vanhmo(v,pt,n,prnt)
      elseif(typ == 'expres') then
         call fparr(card,'amplitude',amp,2,' ')
         call fparr(card,'exponent',expnt,2,' ')
         shift=fpkey(card,'exponent_shift',zero,' ')
         call vres(v,pt,amp,expnt,shift,n,prnt)
      elseif(typ == 'periodic') then
         n_scale=intkey(card,'n_i',10,' ')         
         e_c=fpkey(card,'e_c',.001d0,' ')         
         awell=n_scale/e_c
         call vperiod(v,pt,awell,n,prnt)
      else
         call lnkerr('error in potential')
      endif
  END IF
 1    format(/,1x,'potential type = ',a32)
 2    format(/,1x,'oscillator mass      = ',e15.8, &
             /,1x,'oscillator-frequency = ',e15.8)
 3    format(/,1x,'length scale = ',e15.8,1x,'energy scale = ',e15.8)
END SUBROUTINE v_mat

!***********************************************************************                
                           MODULE Radial_Quadrature
!deck Rangular_Quadrature
!***begin prologue     Radial_Quadrature
!***date written       20140706  (yyyymmdd)                                                
!***revision date                (yyyymmdd)                                                
!***keywords           Radial_Quadrature
!***author             schneider, b. i.(nist)                                            
!***source                                                                        
!***purpose            Compute the points and weights for a radial quadrature
!***                   based on either Gauss or Nwton-Cotes rules.
!description           The radial region is broken down into shells.                     
!***                   In each shell a different size quadrature may be used.
!***                   In most cases the quadrature should be chosen so that the first point of region
!***                   (i+1) and the last point of region i coincide.  This will allow
!***                   a finite element DVR to be applied in a straightforward manner.
!***                   
!***                   
!***                   
!***                   
!***references                                                                          
!***routines called                                                                     
!***end prologue                                       
                           USE input_output
                           USE accuracy
                           USE Data
                           USE Grid_Defined_Types
                           IMPLICIT NONE
!**********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Radial_Grid                       
                       MODULE PROCEDURE Gauss_Radial_Quadrature,          &
                                        Ang_06,                            &
                                        Ang_14,                            &
                                        Ang_26,                            &  
                                        Ang_38,                            &  
                                        Ang_50,                            &  
                                        Ang_74,                            &  
                                        Ang_86,                            &  
                                        Ang_110,                           &  
                                        Ang_146,                           &  
                                        Ang_170,                           &  
                                        Ang_194,                           &   
                                        Ang_230,                           &  
                                        Ang_266,                           &  
                                        Ang_302,                           &  
                                        Ang_350,                           &  
                                        Ang_434,                           &  
                                        Ang_590,                           &  
                                        Ang_770,                           &  
                                        Ang_974,                           &  
                                        Ang_1202,                          &  
                                        Ang_1454,                          &  
                                        Ang_1730,                          &  
                                        Ang_2030,                          &  
                                        Ang_2354,                          &  
                                        Ang_2702,                          &  
                                        Ang_3074,                          &  
                                        Ang_3470,                          &  
                                        Ang_3890,                          &  
                                        Ang_4334,                          &  
                                        Ang_4802,                          &
                                        Ang_5294,                          &  
                                        Ang_5810  
                            END INTERFACE Radial_Grid                       
!
!
!
!
!************************************************************************                
!!***********************************************************************               
                           Contains
!***********************************************************************                
!***********************************************************************                
!deck Radial_Points_Weights.f
!***begin prologue     Radial_Points_Weights
!***date written       9308021   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           shells
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            calculate radial quadrature points
!***                   and weights
!***references
!***routines called    gaussq ( math )
!***end prologue       Radial_Points_Weights

SUBROUTINE Radial_Points_Weights(atom,ang_grd,rad_grd,TYPE,icen,cen,r,nr,eta,rpt,thpt,phpt,wtr,wtsum,  &
    wtth,wtph,angleb,work,wtleb,sthet,sphi,cphi,  &
    yuk,scr,grid,wt,spgrid,ns,nrmax,nthet,nphi,  &
    ngrid,nwts,ncen,ncplus,prnt,nleb,nang,nonsep, yukon,nodisk)
  IMPLICIT NONE
  TYPE(ATOMS)                            :: atom
  TYPE(RAD_PT_WT)                        :: rad_grd
CHARACTER (LEN=*), INTENT(IN)            :: TYPE
INTEGER, INTENT(IN OUT)                  :: icen
REAL, INTENT(IN OUT)                     :: cen(3,ncplus)
REAL, INTENT(IN OUT)                     :: r(ns+1)
INTEGER, INTENT(IN)                      :: nr(ns)
REAL, INTENT(IN)                         :: eta(ncplus)
REAL, INTENT(OUT)                        :: rpt(nrmax,ns)
REAL, INTENT(OUT)                        :: thpt(nthet)
REAL, INTENT(OUT)                        :: phpt(nphi)
REAL, INTENT(OUT)                        :: wtr(*)
REAL, INTENT(IN OUT)                     :: wtsum(*)
REAL, INTENT(OUT)                        :: wtth(nthet)
REAL, INTENT(OUT)                        :: wtph(nphi)
REAL, INTENT(IN OUT)                     :: angleb(*)
REAL, INTENT(IN OUT)                     :: work(*)
REAL, INTENT(IN OUT)                     :: wtleb(*)
REAL, INTENT(OUT)                        :: sthet(nthet)
REAL, INTENT(OUT)                        :: sphi(nphi)
REAL, INTENT(OUT)                        :: cphi(nphi)
REAL, INTENT(IN OUT)                     :: yuk(*)
REAL, INTENT(IN OUT)                     :: scr(*)
REAL, INTENT(IN OUT)                     :: grid(3,ngrid)
REAL, INTENT(IN OUT)                     :: wt(nwts)
REAL, INTENT(IN OUT)                     :: spgrid(3,ngrid)
INTEGER, INTENT(IN)                      :: ns
INTEGER, INTENT(IN)                      :: nrmax
INTEGER, INTENT(IN)                      :: nthet
INTEGER, INTENT(IN)                      :: nphi
INTEGER, INTENT(IN)                      :: ngrid
INTEGER, INTENT(IN)                      :: nwts
INTEGER, INTENT(IN)                      :: ncen
INTEGER, INTENT(IN OUT)                  :: ncplus
LOGICAL, INTENT(IN)                      :: prnt
INTEGER, INTENT(IN OUT)                  :: nleb
INTEGER, INTENT(OUT)                     :: nang
LOGICAL, INTENT(IN)                      :: nonsep
LOGICAL, INTENT(IN)                      :: yukon
LOGICAL, INTENT(IN OUT)                  :: nodisk
INTEGER :: thefac, phifac

CHARACTER (LEN=30) :: rtyp, thtyp, phtyp, str
CHARACTER (LEN=3) :: itoc, chra
DIMENSION  delthe(2), delphi(2)
DIMENSION  dummy(2)
!
!         do the radial quadrature points. they are divided into shells.
!

  ngcnt=0
  nwtcnt=0
  aprvol=0.d0
  yukint=0.d0
  ALLOCATE(atom%_yukawa(1:atom%ngrid) 
  CALL rzero(atom%yukawa,atom%ngrid)
  icount=0
  locwtr=nrmax*(i-1)+1
!     get radial quadrature points and weights for this shell
  IF (atom%rtyp == 'legendre') THEN
    ampr=(r(i+1)-r(i))*.5D0
    abpr=(r(i+1)+r(i))*.5D0
    CALL gaussq(rtyp,nr(i),0.d+00,0.d+00,0,dummy,scr, rpt(1,i),wtr(locwtr))
    nwtsr=nr(i)
    icount=icount+nr(i)
    loccnt=locwtr
    DO  j=1,nr(i)
      numr=numr+1
      rpt(j,i)=ampr*rpt(j,i)+abpr
      wtr(loccnt)=ampr*wtr(loccnt)
      loccnt=loccnt+1
    END DO
  ELSE IF(rtyp == 'newton-cotes') THEN
    locwtr=nrmax*(nrmax-1)*(i-1)+1
    CALL necote(r(i),r(i+1),rpt(1,i),wtr(locwtr),nr(i), .false.)
    nwtsr=(nr(i)-1)*nr(i)
    icount=icount+nr(i)
  ELSE
    CALL lnkerr('error in radial quadrature')
  END IF
  CALL iosys ('write real "radial points '//str//  &
      '" to lamdat without rewinding',nr(i), rpt(1,i),0,' ')
  CALL iosys ('write real "unscaled radial weights '//  &
      str//'" to lamdat without '// 'rewinding',nwtsr,wtr(locwtr),0,' ')
  IF(rtyp == 'newton-cotes') THEN
    CALL sumncw(wtr(locwtr),wtsum,nr(i))
    CALL iosys ('write real "summed radial weights '//  &
        str//'" to lamdat without rewinding', nr(i),wtsum,0,' ')
  END IF
!      the procedure is different if the angular quadrature is separable or
!      non-separable in  (theta,phi)
  CALL mkgr(grid(1,ngcnt+1),spgrid(1,ngcnt+1),rpt(1,i),thpt,  &
      sthet,sphi,cphi,cen(1,icen),nr(i),nthet, nphi,nang,nonsep)
  IF (yukon) THEN
    CALL yukawa(yuk(ngcnt+1),grid(1,ngcnt+1),eta,cen,  &
        nr(i),nthet,nphi,ncen,nang,nonsep)
  END IF
  CALL mkwt(wtr(locwtr),wtr(locwtr),wtth,wtph,wtleb,  &
      wt(nwtcnt+1),work,nr(i),nthet,nphi,rtyp,nang,nonsep)
  IF (yukon) THEN
    CALL mkyunt(rpt(1,i),wt(nwtcnt+1),yuk(ngcnt+1),aprvol,  &
        yukint,work,nr(i),nthet,nphi,rtyp,nang, nonsep)
  END IF
  ngcnt=ngcnt+nr(i)*nang
  nwtcnt=nwtcnt+nwtsr*nang
  CALL scalwt(rpt(1,i),wtr(locwtr),wtr(locwtr),work,nr(i),rtyp)
  CALL iosys ('write real "scaled radial weights '//  &
      str//'" to lamdat without '// 'rewinding',nwtsr,wtr(locwtr),0,' ')
  IF (prnt) THEN
    WRITE(iout,*) '     data on radial points for shell = ',i
    WRITE(iout,*)
    WRITE(iout,*) 'radial points and weights'
    WRITE(iout,2) (rpt(ii,i),ii=1,nr(i))
    WRITE(iout,2) (wtr(ii),ii=locwtr,locwtr+nwtsr-1)
  END IF
END DO
IF (ngcnt /= ngrid) THEN
  CALL lnkerr('error in grid point count')
END IF
IF (nwtcnt /= nwts) THEN
  CALL lnkerr('error in grid weight count')
END IF
CALL iosys ('write integer "total number of radial '//  &
    'points '//str//'" to lamdat',1,icount,0,' ')
IF (.NOT.nodisk) THEN
  CALL iosys ('write real "yukawa potential '//  &
      str//' coordinates" to lamdat',ngcnt, yuk,0,' ')
  CALL iosys ('write real "atomic grid '//  &
      str//'" to lamdat',3*ngcnt,grid,0,' ')
  CALL iosys ('write real "spherical atomic grid '//  &
      str//'" to lamdat',3*ngcnt,spgrid,0,' ')
  CALL iosys ('write real "unscaled atomic weights '//  &
      str//'" to lamdat',nwts,wt,0,' ')
END IF
IF (yukon) THEN
  yukext=0.d0
  DO  nc=1,ncen
    yukext=yukext+4.d+00*pi/eta(nc)**2
  END DO
  exvol=4.d0*pi*r(ns+1)**3/3.d0
  WRITE(iout,3) exvol
  WRITE(iout,4) aprvol
  WRITE(iout,5) yukext
  WRITE(iout,6) yukint
END IF
  END Subroutine Radial_Points_Weights
!*****************************************************************************80
!*****************************************************************************80
  END MODULE Radial_Quadrature
!*****************************************************************************80

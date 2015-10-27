!deck Points_and_Weights
!***begin prologue     Points_and_Weights
!***date written       140601   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Points_and_Weights
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            Compute all information needed for each atomic grid.
!***                   The grids are a tensor product of the radial*theta*phi points
!***                   for each atom EXCEPT when Lebedev quadratures are used.
!***                   Then it is the product of the radial*nleb points,
!***description  
!***             
!***             
!***references
!***routines called    
!***end prologue       Points_and_Weights
!********************************************************************************
!********************************************************************************
                        MODULE Points_and_Weights
  USE Angular_Quadrature
  USE Radial_Quadrature
  USE Grid_Defined_Types
  IMPLICIT NONE
!********************************************************************************
!********************************************************************************
!
                                 Contains
!
!********************************************************************************
!********************************************************************************
!deck shells.f
!***begin prologue     shells
!***date written       9308021   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           shells
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            calculate angular and radial quadrature points
!***                   and weights
!***references
!***routines called    gaussq ( math )
!***end prologue       shells

SUBROUTINE shells(ang_grd,rad_grd,TYPE,icen,cen,r,nr,eta,rpt,thpt,phpt,wtr,wtsum,  &
    wtth,wtph,angleb,work,wtleb,sthet,sphi,cphi,  &
    yuk,scr,grid,wt,spgrid,ns,nrmax,nthet,nphi,  &
    ngrid,nwts,ncen,ncplus,prnt,nleb,nang,nonsep, yukon,nodisk)
  IMPLICIT NONE
  TYPE(ANG_PT_WT)                        :: ang_grd
  TYPE(RAD_PT_WT)                        :: rad_grd
  TYPE(RULE)                             :: leb_rule
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


!  CALL iosys('read character "radial quadrature type '//str//'" from grid',-1,0,0,atom_rad_grd%rtyp)
!  CALL iosys('read real "theta range '// str//'" from grid',2,atom%ang_grd%range,0,' ')
!  CALL iosys('read real "phi range '// str//'" from grid',2,atom%ang_grd%range(3),0,' ')
!  CALL iosys ('read real "radial edges '//str//'" from grid',nedges,atom_rad_grd%r,0,' ')
  CALL iosys ('create real "radial points '//str//'" on grid',atom%rad_grd%ntrad,0,0,' ')
  CALL iosys ('create real "scaled radial weights '// str//'" on grid',atom%rad_grd%ntradwt, 0,0,' ')
  CALL iosys ('create real "unscaled radial weights '//str//'" on grid',atom%rad_grd%ntradwt, 0,0,' ')
  nedges=rad_grd%nshell+1
  IF (rad_grd%rtyp == 'newton-cotes') THEN
      CALL iosys ('create real "summed radial weights '//str//'" on grid',rad_grd%ntradwt,0,0,' ')
  END IF
  IF (.NOT.nodisk) THEN
      CALL iosys ('create real "yukawa potential '//str//' coordinates" on grid',ngrid, 0,0,' ')
      CALL iosys ('create real "atomic grid '//str//' " on grid',3*ngrid,0,0,' ')
      CALL iosys ('create real "spherical atomic grid '//str//'" on grid',3*ngrid,0,0,' ')
      CALL iosys ('create real "unscaled atomic weights '//str//'" on grid',nwts,0,0,' ')
  END IF
  IF (nonsep) THEN
  
!          in a lebedev quadrature the number of theta and phi points are
!          identical since the quadrature is non-separable. the combined
!          angular weight is stored in wtth.
!
      Call Generate_Lebedev_Points_Weights(ang_grd,leb_rule)
!
      ang_grd%sumwt=0.d0
      DO i=1,ang_grd%nang
         ang_grd%w(i)=ang_grd%w(i)*fourpi 
         ang_grd%sumwt=ang_grd%sumwt+ang_grd%w(i)
      END DO
      ang_grd%sumwt=ang_grd%sumwt/fourpi
      write(iout,1) ' sum of lebedev weights = ', ang_grd%sumwt
!
!     Compute the spherical variables based on the Lebedev quadrature
!
      ALLOCATE(ang_grd%cthet(1:ang_grd%nang), ang_grd%sthet(1:ang_grd%nang),   &
               ang_grd%phipt(1:ang_grd%nang), ang_grd%sphi(1:ang_grd%nang),    &
               ang_grd%cphi(1:ang_grd%nang) )
      Call ang(ang_grd)
  ELSE
!      CALL iosys('read character "theta quadrature type '//str//'" from grid',-1,0,0,atom_ang_grd%thtyp)
!      CALL iosys('read character "phi quadrature type '//str//'" from grid',-1,0,0,ang_grd%phtyp)
!
      Call Gauss_Angular_Quadrature (ang_grd,irange=.true.)
      CALL iosys ('write real "theta points '//str//   &
                  '" to grid',ang_grd%nthet,ang_grd%thetpt,0,' ')
      CALL iosys ('write real "theta weights '//str//  &
                  '" to grid',ang_grd%nthet,ang_grd%wtthet,0,' ')
      CALL iosys ('write real "phi points '//str// '" to grid',ang_grd%nphi,ang_grd%phipt,0,' ')
      CALL iosys ('write real "phi weights '//str// '" to grid',ang_grd%nphi,ang_grd%wtphi,0,' ')
      ang_grd%nang=ang_grd%nthet*ang_grd%nphi
  END IF
  IF (prnt) THEN
      WRITE(iout,*)
      WRITE(iout,*) 'cos(theta) points and weights'
      WRITE(iout,2) (ang_grd%thetpt(ii),ii=1,ang_grd%nthet)
      WRITE(iout,2) (ang_grd%wtthet(ii),ii=1,ang_grd%nthet)
      WRITE(iout,*)
      WRITE(iout,*) 'phi points and weights'
      WRITE(iout,2) (ang_grd%phipt(ii),ii=1,ang_grd%nphi)
      WRITE(iout,2) (ang_grd%wtphi(ii),ii=1,ang_grd%nphi)
  END IF

!         do the radial quadrature points. they are divided into shells.


ngcnt=0
nwtcnt=0
aprvol=0.d0
yukint=0.d0
CALL rzero(yuk,ngrid)
icount=0
DO  i=1,ns
  IF (r(i) == 0.d0) THEN
    r(i)=1.d-10
  END IF
  locwtr=nrmax*(i-1)+1
!     get radial quadrature points and weights for this shell
  IF (rtyp == 'legendre') THEN
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
RETURN
1 FORMAT(' multiplicative factors for theta and phi  :',2F10.5)
2 FORMAT(/4(2X,e15.8))
3 FORMAT(/,5X,'exact       value of volume =',e15.8)
4 FORMAT(/,5X,'approximate value of volume =',e15.8)
5 FORMAT(/,5X,'exact     value of yukawa integral = ',e15.8)
6 FORMAT(/,5X,'numerical value of yukawa integral = ',e15.8)
END SUBROUTINE shells
!********************************************************************************
!********************************************************************************
  END MODULE Points_and_Weights
!********************************************************************************
!********************************************************************************

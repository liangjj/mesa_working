!deck @(#)basout.f 1.1 9/7/91
!***begin prologue     basout
!***date written       901214   (yymmdd)
!***revision date               (yymmdd)
!***keywords           mesa data
!***author             schneider, barry (lanl)
!***source             m8000
!***purpose            
!***references
!***routines called    iosys, util and mdutil
!***end prologue       basout

  SUBROUTINE basout
  USE mesa_global  
  IMPLICIT NONE
  REAL *8, DIMENSION(300,5)                   :: eta
  INTEGER, DIMENSION(300)                     :: nstart, nstop, lnew, mnew, nnew, cnt
  CHARACTER (LEN=4)                           :: itoc
  INTEGER                                     :: i, j, jj, ikount, icont, iloc
  INTEGER                                     :: kbf, iatom, itype, mini, maxi
  INTEGER                                     :: imax, nprimi, nconti
  INTEGER                                     :: icount, m, iii, iiirel
  INTEGER                                     :: npr, ncon, nlow, nhi
  ctype(0,0,0)='s'
  ctype(1,0,0)='x'
  ctype(0,1,0)='y'
  ctype(0,0,1)='z'
  ctype(2,0,0)='xx'
  ctype(0,2,0)='yy'
  ctype(0,0,2)='zz'
  ctype(1,1,0)='xy'
  ctype(1,0,1)='xz'
  ctype(0,1,1)='yz'
  ctype(3,0,0)='xxx'
  ctype(0,3,0)='yyy'
  ctype(0,0,3)='zzz'
  ctype(2,1,0)='xxy'
  ctype(2,0,1)='xxz'
  ctype(1,2,0)='xyy'
  ctype(0,2,1)='yyz'
  ctype(1,0,2)='xzz'
  ctype(0,1,2)='yzz'
  ctype(1,1,1)='xyz'
!----------------------------------------------------------------------c
!     retrieve information about the most demanding shell block        c
!----------------------------------------------------------------------c
  CALL iosys('read integer maxprm from rwf',1,maxprm,0,' ')
  CALL iosys('read integer maxcont from rwf',1,mxcont,0,' ')
  CALL iosys('read integer maxl from rwf',1,maxl,0,' ')
  CALL iosys('read integer maxblk from rwf',1,maxblk,0,' ')
  CALL iosys('read integer d1maxblk from rwf',1,dlen,0,' ')
  CALL iosys('read integer dolp from rwf',1,dolp,0,' ')
!----------------------------------------------------------------------c
!     read in basis set information from read-write file               c
!----------------------------------------------------------------------c
  CALL iosys('read real exponents from rwf',-1,ex,0,' ')
  CALL iosys('read real "contraction coefficients" from rwf', -1,cont,0,' ')
  CALL iosys('read real "nuclear charges" from rwf',-1,zan,0,' ')
  CALL iosys('read real coordinates from rwf',-1,coords,0,' ')
  CALL iosys('read integer "pointer to primitives" from rwf', -1,ptprim,0,' ')
  CALL iosys('read integer "number of primitives" from rwf', -1,noprim,0,' ')
  CALL iosys('read integer "pointer to contraction coefficients"'//  &
             ' from rwf',-1,ptcont,0,' ')
  CALL iosys('read integer "number of contraction coefficients" '//  &
             'from rwf',-1,nocont,0,' ')
  CALL iosys('read integer "number of cartesians" from rwf', -1,nocart,0,' ')
  CALL iosys('read integer "number of pure functions" from rwf', -1,nobf,0,' ')
  CALL iosys('read integer "minimum momentum" from rwf', -1,minmom,0,' ')
  CALL iosys('read integer "maximum momentum" from rwf', -1,maxmom,0,' ')
  CALL iosys('read integer "pointer to cartesians" from rwf', -1,mintyp,0,' ')
  CALL iosys('read integer "pointer to first function" from rwf',  &
             -1,start,0,' ')
  CALL iosys('read integer "power of x" from rwf',-1,nx,0,' ')
  CALL iosys('read integer "power of y" from rwf',-1,ny,0,' ')
  CALL iosys('read integer "power of z" from rwf',-1,nz,0,' ')
!----------------------------------------------------------------------c
!             write out the location and charge of the atoms           c
!----------------------------------------------------------------------c
  CALL iosys ('write real "nuclear charges" to rmtrx',nat,zan, 0,' ')
  DO  i=1,nat
      CALL iosys ('write real "x-y-z atom-'//itoc(i)//'" to rmtrx',  &
                   3,coords(1,i),0,' ')
  END DO

  nstart(1) = 1
  ikount = 0
  iloc = 0
  IF(idrop /= 0) THEN
     WRITE(iout,100)
     CALL iosys('read integer "packing index vector" from rwf',  &
                 oldnbf,INDEX,0,' ')
  ELSE
     DO  i=1,oldnbf
         INDEX(i)=i
     END DO
  END IF
  kbf=0
  DO  iatom=1,nat
      WRITE(iout,101) iatom, zan(iatom)
      DO  itype=1,nbtype
          IF (noprim(iatom,itype) > 0) THEN
              mini=mintyp(itype)
              maxi=mintyp(itype)+nocart(itype)-1
              imax=maxmom(itype)
              WRITE(iout,102) itype,imax
              nprimi=noprim(iatom,itype)
              nconti=nocont(iatom,itype)
              WRITE(iout,103) nprimi, nconti
!  build arrays for quad codes: nstart, nstop, lnew, mnew, nnew, eta
!  the way the  indexing for quad codes works is that there are
!  ncontra contracted functions.  for the ith contracted function nstart(i)
!  and nstop(i) give the first and last locations of the relevant quantities
!  in the arrays lnew(i), mnew(i), nnew(i), and eta(i,j=1,5)
!  lnew, mnew, and new are the powers of x, y, and z of the primitive
!  and eta(i,j=1,5) contains the center, exponent and contraction coefficient
      
              DO  icont=1,nconti
                  DO  m=mini,maxi
                      kbf=kbf+1
                      IF(INDEX(kbf) /= 0) THEN
                         ikount = ikount + 1
                         DO  iii=1,nprimi
                             iiirel = iii-1+ (icont-1)*nprimi
                             IF (cont(ptcont(iatom,itype)+iiirel) > 1.e-20) THEN
                                 iloc = iloc+1
                                 cnt(iloc)=iatom
                                 eta(iloc,1) = coords(1,iatom)
                                 eta(iloc,2) = coords(2,iatom)
                                 eta(iloc,3) = coords(3,iatom)
                                 eta(iloc,4) = ex(ptprim(iatom,itype) +iii - 1)
                                 eta(iloc,5) = cont(ptcont(iatom,itype) +iiirel)
                                 lnew(iloc) = nx(m)
                                 mnew(iloc) = ny(m)
                                 nnew(iloc) = nz(m)
                             END IF
                         END DO
                         nstop(ikount) = iloc
                         nstart(ikount+1) = nstop(ikount)+ 1
                      END IF
                  END DO
              END DO
              IF(prnloc(1)) THEN
                 WRITE(iout,104)
                 CALL matout(ex(ptprim(iatom,itype)),nprimi,1, nprimi,1,iout)
              END IF
              IF(prnloc(1)) THEN
                 WRITE(iout,105)
                 CALL matout(cont(ptcont(iatom,itype)),nprimi,nconti,  &
                             nprimi,nconti,iout)
              END IF
              IF(prnloc(1)) THEN
                 WRITE(iout,106)
                 DO  m=mini,maxi
                     WRITE(iout,107) nx(m), ny(m), nz(m)
                 END DO
              END IF
          END IF
      END DO
  END DO
  npr = iloc
  ncon = ikount
! write the newly indexed basis set information to kohndt.
  CALL iosys ('write integer "no. primitives" to rmtrx',1, npr,0,' ')
  CALL iosys ('write integer "no. contracted" to rmtrx',1, ncon,0,' ')
  CALL iosys ('write integer start to rmtrx',ncon,nstart, 0,' ')
  CALL iosys ('write integer stop to rmtrx',ncon,nstop, 0,' ')
  CALL iosys ('write integer "x power" to rmtrx',npr,lnew,0,' ')
  CALL iosys ('write integer "y power" to rmtrx',npr,mnew,0,' ')
  CALL iosys ('write integer "z power" to rmtrx',npr,nnew,0,' ')
  CALL iosys ('write real eta to rmtrx',1500,eta,0,' ')
  CALL iosys ('write integer atom to rmtrx',npr,cnt,0,' ')
  WRITE(iout,108)
  IF(prnloc(2)) THEN
     WRITE (iout,109)
     WRITE (iout,110)
     DO  i=1,ncon
         nlow = nstart(i)
         nhi = nstop(i)
         DO  j=nlow,nhi
             WRITE(iout,111) i, j, ctype(lnew(j),mnew(j),nnew(j)),  &
             cnt(j), (eta(j,jj),jj=4,5)
         END DO
     END DO
  END IF
100 FORMAT (/,15X,'drop index vector read from rwf')
101 FORMAT (/,5X,'atom number',1X,i4,2X,'charge',f15.8)
102 FORMAT (/,5X,'type primitive',1X,i2,2X,'maximum angular momentum', 1X,i2)
103 FORMAT (/,5X,'no. primitives',1X,i3,2X,'no. contracted',1X,i3)
104 FORMAT(/,5X,'primitive exponents')
105 FORMAT(/,5X,'contraction coefficients')
106 FORMAT (/,5X,'xyz convention')
107 FORMAT (/,5X,'total no. primitive functions',1X,i4,1X,'total no. c  &
    ontracted functions',1X,i4)
108 FORMAT (/,5X,'basis set information written to kohndat')
109 FORMAT (/,5X,'basis set in new format')
110 FORMAT(/,5X,'con fn',3X,'ao',2X,'sym',4X,'cen',9X,'exp',11X,'coef' )
111 FORMAT(/,7X,i3,3X,i3,3X,a3,4X,i2,3X,f12.5,2X,f12.5)
END SUBROUTINE basout

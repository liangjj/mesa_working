!********************************************************************************
!********************************************************************************
                        MODULE Points_Weights
  USE Data
  USE Grid_Defined_Types
  IMPLICIT NONE
!********************************************************************************
!********************************************************************************
!
                                 Contains
!
!********************************************************************************
!********************************************************************************
!deck lebdata.f
!***begin prologue     lebdata
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           xyz, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            
!***references

!***routines called
!***end prologue       lebdata
  SUBROUTINE lebdata(grd)
  IMPLICIT NONE
  TYPE(GRID                     :: grd
  INTEGER, PARAMETER            :: dim_1=26
  INTEGER, PARAMETER            :: dim_2=96
  INTEGER, PARAMETER            :: dim_3=24
  INTEGER, PARAMETER            :: dim_4=48
  INTEGER                       :: ii
  INTEGER                       :: iend
  REAL(idp), DIMENSION(4)       :: l
  REAL(idp), DIMENSION(4)       :: m
  REAL(idp), DIMENSION(4)       :: biga
  REAL(idp), DIMENSION(4)       :: bigb
  REAL(idp), DIMENSION(4)       :: bigc
  REAL(idp), DIMENSION(4)       :: bigd
  REAL(idp), DIMENSION(4)       :: x
  REAL(idp), DIMENSION(1)       :: p
  REAL(idp), DIMENSION(1)       :: q
  REAL(idp), DIMENSION(1)       :: r
  REAL(idp), DIMENSION(1)       :: u
  REAL(idp), DIMENSION(1)       :: w
  REAL(idp), DIMENSION(3,dim_1) :: v
  REAL(idp), DIMENSION(3,dim_2) :: lm
  REAL(idp), DIMENSION(3,dim_3) :: pq
  REAL(idp), DIMENSION(3,dim_4) :: ruw
!
  CALL Get_V_Table(v, ii, iend)
  grd%angleb(:,1:iend) = v(:,1:iend)
!
  m(1:4)=0.D0
  l(1:4)=0.D0
  IF (grd%nleb == 9) THEN
      biga(1)=.00952380952387D0
      biga(3)=.0321428571429D0
      bigc(1)=.0285714285714D0
      p(1)=.888073833977D0
      q(1)=.459700843381D0
      grd%wtleb(1:6)=biga(1)
      grd%wtleb(7:14)=biga(3)
      CALL Get_PQ_Table (pq  p, q )
      grd%angleb(:,ii+9:ii+32) = pq(:,1:24)
      grd%wtleb(ii+9,ii+32) = bigc(1)
  END IF
  IF (grd%nleb >= 11 ) THEN
      CALL Get_LM_Table(lm, l, m, grd%nleb)
  END IF  
  IF (grd%nleb >= 23 ) THEN
      CALL Get_RUW_Table(ruw, r, u, w)
  END IF  
  IF (grd%nleb == 11 ) THEN
      m(1)=.90453403373329086794D0
      l(1)=.30151134457800000000D0
      bigc(1)=0.d0
      biga(1)=.01269841269841269841D0
      biga(2)=.02257495590828924162D0
      biga(3)=.02109375000000000000D0
      bigb(1)=.02017333553791887125D0
      bigb(2)=0.d0
      bigb(3)=0.d0
      bigb(4)=0.d0
      grd%wtleb(1:6)=biga(1)
      grd%wtleb(7:18)=biga(2)
      grd%wtleb(19:26)=biga(3)
      grd%angleb(:,ii+9:ii+32) = lm(:,1:24)
      grd%wtleb(ii+9:ii+32) = bigb(1)
  ELSE IF ( grd%nleb == 13 ) THEN
      m(1)=.73379938570534280703D0
      l(1)=.48038446141526140045D0
      bigc(1)=.01652217099371570916D0
      biga(1)=.00051306717973384640D0
      biga(2)=.01660406956574203960D0
      biga(3)=-.02958603896103896103D0
      bigb(1)=.02657620708215946311D0
      p(1)=.94715622136258788774D0
      q(1)=.32077264898077643273D0
      grd%angleb(:,ii+9:ii+32) = lm(:,1:24)
      grd%wtleb(ii+9:ii+32) = bigb(1)
      grd%angleb(:,ii+33:ii+56) = pq(:,1:24)
      grd%wtleb(:,ii+33:ii+56) = bigc(1)
  ELSE IF (grd%nleb == 17) THEN
      n1=3
      bigc(1)=.00969499636166D0
      biga(1)=.00382827049494D0
      biga(2)=0.d0
      biga(3)=.00979373751249D0
      bigb(1)=.00821173728319D0
      bigb(2)=.00959547133607D0
      bigb(3)=.00994281489118D0
      bigb(4)=0.d0
      m(1)=.965124035087D0
      l(1)=.185115635345D0
      m(2)=.828769981253D0
      l(2)=.395689473056D0
      m(3)=.215957291846D0
      l(3)=.690421048382D0
      p(1)=.878158910604D0
      q(1)=.478369028812D0
      grd%angleb(:,ii+9:ii+80) = lm(:,1:72)
      grd%wtleb(ii+9:ii+32) = bigb(1)
      grd%wtleb(ii+33:ii+56) = bigb(2)
      grd%wtleb(ii+57:ii+80) = bigb(3)
      grd%angleb(:,ii+81:ii+104) = pq(:,1:24)
      grd%wtleb(:,ii+81:ii+104) = bigc(1)
  ELSE IF (grd%nleb == 23) THEN
      n1=4
      n2=1
      n3=1
      bigc(1)=5.05184606462D-03
      p(1)=.938319218138D0
      q(1)=.345770219761D0
      bigd(1)=5.53024891623D-03
      r(1)=.836036015482D0
      u(1)=.159041710538D0
      w(1)=.525118572443D0
      biga(1)=1.78234044724D-03
      biga(2)=5.71690594998D-03
      biga(3)=5.57338317884D-03
      bigb(1)=5.51877146727D-03
      bigb(2)=5.15823771181D-03
      bigb(3)=5.60870408259D-03
      bigb(4)=4.10677702817D-03
      m(1)=.777493219315D0
      l(1)=.444693317871D0
      m(2)=.912509096867D0
      l(2)=.289246562758D0
      m(3)=.314196994183D0
      l(3)=.671297344270D0
      m(4)=.982972302707D0
      l(4)=.129933544765D0
      grd%angleb(:,ii+9:ii+104) = lm(:,1:96)
      grd%wtleb(ii+9:ii+32) = bigb(1)
      grd%wtleb(ii+33:ii+56) = bigb(2)
      grd%wtleb(ii+57:ii+80) = bigb(3)
      grd%wtleb(ii+81:ii+104) = bigb(4)
      grd%angleb(:,ii+105:ii+128) = pq(:,1:24)
      grd%wtleb(:,ii+129:ii+176) = bigd(1)
  END IF
  ELSE IF (grd%nleb == 11) THEN

      grd%wtleb(1:6)=biga(1)
      grd%wtleb(7:18)=biga(2)
      grd%wtleb(19:26)=biga(3)


      getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
  ELSE IF (grd%nleb == 13) THEN

           DO i=1,6
              grd%wtleb(i)=biga(1)
              grd%wtleb(i+6)=biga(2)
              grd%wtleb(i+12)=biga(2)
              grd%wtleb(i+18)=biga(3)
          END DO
          grd%wtleb(25)=biga(3)
          grd%wtleb(26)=biga(3)
          CALL getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
          CALL getcnd(grd%angleb(1,ii+33),grd%wtleb(ii+33),p,q,bigc)
  ELSE IF (grd%nleb == 17) THEN
           n1=3
           bigc(1)=.00969499636166D0
           biga(1)=.00382827049494D0
           biga(2)=0.d0
           biga(3)=.00979373751249D0
           bigb(1)=.00821173728319D0
           bigb(2)=.00959547133607D0
           bigb(3)=.00994281489118D0
           bigb(4)=0.d0
           m(1)=.965124035087D0
           l(1)=.185115635345D0
           m(2)=.828769981253D0
           l(2)=.395689473056D0
           m(3)=.215957291846D0
           l(3)=.690421048382D0
           p(1)=.878158910604D0
           q(1)=.478369028812D0
           DO i=1,6
              grd%wtleb(i)=biga(1)
           END DO
           DO i=7,14
              grd%wtleb(i)=biga(3)
           END DO
           CALL getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
           CALL getcnd(grd%angleb(1,ii+81),grd%wtleb(ii+81),p,q,bigc)
  ELSE IF (grd%nleb == 23) THEN
           n1=4
           n2=1
           n3=1
           bigc(1)=5.05184606462D-03
           p(1)=.938319218138D0
           q(1)=.345770219761D0
           bigd(1)=5.53024891623D-03
           r(1)=.836036015482D0
           u(1)=.159041710538D0
           w(1)=.525118572443D0
           biga(1)=1.78234044724D-03
           biga(2)=5.71690594998D-03
           biga(3)=5.57338317884D-03
           bigb(1)=5.51877146727D-03
           bigb(2)=5.15823771181D-03
           bigb(3)=5.60870408259D-03
           bigb(4)=4.10677702817D-03
           m(1)=.777493219315D0
           l(1)=.444693317871D0
           m(2)=.912509096867D0
           l(2)=.289246562758D0
           m(3)=.314196994183D0
           l(3)=.671297344270D0
           m(4)=.982972302707D0
           l(4)=.129933544765D0
           DO i=1,6
              grd%wtleb(i)=biga(1)
           END DO
           DO i=7,18
              grd%wtleb(i)=biga(2)
           END DO
           DO i=19,26
              grd%wtleb(i)=biga(3)
           END DO
           CALL getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
           CALL getcnd(grd%angleb(1,ii+105),grd%wtleb(ii+105),p,q,bigc)
           CALL getdnd(grd%angleb(1,ii+129),grd%wtleb(ii+129),r,u,w,bigd)
  END IF
  IF (grd%nleb == 9) THEN

!********************************************************************************
  END SUBROUTINE lebdata
!********************************************************************************
!deck Get_V_Table
!***begin prologue     Get_V_Table
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           xyz, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            
!***references

!***routines called
!***end prologue       Get_V_Table
  SUBROUTINE Get_V_Table(v,ii,iend
  IMPLICIT NONE
  INTEGER                       :: ii
  INTEGER                       :: iend
  REAL(idp), DIMENSION(:,:)     ::v
!
  v(:,1)=[0.d0,0.d0,1.d0]
  v(:,2)=[0.d0,0.d0,-1.d0]
  v(:,3)=[0.d0,1.d0,0.d0]
  v(:,4)=[0.d0,-1.d0,0.d0]
  v(:,5)=[1.d0,0.d0,0.d0]
  v(:,6)=[-1.d0,0.d0,0.d0]
  iend = 6
!
  IF (grd%nleb == 11.OR.grd%nleb == 13.OR.grd%nleb == 23) THEN
      v(:,7)=[sqrt_half,sqrt_half,0.d0]
      v(:,8)=[sqrt_half,-sqrt_half,0.d0]
      v(:,9)=[-sqrt_half,sqrt_half,0.d0]
      v(:,10)=[-sqrt_half,-sqrt_half,0.d0]
      v(:,11)=[sqrt_half,0.d0,sqrt_half]
      v(:,12)=[sqrt_half,0.d0,-sqrt_half]
      v(:,13)=[-sqrt_half,0.d0,sqrt_half]
      v(:,14)=[-sqrt_half,0.d0,-sqrt_half]
      v(:,15)=[0.d0,sqrt_half,sqrt_half]
      v(:,16)=[0.d0,sqrt_half,-sqrt_half]
      v(:,17)=[0.d0,-sqrt_half,sqrt_half]
      v(:,18)=[0.d0,-sqrt_half,-sqrt_half]
      ii = 18
!
  ELSE
      ii = 6
  END IF
  v(:,ii+1)=[sqrt_third,sqrt_third,sqrt_third]
  v(:,ii+2)=[sqrt_third,sqrt_third,-sqrt_third]
  v(:,ii+3)=[sqrt_third,-sqrt_third,sqrt_third]
  v(:,ii+4)=[sqrt_third,-sqrt_third,-sqrt_third]
  v(:,ii+5)=[-sqrt_third,sqrt_third,sqrt_third]
  v(:,ii+6)=[-sqrt_third,sqrt_third,-sqrt_third]
  v(:,ii+7)=[-sqrt_third,-sqrt_third,sqrt_third]
  v(:,ii+8)=[-sqrt_third,-sqrt_third,-sqrt_third]
  iend = ii + 8
!********************************************************************************
  END SUBROUTINE Get_V_Table
!********************************************************************************
!deck Get_PQ_Table
!***begin prologue     Get_PQ_Table
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           xyz, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            
!***references

!***routines called
!***end prologue       Get_PQ_Table
  SUBROUTINE Get_PQ_Table(pq,p,q)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)  ::pq
  REAL(idp), DIMENSION(:)    ::p
  REAL(idp), DIMENSION(:)    ::q
!
  pq(:,1) = [  p(1), q(1), 0.d0 ]
  pq(:,2) = [  p(1),-q(1), 0.d0 ]
  pq(:,3) = [ -p(1), q(1), 0.d0 ]
  pq(:,4) = [ -p(1),-q(1), 0.d0 ]
  pq(:,5) = [  p(1), 0.d0, q(1) ]
  pq(:,6) = [  p(1), 0.d0,-q(1) ]
  pq(:,7) = [ -p(1), 0.d0, q(1) ]
  pq(:,8) = [ -p(1), 0.d0,-q(1) ]
!
  pq(:,9) =  [ 0.d0, p(1), q(1) ]
  pq(:,10) = [ 0.d0  p(1),-q(1) ]
  pq(:,11) = [ 0.d0,-p(1), q(1) ]
  pq(:,12) = [ 0.d0,-p(1),-q(1) ]
  pq(:,13) = [ q(1), p(1), 0.d0 ]
  pq(:,14)  =[ q(1),-p(1), 0.d0 ]
  pq(:,15)  =[-q(1), p(1), 0.d0 ]
  pq(:,16) = [-q(1),-p(1), 0.d0 ]
!
  pq(:,17) = [ q(1), 0.d0, p(1) ]
  pq(:,18) = [ q(1), 0.d0,-p(1) ]
  pq(:,19) = [-q(1), 0.d0, p(1) ]
  pq(:,20) = [-q(1), 0.d0,-p(1) ]
  pq(:,21) = [ 0.d0, q(1), p(1) ]
  pq(:,22) = [ 0.d0, q(1),-p(1) ]
  pq(:,23) = [ 0.d0,-q(1),p(1) ]
  pq(:,24) = [ 0.d0,-q(1),-p(1) ]
!********************************************************************************
  END SUBROUTINE Get_PQ_Table
!********************************************************************************
!********************************************************************************
!deck Get_LM_Table
!***begin prologue     Get_LM_Table
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           xyz, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            
!***references

!***routines called
!***end prologue       Get_LM_Table
  SUBROUTINE Get_LM_Table(lm,l,m)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)  ::y
  REAL(idp), DIMENSION(:)    ::l
  REAL(idp), DIMENSION(:)    ::m
!
!
  lm(:,1) = [  l(1), l(1), m(1) ]
  lm(:,2) = [  l(1), l(1),-m(1) ]
  lm(:,3) = [  l(1),-l(1), m(1) ]
  lm(:,4) = [  l(1),-l(1),-m(1) ]
  lm(:,5) = [ -l(1), l(1), m(1) ]
  lm(:,6) = [ -l(1), l(1),-m(1) ]
  lm(:,7) = [ -l(1),-l(1), m(1) ]
  lm(:,8) = [ -l(1),-l(1),-m(1) ]
!
!       eight=(lml)
!
  lm(:,9) =  [ l(1), m(1), l(1) ]
  lm(:,10) = [ l(1), m(1),-l(1) ]
  lm(:,11) = [ l(1),-m(1), l(1) ]
  lm(:,12) = [ l(1),-m(1),-l(1) ]
  lm(:,13) = [-l(1), m(1), l(1) ]
  lm(:,14) = [-l(1), m(1),-l(1) ]
  lm(:,15) = [-l(1),-m(1), l(1) ]
  lm(:,16) = [-l(1),-m(1),-l(1) ]
!
!       eight=(mll)
!
  lm(:,17) = [ m(1), l(1), l(1) ]
  lm(:,18) = [ m(1), l(1),-l(1) ]
  lm(:,19) = [ m(1),-l(1), l(1) ]
  lm(:,20) = [ m(1),-l(1),-l(1) ]
  lm(:,21) = [-m(1), l(1), l(1) ]
  lm(:,22) = [-m(1), l(1),-l(1) ]
  lm(:,23) = [-m(1),-l(1), l(1) ]
  lm(:,24) = [-m(1),-l(1),-l(1) ]
!
!         eight (l,l,m)
!
  lm(:,25) = [ l(2), l(2), m(2) ]
  lm(:,26) = [ l(2), l(2),-m(2) ]
  lm(:,27) = [ l(2),-l(2), m(2) ]
  lm(:,28) = [ l(2),-l(2),-m(2) ]
  lm(:,29) = [-l(2), l(2), m(2) ]
  lm(:,30) = [-l(2), l(2),-m(2) ]
  lm(:,31) = [-l(2),-l(2), m(2) ]
  lm(:,32) = [-l(2),-l(2),-m(2) ]
!
!
!         eight (l,m,l)
!
  lm(:,33) = [ l(2), m(2), l(2) ]
  lm(:,34) = [ l(2), m(2),-l(2) ]
  lm(:,35) = [ l(2),-m(2), l(2) ]
  lm(:,36) = [ l(2),-m(2),-l(2) ]
  lm(:,37) = [-l(2), m(2), l(2) ]
  lm(:,38) = [-l(2), m(2),-l(2) ]
  lm(:,39) = [-l(2),-m(2), l(2) ]
  lm(:,40) = [-l(2),-m(2),-l(2) ]
!
!         eight (m,l,l)
!
  lm(:,41) = [ m(2), l(2), l(2) ]
  lm(:,42) = [ m(2), l(2),-l(2) ]
  lm(:,43) = [ m(2),-l(2), l(2) ]
  lm(:,44) = [ m(2),-l(2),-l(2) ]
  lm(:,45) = [-m(2), l(2), l(2) ]
  lm(:,46) = [-m(2), l(2),-l(2) ]
  lm(:,47) = [-m(2), -l(2), l(2)]
  lm(:,48) = [-m(2), -l(2),-l(2)]
!
!         eight (l,l,m)
!
  lm(:,49) = [ l(3), l(3), m(3) ]
  lm(:,50) = [ l(3), l(3),-m(3) ]
  lm(:,51) = [ l(3),-l(3), m(3) ]
  lm(:,52) = [ l(3),-l(3),-m(3) ]
  lm(:,53) = [-l(3), l(3), m(3) ]
  lm(:,54) = [-l(3), l(3),-m(3) ]
  lm(:,55) = [-l(3),-l(3), m(3) ]
  lm(:,56) = [-l(3),-l(3),-m(3) ]
!
!         eight (l,m,l)
!
  lm(:,57) = [ l(3), m(3), l(3) ]
  lm(:,58) = [ l(3), m(3),-l(3) ]
  lm(:,59) = [ l(3),-m(3), l(3) ]
  lm(:,60) = [ l(3),-m(3),-l(3) ]
  lm(:,61) = [-l(3), m(3), l(3) ]
  lm(:,62) = [-l(3), m(3),-l(3) ]
  lm(:,63) = [-l(3),-m(3), l(3) ]
  lm(:,64) = [-l(3),-m(3),-l(3) ]
!
!         eight (m,l,l)
!
  lm(:,65) = [ m(3), l(3), l(3) ]
  lm(:,66) = [ m(3), l(3),-l(3) ]
  lm(:,67) = [ m(3),-l(3), l(3) ]
  lm(:,68) = [ m(3),-l(3),-l(3) ]
  lm(:,69) = [-m(3), l(3), l(3) ]
  lm(:,70) = [-m(3), l(3),-l(3) ]
  lm(:,71) = [-m(3),-l(3), l(3) ]
  lm(:,72) = [-m(3),-l(3),-l(3) ]
!
!
!         eight (l,l,m)

  lm(:,73) = [ l(4), l(4), m(4) ]
  lm(:,74) = [ l(4), l(4),-m(4) ]
  lm(:,75) = [ l(4),-l(4), m(4) ]
  lm(:,76) = [ l(4),-l(4),-m(4) ]
  lm(:,77) = [-l(4), l(4), m(4) ]
  lm(:,78) = [-l(4), l(4),-m(4) ]
  lm(:,79) = [-l(4),-l(4), m(4) ]
  lm(:,80) = [-l(4),-l(4),-m(4) ]
!
!         eight (l,m,l)
!
  lm(:,81) = [ l(4), m(4), l(4) ]
  lm(:,82) = [ l(4), m(4),-l(4) ]
  lm(:,83) = [ l(4),-m(4), l(4) ]
  lm(:,84) = [ l(4),-m(4),-l(4) ]
  lm(:,85) = [-l(4), m(4), l(4) ]
  lm(:,86) = [-l(4), m(4),-l(4) ]
  lm(:,87) = [-l(4),-m(4), l(4) ]
  lm(:,88) = [-l(4),-m(4),-l(4) ]
!
!         eight (m,l,l)
!
  lm(:,89) = [ m(4), l(4), l(4) ]
  lm(:,90) = [ m(4), l(4),-l(4) ]
  lm(:,91) = [ m(4),-l(4), l(4) ]
  lm(:,92) = [ m(4),-l(4),-l(4) ]
  lm(:,93) = [-m(4), l(4), l(4) ]
  lm(:,94) = [-m(4), l(4),-l(4) ]
  lm(:,95) = [-m(4),-l(4), l(4) ]
  lm(:,96) = [-m(4),-l(4),-l(4) ]
!********************************************************************************
  END SUBROUTINE Get_LM_Table
!********************************************************************************
!deck Quadrature.f
!***begin prologue     Quadrature
!***date written       9308021   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           shells
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            calculate angular and radial quadrature points
!***                   and weights
!***references
!***routines called    gaussq ( math )
!***end prologue       Quadrature

  SUBROUTINE Quadrature(cent,grd,sh)
!, TYPE,icen,cen,r,nr,eta,rpt,thpt,phpt,wtr,wtsum,  &
!                        wtth,wtph,angleb,work,wtleb,sthet,sphi,cphi,     &
!                        yuk,scr,cart_grid,wt,sph_grid,ns,nrmax,nthet,nphi,      &
!                        ngrid,nwts,ncen,ncplus,prnt,nleb,nang,nonsep,    &
!                        yukon,nodisk)
  IMPLICIT NONE
  TYPE(CENTERS)                                :: cent
  TYPE(GRID)                                   :: grd
  TYPE(SHELL), DIMENSION(:), ALLOCATABLE       :: sh
  INTEGER                      :: ntrad
  INTEGER                      :: ntwts
  INTEGER                      :: icen
  INTEGER                      :: lenth
  REAL(idp)                    :: cen(3,ncplus)
  REAL(idp)                    :: r(ns+1)
  INTEGER                      :: nr(ns)
  REAL(idp)                    :: eta(ncplus)
  REAL(idp)                    :: rpt(nrmax,ns)
  REAL(idp)                    :: thpt(nthet)
  REAL(idp)                    :: phpt(nphi)
  REAL(idp)                    :: wtr(*)
  REAL(idp)                    :: wtsum(*)
  REAL(idp)                    :: wtth(nthet)
  REAL(idp)                    :: wtph(nphi)
  REAL(idp)                    :: angleb(*)
  REAL(idp)                    :: work(*)
  REAL(idp)                    :: wtleb(*)
  REAL(idp)                    :: sthet(nthet)
  REAL(idp)                    :: sphi(nphi)
  REAL(idp)                    :: cphi(nphi)
  REAL(idp)                    :: yuk(*)
  REAL(idp)                    :: scr(*)
  REAL(idp)                    :: cart_grid(3,ngrid)
  REAL(idp)                    :: wt(nwts)
  REAL(idp)                    :: spgrid(3,ngrid)
  INTEGER                      :: ns
  INTEGER                      :: nrmax
  INTEGER                      :: nthet
  INTEGER                      :: nphi
  INTEGER                      :: ngrid
  INTEGER                      :: nwts
  INTEGER                      :: ncen
  INTEGER                      :: ncplus
  LOGICAL                      :: prnt
  INTEGER                      :: nleb
  INTEGER                      :: nang
  LOGICAL                      :: nonsep
  LOGICAL                      :: yukon
  LOGICAL                      :: nodisk

  INTEGER :: thefac, phifac

  CHARACTER (LEN=30) :: rtyp, thtyp, phtyp, str
  CHARACTER (LEN=3) :: itoc, chra

  DIMENSION  delthe(2), delphi(2)

  DIMENSION  dummy(2)

  lstrng=lenth(str)
  WRITE(iout,*)
  WRITE(iout,*) 'generating grid for '//str(1:lstrng)
!              Total size of radial points for this grid
  ALLOCATE(grd%sh(1:ns))
  DO i = 1, grd%nshell
     ALLOCATE(grd%sh(i)%rpt(1:grd%nr(i)), grd%sh(i)%wtr(1:grd%nrwts(i)))
  END DO
  nedges=grid%nshell+1
!
!                     create some files to hold data
!
  CALL iosys ('create real "radial points '//str//'" on grid',grd%ntrad,0,0,' ')
  CALL iosys ('create real "scaled radial weights '//str//'" on grid',grd%ntradwt, 0,0,' ')
  CALL iosys ('create real "unscaled radial weights '//str//'" on grid',grd%ntradwt, 0,0,' ')
  IF (grd%rtyp == 'newton-cotes') THEN
      CALL iosys ('create real "summed radial weights '//str//'" on grid',grd%ntrad,0,0,' ')
  END IF
  IF (.NOT.nodisk) THEN
      CALL iosys ('create real "yukawa potential '//str//' coordinates" on grid',grd%ngrid, 0,0,' ')
      CALL iosys ('create real "atomic grid '//str//'" on grid',3*grd%ngrid,0,0,' ')
      CALL iosys ('create real "spherical atomic grid '//str//'" on grid',3*grd%ngrid,0,0,' ')
      CALL iosys ('create real "unscaled atomic weights '//str//'" on grid',grd%nwts,0,0,' ')
  END IF
  IF (nonsep) THEN
!
!          in a lebedev quadrature the number of theta and phi points are
!          identical since the quadrature is non-separable. the combined
!          angular weight is stored in wtth.
!
      ALLOCATE(grd%angleb(1:3,1:grd%nang),grd%cthet(1:grd%nang),grd%sthet(1:grd%nang),     &
               grd%phpt(1:grd%nang),grd%sphi(1:grd%nang),grd%cphi(1:grd%nang),             &
               grd%wtleb(1:grd%nang))   
      CALL lebdev(angleb,angleb,thpt,sthet,phpt,sphi,cphi,wtleb, grd%nleb,grd%nang,str)
  
  ELSE
      CALL iosys('read character "theta quadrature type '//str//'" from gridt',-1,0,0,grd%thtyp)
      CALL iosys('read character "phi quadrature type '//str//'" from grid',-1,0,0,grd%phtyp)
      difth=grd%range(2)-grd%range(1)
      difph=grd%range(4)-grd%range(3)
!     multiplicative factors for theta and phi integration if symmetry
!     allows integrals to be calculated on subdomains.
      thefac=pi/difth
      phifac=2.d0*pi/difph
      WRITE(iout,1) thefac, phifac
      tmp=COS(grd%range(1))
      ampt=tmp
      abpt=COS(grd%range(2))
      ampt=(ampt-abpt)*.5D0
      abpt=(abpt+tmp)*.5D0
      CALL gaussq(grd%thtyp,grd%nthet,0.d+00,0.d+00,0,dummy,scr,thpt,grd%wtth)
  DO  i=1,grd%nthet
    thpt(i)=ampt*thpt(i)+abpt
    sthet(i)=SQRT((1.d0-thpt(i)*thpt(i)))
    grd%wtth(i)=thefac*ampt*grd%wtth(i)
  END DO
  ampp=(grd%range(4)-grd%range(3))*.5D0
  abpp=(grd%range(4)+grd%range(31))*.5D0
  CALL gaussq(grd%phtyp,grd%nphi,0.d+00,0.d+00,0,dummy,scr, phpt,grd%wtph)
  DO  i=1,grd%nphi
    phpt(i)=ampp*phpt(i)+abpp
    sphi(i)=SIN(phpt(i))
    cphi(i)=COS(phpt(i))
    grd%wtph(i)=phifac*ampp*grd%wtph(i)
  END DO
  CALL iosys ('write real "theta points '//str//  &
      '" to lamdat',grd%nthet,thpt,0,' ')
  CALL iosys ('write real "theta weights '//str//  &
      '" to lamdat',grd%nthet,grd%wtth,0,' ')
  CALL iosys ('write real "phi points '//str// '" to lamdat',grd%nphi,phpt,0,' ')
  CALL iosys ('write real "phi weights '//str// '" to lamdat',grd%nphi,grd%wtph,0,' ')
  grd%nang=grd%nthet*grd%nphi
END IF
IF (prnt) THEN
  WRITE(iout,*)
  WRITE(iout,*) 'cos(theta) points and weights'
  WRITE(iout,2) (thpt(ii),ii=1,grd%nthet)
  WRITE(iout,2) (grd%wtth(ii),ii=1,grd%nthet)
  WRITE(iout,*)
  WRITE(iout,*) 'phi points and weights'
  WRITE(iout,2) (phpt(ii),ii=1,grd%nphi)
  WRITE(iout,2) (grd%wtph(ii),ii=1,grd%nphi)
END IF

!         do the radial quadrature points. they are divided into shells.


ngcnt=0
nwtcnt=0
aprvol=0.d0
yukint=0.d0
CALL rzero(yuk,grd%ngrid)
icount=0
DO  i=1,ns
  IF (r(i) == 0.d0) THEN
    r(i)=1.d-10
  END IF
  locwtr=nrmax*(i-1)+1
!     get radial quadrature points and weights for this shell
  IF (grd%rtyp == 'legendre') THEN
    ampr=(r(i+1)-r(i))*.5D0
    abpr=(r(i+1)+r(i))*.5D0
    CALL gaussq(grd%rtyp,grd%nr(i),0.d+00,0.d+00,0,dummy,scr, rpt(1,i),wtr(locwtr))
    grd%nwtsr=grd%nr(i)
    icount=icount+grd%nr(i)
    loccnt=locwtr
    DO  j=1,grd%nr(i)
      numr=numr+1
      rpt(j,i)=ampr*rpt(j,i)+abpr
      wtr(loccnt)=ampr*wtr(loccnt)
      loccnt=loccnt+1
    END DO
  ELSE IF(grd%rtyp == 'newton-cotes') THEN
    locwtr=nrmax*(nrmax-1)*(i-1)+1
    CALL necote(r(i),r(i+1),rpt(1,i),wtr(locwtr),grd%nr(i), .false.)
    grd%nwtsr=(grd%nr(i)-1)*grd%nr(i)
    icount=icount+grd%nr(i)
  ELSE
    CALL lnkerr('error in radial quadrature')
  END IF
  CALL iosys ('write real "radial points '//str//  &
      '" to lamdat without rewinding',grd%nr(i), rpt(1,i),0,' ')
  CALL iosys ('write real "unscaled radial weights '//  &
      str//'" to lamdat without '// 'rewinding',grd%nwtsr,wtr(locwtr),0,' ')
  IF(grd%rtyp == 'newton-cotes') THEN
    CALL sumncw(wtr(locwtr),wtsum,grd%nr(i))
    CALL iosys ('write real "summed radial weights '//  &
        str//'" to lamdat without rewinding', grd%nr(i),wtsum,0,' ')
  END IF
!      the procedure is different if the angular quadrature is separable or
!      non-separable in  (theta,phi)
  CALL mkgr(cart_grid(1,ngcnt+1),sph_grid(1,ngcnt+1),rpt(1,i),thpt,  &
      sthet,sphi,cphi,cen(1,icen),grd%nr(i),grd%nthet, grd%nphi,grd%nang,nonsep)
  IF (yukon) THEN
    CALL yukawa(yuk(ngcnt+1),cart_grid(1,ngcnt+1),eta,cen,  &
        grd%nr(i),grd%nthet,grd%nphi,ncen,grd%nang,nonsep)
  END IF
  CALL mkwt(wtr(locwtr),wtr(locwtr),grd%wtth,grd%wtph,wtleb,  &
      wt(nwtcnt+1),work,grd%nr(i),grd%nthet,grd%nphi,grd%rtyp,grd%nang,nonsep)
  IF (yukon) THEN
    CALL mkyunt(rpt(1,i),wt(nwtcnt+1),yuk(ngcnt+1),aprvol,  &
        yukint,work,grd%nr(i),grd%nthet,grd%nphi,grd%rtyp,grd%nang, nonsep)
  END IF
  ngcnt=ngcnt+grd%nr(i)*grd%nang
  nwtcnt=nwtcnt+grd%nwtsr*grd%nang
  CALL scalwt(rpt(1,i),wtr(locwtr),wtr(locwtr),work,grd%nr(i),grd%rtyp)
  CALL iosys ('write real "scaled radial weights '//  &
      str//'" to lamdat without '// 'rewinding',grd%nwtsr,wtr(locwtr),0,' ')
  IF (prnt) THEN
    WRITE(iout,*) '     data on radial points for shell = ',i
    WRITE(iout,*)
    WRITE(iout,*) 'radial points and weights'
    WRITE(iout,2) (rpt(ii,i),ii=1,grd%nr(i))
    WRITE(iout,2) (wtr(ii),ii=locwtr,locwtr+grd%nwtsr-1)
  END IF
END DO
IF (ngcnt /= grd%ngrid) THEN
  CALL lnkerr('error in grid point count')
END IF
IF (nwtcnt /= grd%nwts) THEN
  CALL lnkerr('error in grid weight count')
END IF
CALL iosys ('write integer "total number of radial '//  &
    'points '//str//'" to lamdat',1,icount,0,' ')
IF (.NOT.nodisk) THEN
  CALL iosys ('write real "yukawa potential '//  &
      str//' coordinates" to lamdat',ngcnt, yuk,0,' ')
  CALL iosys ('write real "atomic grid '//  &
      str//'" to lamdat',3*ngcnt,cart_grid,0,' ')
  CALL iosys ('write real "spherical atomic grid '//  &
      str//'" to lamdat',3*ngcnt,sph_grid,0,' ')
  CALL iosys ('write real "unscaled atomic weights '//  &
      str//'" to lamdat',grd%nwts,wt,0,' ')
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
1 FORMAT(' multiplicative factors for theta and phi  :',2F10.5)
2 FORMAT(/4(2X,e15.8))
3 FORMAT(/,5X,'exact       value of volume =',e15.8)
4 FORMAT(/,5X,'approximate value of volume =',e15.8)
5 FORMAT(/,5X,'exact     value of yukawa integral = ',e15.8)
6 FORMAT(/,5X,'numerical value of yukawa integral = ',e15.8)
!********************************************************************************
!********************************************************************************
  END SUBROUTINE Quadrature
!********************************************************************************
!********************************************************************************
!deck Lebdev.f
!***begin prologue     Lebdev
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           lebdev, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            lebedev quadrature points and weights
!***routines called
!***end prologue       lebdev
!
  SUBROUTINE lebdev(pt,angle,cthet,sthet,phpt,sphi,cphi,wt,nleb, npts,str)
  IMPLICIT NONE
  TYPE(CENTERS)                                :: cent
  TYPE(GRID)                                   :: grd
  TYPE(SHELL), DIMENSION(:), ALLOCATABLE       :: sh
  REAL(idp), DIMENSION(:,:)                  :: pt(3,npts)
  REAL(idp), DIMENSION(:,:)                  :: angle(npts,3)
  REAL(idp), DIMENSION(:)                    :: cthet
  REAL(idp), DIMENSION(:)                    :: sthet
  REAL(idp), DIMENSION(:)                    :: phpt
  REAL(idp), DIMENSION(:)                    :: sphi
  REAL(idp), DIMENSION(:)                    :: cphi
  REAL(idp), DIMENSION(:)                    :: wt
  REAL(idp), DIMENSION(4)                    :: l
  REAL(idp), DIMENSION(4)                    :: m
  REAL(idp), DIMENSION(4)                    :: biga
  REAL(idp), DIMENSION(4)                    :: bigb
  REAL(idp), DIMENSION(4)                    :: bigc
  REAL(idp), DIMENSION(4)                    :: bigd
  REAL(idp), DIMENSION(1)                    :: p
  REAL(idp), DIMENSION(1)                    :: q
  REAL(idp), DIMENSION(1)                    :: r
  REAL(idp), DIMENSION(1)                    :: u
  REAL(idp), DIMENSION(1)                    :: w
  INTEGER                                    :: nleb
  INTEGER                                    :: npts
  REAL(idp)                                  :: sumwt
!
!
  ALLOCATE(grd%angleb(1:3,1:grd%nang),grd%cthet(1:grd%nang),grd%sthet(1:grd%nang),     &
           grd%phpt(1:grd%nang),grd%sphi(1:grd%nang),grd%cphi(1:grd%nang),             &
           grd%wtleb(1:grd%nang))   
!
!              the a1 nodes-these are independent of the quadrature
  CALL xyz(grd,0.d0,0.d0,1.d0,1)
  CALL xyz(grd,0.d0,0.d0,-1.d0,2)
  CALL xyz(grd,0.d0,1.d0,0.d0,3)
  CALL xyz(grd,0.d0,-1.d0,0.d0,4)
  CALL xyz(grd,1.d0,0.d0,0.d0,5)
  CALL xyz(grd,-1.d0,0.d0,0.d0,6)
!              the a2 nodes-these are independent of the quadrature
  IF (grd%nleb == 11.OR.grd%nleb == 13.OR.grd%nleb == 23) THEN
     CALL xyz(grd,sqrt_half,sqrt_half,0.d0,7)
     CALL xyz(grd,sqrt_half,-sqrt_half,0.d0,8)
     CALL xyz(grd,-sqrt_half,sqrt_half,0.d0,9)
     CALL xyz(grd,-sqrt_half,-sqrt_half,0.d0,10)
     CALL xyz(grd,sqrt_half,0.d0,sqrt_half,11)
     CALL xyz(grd,sqrt_half,0.d0,-sqrt_half,12)
     CALL xyz(grd,-sqrt_half,0.d0,sqrt_half,13)
     CALL xyz(grd,-sqrt_half,0.d0,-sqrt_half,14)
     CALL xyz(grd,0.d0,sqrt_half,sqrt_half,15)
     CALL xyz(grd,0.d0,sqrt_half,-sqrt_half,16)
     CALL xyz(grd,0.d0,-sqrt_half,sqrt_half,17)
     CALL xyz(grd,0.d0,-sqrt_half,-sqrt_half,18)
     ii=18
  ELSE
     ii=6
  END IF
!              the a3 nodes-these are independent of the quadrature
  CALL xyz(grd,sqrt_third,sqrt_third,sqrt_third,ii+1)
  CALL xyz(grd,sqrt_third,sqrt_third,-sqrt_third,ii+2)
  CALL xyz(grd,sqrt_third,-sqrt_third,sqrt_third,ii+3)
  CALL xyz(grd,sqrt_third,-sqrt_third,-sqrt_third,ii+4)
  CALL xyz(grd,-sqrt_third,sqrt_third,sqrt_third,ii+5)
  CALL xyz(grd,-sqrt_third,sqrt_third,-sqrt_third,ii+6)
  CALL xyz(grd,-sqrt_third,-sqrt_third,sqrt_third,ii+7)
  CALL xyz(grd,-sqrt_third,-sqrt_third,-sqrt_third,ii+8)
  DO i=1,4
     m(i)=0.d0
     l(i)=0.d0
  END DO
  IF (grd%nleb == 9) THEN
      biga(1)=.00952380952387D0
      biga(3)=.0321428571429D0
      bigc(1)=.0285714285714D0
      p(1)=.888073833977D0
      q(1)=.459700843381D0
      DO i=1,6
         grd%wtleb(i)=biga(1)
      END DO
      DO i=7,14
         grd%wtleb(i)=biga(3)
      END DO
      CALL getcnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),p,q,bigc)
  ELSE IF (grd%nleb == 11) THEN
      m(1)=.90453403373329086794D0
      l(1)=.30151134457800000000D0
      bigc(1)=0.d0
      biga(1)=.01269841269841269841D0
      biga(2)=.02257495590828924162D0
      biga(3)=.02109375000000000000D0
      bigb(1)=.02017333553791887125D0
      bigb(2)=0.d0
      bigb(3)=0.d0
      bigb(4)=0.d0
      DO i=1,6
         grd%wtleb(i)=biga(1)
         grd%wtleb(i+6)=biga(2)
         grd%wtleb(i+12)=biga(2)
         grd%wtleb(i+18)=biga(3)
      END DO
      grd%wtleb(25)=biga(3)
      grd%wtleb(26)=biga(3)
      CALL getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
  ELSE IF (grd%nleb == 13) THEN
           m(1)=.73379938570534280703D0
           l(1)=.48038446141526140045D0
           bigc(1)=.01652217099371570916D0
           biga(1)=.00051306717973384640D0
           biga(2)=.01660406956574203960D0
           biga(3)=-.02958603896103896103D0
           bigb(1)=.02657620708215946311D0
           p(1)=.94715622136258788774D0
           q(1)=.32077264898077643273D0
           DO i=1,6
              grd%wtleb(i)=biga(1)
              grd%wtleb(i+6)=biga(2)
              grd%wtleb(i+12)=biga(2)
              grd%wtleb(i+18)=biga(3)
          END DO
          grd%wtleb(25)=biga(3)
          grd%wtleb(26)=biga(3)
          CALL getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
          CALL getcnd(grd%angleb(1,ii+33),grd%wtleb(ii+33),p,q,bigc)
  ELSE IF (grd%nleb == 17) THEN
           n1=3
           bigc(1)=.00969499636166D0
           biga(1)=.00382827049494D0
           biga(2)=0.d0
           biga(3)=.00979373751249D0
           bigb(1)=.00821173728319D0
           bigb(2)=.00959547133607D0
           bigb(3)=.00994281489118D0
           bigb(4)=0.d0
           m(1)=.965124035087D0
           l(1)=.185115635345D0
           m(2)=.828769981253D0
           l(2)=.395689473056D0
           m(3)=.215957291846D0
           l(3)=.690421048382D0
           p(1)=.878158910604D0
           q(1)=.478369028812D0
           DO i=1,6
              grd%wtleb(i)=biga(1)
           END DO
           DO i=7,14
              grd%wtleb(i)=biga(3)
           END DO
           CALL getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
           CALL getcnd(grd%angleb(1,ii+81),grd%wtleb(ii+81),p,q,bigc)
  ELSE IF (grd%nleb == 23) THEN
           n1=4
           n2=1
           n3=1
           bigc(1)=5.05184606462D-03
           p(1)=.938319218138D0
           q(1)=.345770219761D0
           bigd(1)=5.53024891623D-03
           r(1)=.836036015482D0
           u(1)=.159041710538D0
           w(1)=.525118572443D0
           biga(1)=1.78234044724D-03
           biga(2)=5.71690594998D-03
           biga(3)=5.57338317884D-03
           bigb(1)=5.51877146727D-03
           bigb(2)=5.15823771181D-03
           bigb(3)=5.60870408259D-03
           bigb(4)=4.10677702817D-03
           m(1)=.777493219315D0
           l(1)=.444693317871D0
           m(2)=.912509096867D0
           l(2)=.289246562758D0
           m(3)=.314196994183D0
           l(3)=.671297344270D0
           m(4)=.982972302707D0
           l(4)=.129933544765D0
           DO i=1,6
              grd%wtleb(i)=biga(1)
           END DO
           DO i=7,18
              grd%wtleb(i)=biga(2)
           END DO
           DO i=19,26
              grd%wtleb(i)=biga(3)
           END DO
           CALL getbnd(grd%angleb(1,ii+9),grd%wtleb(ii+9),l,m,bigb,grd%nleb)
           CALL getcnd(grd%angleb(1,ii+105),grd%wtleb(ii+105),p,q,bigc)
           CALL getdnd(grd%angleb(1,ii+129),grd%wtleb(ii+129),r,u,w,bigd)
  END IF
  sum=0.d0
  DO i=1,grd%nang
     grd%wtleb(i)=grd%wtleb(i)*fourpi
     sum=sum+grd%wtleb(i)
  END DO
  sum=sum/fourpi
  WRITE(iout,1) sum
1 FORMAT (/,' sum of lebedev weights = ',e15.8)
  CALL ang(grd%angleb,grd%cthet,grd%sthet,grd%cphi,grd%sphi,grd%phpt,grd%nang)
  DO i=1,npts
     angle(i,1)=grd%cthet(i)
     angle(i,2)=grd%phpt(i)
  END DO
  CALL iosys ('write real "lebedev angular points '//str//'" to grid',2*npts,angle,0,' ')
CALL iosys ('write real "lebedev angular weights '//str//'" to grid',npts,grd%wtleb,0,' ')
!********************************************************************************
!********************************************************************************
  END SUBROUTINE Lebdev
!********************************************************************************
!********************************************************************************
!deck aij.f
!***begin prologue     aij
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           size, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            calculate atomic sizes based on bragg-slater radii
!***references         becke papers in jcp on lda and poisson equation.

!***routines called
!***end prologue       size
  FUNCTION aij (za,zb)
  IMPLICIT NONE
  REAL(idp)                       :: za
  REAL(idp)                       :: zb
  REAL(idp)                       :: aij
  REAL(idp), DIMENSION(10)        :: r
  REAL(idp)                       :: chi
  REAL(idp)                       :: uab
  REAL(idp)                       :: na
  REAL(idp)                       :: nb
  DATA r /10*1.d0 /
  na=za
  nb=zb
  chi=r(na)/r(nb)
  uab=(chi-1.d0)/(chi+1.d0)
  aij=uab/(uab*uab-1.d0)
  IF(ABS(aij) > .5D0) THEN
     aij=SIGN(.5D0,aij)
  END IF
!********************************************************************************
  END FUNCTION aij
!********************************************************************************
!********************************************************************************
!deck getbnd.f
!***begin prologue     getbnd
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           size, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            data for b nodes in lebedev quadrature
!***references         lebedev paper

!***routines called
!***end prologue       size
  SUBROUTINE getbnd (pt,wt,l,m,bigb,nleb)
  SUBROUTINE getbnd (grd,l,m,bigb,nleb)
  IMPLICIT NONE
  REAL(idp), DIMENSION(3,:)                :: pt
  REAL(idp), DIMENSION(:)                  :: wt
  REAL(idp), DIMENSION(:)                  :: l
  REAL(idp), DIMENSION(:)                  :: m
  REAL(idp), DIMENSION(:)                  :: bigb
  INTEGER                                  :: nleb
!
!          eight (l,l,m)
  ii=0
  CALL xyzw(pt(1,ii+1),wt(ii+1),l(1),l(1),m(1),bigb(1))
  CALL xyzw(pt(1,ii+2),wt(ii+2),l(1),l(1),-m(1),bigb(1))
  CALL xyzw(pt(1,ii+3),wt(ii+3),l(1),-l(1),m(1),bigb(1))
  CALL xyzw(pt(1,ii+4),wt(ii+4),l(1),-l(1),-m(1),bigb(1))
  CALL xyzw(pt(1,ii+5),wt(ii+5),-l(1),l(1),m(1),bigb(1))
  CALL xyzw(pt(1,ii+6),wt(ii+6),-l(1),l(1),-m(1),bigb(1))
  CALL xyzw(pt(1,ii+7),wt(ii+7),-l(1),-l(1),m(1),bigb(1))
  CALL xyzw(pt(1,ii+8),wt(ii+8),-l(1),-l(1),-m(1),bigb(1))
!
!         eight (l,m,l)
  ii=8
  CALL xyzw(pt(1,ii+1),wt(ii+1),l(1),m(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+2),wt(ii+2),l(1),m(1),-l(1),bigb(1))
  CALL xyzw(pt(1,ii+3),wt(ii+3),l(1),-m(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+4),wt(ii+4),l(1),-m(1),-l(1),bigb(1))
  CALL xyzw(pt(1,ii+5),wt(ii+5),-l(1),m(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+6),wt(ii+6),-l(1),m(1),-l(1),bigb(1))
  CALL xyzw(pt(1,ii+7),wt(ii+7),-l(1),-m(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+8),wt(ii+8),-l(1),-m(1),-l(1),bigb(1))
!
!         eight (m,l,l)
ii=16
  CALL xyzw(pt(1,ii+1),wt(ii+1),m(1),l(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+2),wt(ii+2),m(1),l(1),-l(1),bigb(1))
  CALL xyzw(pt(1,ii+3),wt(ii+3),m(1),-l(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+4),wt(ii+4),m(1),-l(1),-l(1),bigb(1))
  CALL xyzw(pt(1,ii+5),wt(ii+5),-m(1),l(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+6),wt(ii+6),-m(1),l(1),-l(1),bigb(1))
  CALL xyzw(pt(1,ii+7),wt(ii+7),-m(1),-l(1),l(1),bigb(1))
  CALL xyzw(pt(1,ii+8),wt(ii+8),-m(1),-l(1),-l(1),bigb(1))
  IF (nleb == 13) THEN
      RETURN
  END IF
  IF(nleb > 11) THEN
!
!         eight (l,l,m)
!
     ii=24
     CALL xyzw(pt(1,ii+1),wt(ii+1),l(2),l(2),m(2),bigb(2))
     CALL xyzw(pt(1,ii+2),wt(ii+2),l(2),l(2),-m(2),bigb(2))
     CALL xyzw(pt(1,ii+3),wt(ii+3),l(2),-l(2),m(2),bigb(2))
     CALL xyzw(pt(1,ii+4),wt(ii+4),l(2),-l(2),-m(2),bigb(2))
     CALL xyzw(pt(1,ii+5),wt(ii+5),-l(2),l(2),m(2),bigb(2))
     CALL xyzw(pt(1,ii+6),wt(ii+6),-l(2),l(2),-m(2),bigb(2))
     CALL xyzw(pt(1,ii+7),wt(ii+7),-l(2),-l(2),m(2),bigb(2))
     CALL xyzw(pt(1,ii+8),wt(ii+8),-l(2),-l(2),-m(2),bigb(2))
!
!         eight (l,m,l)
!
     ii=32
     CALL xyzw(pt(1,ii+1),wt(ii+1),l(2),m(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+2),wt(ii+2),l(2),m(2),-l(2),bigb(2))
     CALL xyzw(pt(1,ii+3),wt(ii+3),l(2),-m(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+4),wt(ii+4),l(2),-m(2),-l(2),bigb(2))
     CALL xyzw(pt(1,ii+5),wt(ii+5),-l(2),m(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+6),wt(ii+6),-l(2),m(2),-l(2),bigb(2))
     CALL xyzw(pt(1,ii+7),wt(ii+7),-l(2),-m(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+8),wt(ii+8),-l(2),-m(2),-l(2),bigb(2))
!
!         eight (m,l,l)
!
     ii=40
     CALL xyzw(pt(1,ii+1),wt(ii+1),m(2),l(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+2),wt(ii+2),m(2),l(2),-l(2),bigb(2))
     CALL xyzw(pt(1,ii+3),wt(ii+3),m(2),-l(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+4),wt(ii+4),m(2),-l(2),-l(2),bigb(2))
     CALL xyzw(pt(1,ii+5),wt(ii+5),-m(2),l(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+6),wt(ii+6),-m(2),l(2),-l(2),bigb(2))
     CALL xyzw(pt(1,ii+7),wt(ii+7),-m(2),-l(2),l(2),bigb(2))
     CALL xyzw(pt(1,ii+8),wt(ii+8),-m(2),-l(2),-l(2),bigb(2))
!
!         eight (l,l,m)
!
     ii=48
     CALL xyzw(pt(1,ii+1),wt(ii+1),l(3),l(3),m(3),bigb(3))
     CALL xyzw(pt(1,ii+2),wt(ii+2),l(3),l(3),-m(3),bigb(3))
     CALL xyzw(pt(1,ii+3),wt(ii+3),l(3),-l(3),m(3),bigb(3))
     CALL xyzw(pt(1,ii+4),wt(ii+4),l(3),-l(3),-m(3),bigb(3))
     CALL xyzw(pt(1,ii+5),wt(ii+5),-l(3),l(3),m(3),bigb(3))
     CALL xyzw(pt(1,ii+6),wt(ii+6),-l(3),l(3),-m(3),bigb(3))
     CALL xyzw(pt(1,ii+7),wt(ii+7),-l(3),-l(3),m(3),bigb(3))
     CALL xyzw(pt(1,ii+8),wt(ii+8),-l(3),-l(3),-m(3),bigb(3))
!
!         eight (l,m,l)
!
     ii=56
     CALL xyzw(pt(1,ii+1),wt(ii+1),l(3),m(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+2),wt(ii+2),l(3),m(3),-l(3),bigb(3))
     CALL xyzw(pt(1,ii+3),wt(ii+3),l(3),-m(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+4),wt(ii+4),l(3),-m(3),-l(3),bigb(3))
     CALL xyzw(pt(1,ii+5),wt(ii+5),-l(3),m(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+6),wt(ii+6),-l(3),m(3),-l(3),bigb(3))
     CALL xyzw(pt(1,ii+7),wt(ii+7),-l(3),-m(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+8),wt(ii+8),-l(3),-m(3),-l(3),bigb(3))
!
!         eight (m,l,l)
!
     ii=64
     CALL xyzw(pt(1,ii+1),wt(ii+1),m(3),l(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+2),wt(ii+2),m(3),l(3),-l(3),bigb(3))
     CALL xyzw(pt(1,ii+3),wt(ii+3),m(3),-l(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+4),wt(ii+4),m(3),-l(3),-l(3),bigb(3))
     CALL xyzw(pt(1,ii+5),wt(ii+5),-m(3),l(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+6),wt(ii+6),-m(3),l(3),-l(3),bigb(3))
     CALL xyzw(pt(1,ii+7),wt(ii+7),-m(3),-l(3),l(3),bigb(3))
     CALL xyzw(pt(1,ii+8),wt(ii+8),-m(3),-l(3),-l(3),bigb(3))
  END IF
!
  IF (nleb > 17) THEN
!
!         eight (l,l,m)
      ii=72
      CALL xyzw(pt(1,ii+1),wt(ii+1),l(4),l(4),m(4),bigb(4))
      CALL xyzw(pt(1,ii+2),wt(ii+2),l(4),l(4),-m(4),bigb(4))
      CALL xyzw(pt(1,ii+3),wt(ii+3),l(4),-l(4),m(4),bigb(4))
      CALL xyzw(pt(1,ii+4),wt(ii+4),l(4),-l(4),-m(4),bigb(4))
      CALL xyzw(pt(1,ii+5),wt(ii+5),-l(4),l(4),m(4),bigb(4))
      CALL xyzw(pt(1,ii+6),wt(ii+6),-l(4),l(4),-m(4),bigb(4))
      CALL xyzw(pt(1,ii+7),wt(ii+7),-l(4),-l(4),m(4),bigb(4))
      CALL xyzw(pt(1,ii+8),wt(ii+8),-l(4),-l(4),-m(4),bigb(4))
!
!         eight (l,m,l)
!
      ii=80
      CALL xyzw(pt(1,ii+1),wt(ii+1),l(4),m(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+2),wt(ii+2),l(4),m(4),-l(4),bigb(4))
      CALL xyzw(pt(1,ii+3),wt(ii+3),l(4),-m(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+4),wt(ii+4),l(4),-m(4),-l(4),bigb(4))
      CALL xyzw(pt(1,ii+5),wt(ii+5),-l(4),m(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+6),wt(ii+6),-l(4),m(4),-l(4),bigb(4))
      CALL xyzw(pt(1,ii+7),wt(ii+7),-l(4),-m(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+8),wt(ii+8),-l(4),-m(4),-l(4),bigb(4))
!
!         eight (m,l,l)
!
      ii=88
      CALL xyzw(pt(1,ii+1),wt(ii+1),m(4),l(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+2),wt(ii+2),m(4),l(4),-l(4),bigb(4))
      CALL xyzw(pt(1,ii+3),wt(ii+3),m(4),-l(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+4),wt(ii+4),m(4),-l(4),-l(4),bigb(4))
      CALL xyzw(pt(1,ii+5),wt(ii+5),-m(4),l(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+6),wt(ii+6),-m(4),l(4),-l(4),bigb(4))
      CALL xyzw(pt(1,ii+7),wt(ii+7),-m(4),-l(4),l(4),bigb(4))
      CALL xyzw(pt(1,ii+8),wt(ii+8),-m(4),-l(4),-l(4),bigb(4))
  END IF
!
END SUBROUTINE getbnd
!********************************************************************************
!********************************************************************************
!deck getcnd.f
!***begin prologue     getcnd
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           size, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            data for c nodes in lebedev quadrature
!***references         lebedev paper

!***routines called
!***end prologue       size

SUBROUTINE getcnd (pt,wt,p,q,bigc)

REAL*8, INTENT(IN OUT)                   :: pt(3,*)
REAL*8, INTENT(IN OUT)                   :: wt(*)
REAL*8, INTENT(IN OUT)                   :: p(1)
REAL*8, INTENT(IN OUT)                   :: q(1)
REAL*8, INTENT(IN OUT)                   :: bigc(1)
IMPLICIT INTEGER (a-z)


COMMON /io/ inp, iout

CALL xyzw(pt(1,1),wt(1),p(1),q(1),0.d0,bigc(1))
CALL xyzw(pt(1,2),wt(2),p(1),-q(1),0.d0,bigc(1))
CALL xyzw(pt(1,3),wt(3),-p(1),q(1),0.d0,bigc(1))
CALL xyzw(pt(1,4),wt(4),-p(1),-q(1),0.d0,bigc(1))
CALL xyzw(pt(1,5),wt(5),p(1),0.d0,q(1),bigc(1))
CALL xyzw(pt(1,6),wt(6),p(1),0.d0,-q(1),bigc(1))
CALL xyzw(pt(1,7),wt(7),-p(1),0.d0,q(1),bigc(1))
CALL xyzw(pt(1,8),wt(8),-p(1),0.d0,-q(1),bigc(1))
CALL xyzw(pt(1,9),wt(9),0.d0,p(1),q(1))
CALL xyzw(pt(1,10),wt(10),0.d0,p(1),-q(1))
CALL xyzw(pt(1,11),wt(11),0.d0,-p(1),q(1))
CALL xyzw(pt(1,12),wt(12),0.d0,-p(1),-q(1))
CALL xyzw(pt(1,13),wt(13),q(1),p(1),0.d0,bigc(1))
CALL xyzw(pt(1,14),wt(14),q(1),-p(1),0.d0,bigc(1))
CALL xyzw(pt(1,15),wt(15),-q(1),p(1),0.d0,bigc(1))
CALL xyzw(pt(1,16),wt(16),-q(1),-p(1),0.d0,bigc(1))
CALL xyzw(pt(1,17),wt(17),q(1),0.d0,p(1),bigc(1))
CALL xyzw(pt(1,18),wt(18),q(1),0.d0,-p(1),bigc(1))
CALL xyzw(pt(1,19),wt(19),-q(1),0.d0,p(1),bigc(1))
CALL xyzw(pt(1,20),wt(20),-q(1),0.d0,-p(1),bigc(1))
CALL xyzw(pt(1,21),wt(21),0.d0,q(1),p(1),bigc(1))
CALL xyzw(pt(1,22),wt(22),0.d0,q(1),-p(1),bigc(1))
CALL xyzw(pt(1,23),wt(23),0.d0,-q(1),p(1),bigc(1))
CALL xyzw(pt(1,24),wt(24),0.d0,-q(1),-p(1),bigc(1))
RETURN
END SUBROUTINE getcnd
!********************************************************************************
!********************************************************************************
!deck getdnd.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2014-06-18  Time: 08:28:30
 
!***begin prologue     getdnd
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           size, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            data for d nodes in lebedev quadrature
!***references         lebedev paper

!***routines called
!***end prologue       size

SUBROUTINE getdnd (pt,wt,r,u,w,bigd)

REAL*8, INTENT(IN OUT)                   :: pt(3,*)
REAL*8, INTENT(IN OUT)                   :: wt(*)
REAL*8, INTENT(IN OUT)                   :: r(1)
REAL*8, INTENT(IN OUT)                   :: u(1)
REAL*8, INTENT(IN OUT)                   :: w(1)
REAL*8, INTENT(IN OUT)                   :: bigd(1)
IMPLICIT INTEGER (a-z)


COMMON /io/ inp, iout

CALL xyzw(pt(1,1),wt(1),r(1),u(1),w(1),bigd(1))
CALL xyzw(pt(1,2),wt(2),r(1),-u(1),w(1),bigd(1))
CALL xyzw(pt(1,3),wt(3),-r(1),u(1),w(1),bigd(1))
CALL xyzw(pt(1,4),wt(4),-r(1),-u(1),w(1),bigd(1))
CALL xyzw(pt(1,5),wt(5),r(1),u(1),-w(1),bigd(1))
CALL xyzw(pt(1,6),wt(6),r(1),-u(1),-w(1),bigd(1))
CALL xyzw(pt(1,7),wt(7),-r(1),u(1),-w(1),bigd(1))
CALL xyzw(pt(1,8),wt(8),-r(1),-u(1),-w(1),bigd(1))
CALL xyzw(pt(1,9),wt(9),r(1),w(1),u(1),bigd(1))
CALL xyzw(pt(1,10),wt(10),r(1),-w(1),u(1),bigd(1))
CALL xyzw(pt(1,11),wt(11),-r(1),w(1),u(1),bigd(1))
CALL xyzw(pt(1,12),wt(12),-r(1),-w(1),u(1),bigd(1))
CALL xyzw(pt(1,13),wt(13),r(1),w(1),-u(1),bigd(1))
CALL xyzw(pt(1,14),wt(14),r(1),-w(1),-u(1),bigd(1))
CALL xyzw(pt(1,15),wt(15),-r(1),w(1),-u(1),bigd(1))
CALL xyzw(pt(1,16),wt(16),-r(1),-w(1),-u(1),bigd(1))
CALL xyzw(pt(1,17),wt(17),u(1),r(1),w(1),bigd(1))
CALL xyzw(pt(1,18),wt(18),u(1),-r(1),w(1),bigd(1))
CALL xyzw(pt(1,19),wt(19),-u(1),r(1),w(1),bigd(1))
CALL xyzw(pt(1,20),wt(20),-u(1),-r(1),w(1),bigd(1))
CALL xyzw(pt(1,21),wt(21),u(1),r(1),-w(1),bigd(1))
CALL xyzw(pt(1,22),wt(22),u(1),-r(1),-w(1),bigd(1))
CALL xyzw(pt(1,23),wt(23),-u(1),r(1),-w(1),bigd(1))
CALL xyzw(pt(1,24),wt(24),-u(1),-r(1),-w(1),bigd(1))
CALL xyzw(pt(1,25),wt(25),u(1),w(1),r(1),bigd(1))
CALL xyzw(pt(1,26),wt(26),u(1),-w(1),r(1),bigd(1))
CALL xyzw(pt(1,27),wt(27),-u(1),w(1),r(1),bigd(1))
CALL xyzw(pt(1,28),wt(28),-u(1),-w(1),r(1),bigd(1))
CALL xyzw(pt(1,29),wt(29),u(1),w(1),-r(1),bigd(1))
CALL xyzw(pt(1,30),wt(30),u(1),-w(1),-r(1),bigd(1))
CALL xyzw(pt(1,31),wt(31),-u(1),w(1),-r(1),bigd(1))
CALL xyzw(pt(1,32),wt(32),-u(1),-w(1),-r(1),bigd(1))
CALL xyzw(pt(1,33),wt(33),w(1),u(1),r(1),bigd(1))
CALL xyzw(pt(1,34),wt(34),w(1),-u(1),r(1),bigd(1))
CALL xyzw(pt(1,35),wt(35),-w(1),u(1),r(1),bigd(1))
CALL xyzw(pt(1,36),wt(36),-w(1),-u(1),r(1),bigd(1))
CALL xyzw(pt(1,37),wt(37),w(1),u(1),-r(1),bigd(1))
CALL xyzw(pt(1,38),wt(38),w(1),-u(1),-r(1),bigd(1))
CALL xyzw(pt(1,39),wt(39),-w(1),u(1),-r(1),bigd(1))
CALL xyzw(pt(1,40),wt(40),-w(1),-u(1),-r(1),bigd(1))
CALL xyzw(pt(1,41),wt(41),w(1),r(1),u(1),bigd(1))
CALL xyzw(pt(1,42),wt(42),w(1),-r(1),u(1),bigd(1))
CALL xyzw(pt(1,43),wt(43),-w(1),r(1),u(1),bigd(1))
CALL xyzw(pt(1,44),wt(44),-w(1),-r(1),u(1),bigd(1))
CALL xyzw(pt(1,45),wt(45),w(1),r(1),-u(1),bigd(1))
CALL xyzw(pt(1,46),wt(46),w(1),-r(1),-u(1),bigd(1))
CALL xyzw(pt(1,47),wt(47),-w(1),r(1),-u(1),bigd(1))
CALL xyzw(pt(1,48),wt(48),-w(1),-r(1),-u(1),bigd(1))
RETURN
END SUBROUTINE getdnd
!********************************************************************************
!********************************************************************************
!deck xyz.f
!***begin prologue     xyz
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           xyz, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            fill a vector
!***references

!***routines called
!***end prologue       xyz
  SUBROUTINE xyz (grd,x,y,z,i)
  IMPLICIT NONE
  Type(GRID)                               :: grd
  REAL(idp)                                :: x
  REAL(idp)                                :: y
  REAL(idp)                                :: z
  INTEGER                                  :: i
  REAL(idp), DIMENSION(3,18), PARAMETER    ::v
  grd%angleb(:,i)=v(1)
!********************************************************************************
  END SUBROUTINE xyz
!********************************************************************************
!deck xyzw.f
!***begin prologue     xyzw
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           xyzw, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            fill a vector and weight
!*o**references

!***routines called
!***end prologue       xyz
!
  SUBROUTINE xyzw (pt,wt,x,y,z,wtval)
  IMPLICIT NONE
  REAL(idp)                      :: pt(3)
  REAL(idp)                      :: wt
  REAL(idp)                      :: x
  REAL(idp)                      :: y
  REAL(idp)                      :: z
  REAL(idp)                      :: wtval
  pt(1)=x
  pt(2)=y
  pt(3)=z
  wt=wtval
!********************************************************************************
!********************************************************************************
  END SUBROUTINE xyzw
!********************************************************************************
!********************************************************************************
  END MODULE Points_Weights
!********************************************************************************
!********************************************************************************


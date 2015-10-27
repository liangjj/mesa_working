*deck speli4
      subroutine speli4 (iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta,
     +   c, d, n, nbdcnd, bdc, bdd, cofx, an, bn, cn, dn, un, zn, am,
     +   bm, cm, dm, um, zm, grhs, usol, idmn, w, pertrb, ierror)
c***begin prologue  speli4
c***subsidiary
c***purpose  subsidiary to sepx4
c***library   slatec
c***type      single precision (speli4-s)
c***author  (unknown)
c***description
c
c     speli4 sets up vectors and arrays for input to blktri
c     and computes a second order solution in usol.  a return jump to
c     sepx4 occurs if iorder=2.  if iorder=4 a fourth order
c     solution is generated in usol.
c
c***see also  sepx4
c***routines called  chksn4, defe4, genbun, minso4, ortho4, tris4
c***common blocks    spl4
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  speli4
c
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                w(*)
      dimension       grhs(idmn,*)           ,usol(idmn,*)
      dimension       an(*)      ,bn(*)      ,cn(*)      ,dn(*)      ,
     1                un(*)      ,zn(*)
      dimension       am(*)      ,bm(*)      ,cm(*)      ,dm(*)      ,
     1                um(*)      ,zm(*)
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     1                ait        ,bit        ,cit        ,dit        ,
     2                mit        ,nit        ,is         ,ms         ,
     3                js         ,ns         ,dlx        ,dly        ,
     4                tdlx3      ,tdly3      ,dlx4       ,dly4
      logical         singlr
      external cofx
c***first executable statement  speli4
      kswx = mbdcnd+1
      kswy = nbdcnd+1
      k = m+1
      l = n+1
      ait = a
      bit = b
      cit = c
      dit = d
      dly=(dit-cit)/n
c
c     set right hand side values from grhs in usol on the interior
c     and non-specified boundaries.
c
      do  20 i=2,m
         do  10 j=2,n
      usol(i,j)=dly**2*grhs(i,j)
   10    continue
   20 continue
      if (kswx.eq.2 .or. kswx.eq.3) go to  40
      do  30 j=2,n
      usol(1,j)=dly**2*grhs(1,j)
   30 continue
   40 continue
      if (kswx.eq.2 .or. kswx.eq.5) go to  60
      do  50 j=2,n
      usol(k,j)=dly**2*grhs(k,j)
   50 continue
   60 continue
      if (kswy.eq.2 .or. kswy.eq.3) go to  80
      do  70 i=2,m
      usol(i,1)=dly**2*grhs(i,1)
   70 continue
   80 continue
      if (kswy.eq.2 .or. kswy.eq.5) go to 100
      do  90 i=2,m
      usol(i,l)=dly**2*grhs(i,l)
   90 continue
  100 continue
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.3)
     1usol(1,1)=dly**2*grhs(1,1)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.3)
     1usol(k,1)=dly**2*grhs(k,1)
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.5)
     1usol(1,l)=dly**2*grhs(1,l)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.5)
     1usol(k,l)=dly**2*grhs(k,l)
c
c     set switches for periodic or non-periodic boundaries
c
      mp=1
      if(kswx.eq.1) mp=0
      np=nbdcnd
c
c     set dlx,dly and size of block tri-diagonal system generated
c     in nint,mint
c
      dlx = (bit-ait)/m
      mit = k-1
      if (kswx .eq. 2) mit = k-2
      if (kswx .eq. 4) mit = k
      dly = (dit-cit)/n
      nit = l-1
      if (kswy .eq. 2) nit = l-2
      if (kswy .eq. 4) nit = l
      tdlx3 = 2.0*dlx**3
      dlx4 = dlx**4
      tdly3 = 2.0*dly**3
      dly4 = dly**4
c
c     set subscript limits for portion of array to input to blktri
c
      is = 1
      js = 1
      if (kswx.eq.2 .or. kswx.eq.3) is = 2
      if (kswy.eq.2 .or. kswy.eq.3) js = 2
      ns = nit+js-1
      ms = mit+is-1
c
c     set x - direction
c
      do 110 i=1,mit
         xi = ait+(is+i-2)*dlx
         call cofx (xi,ai,bi,ci)
         axi = (ai/dlx-0.5*bi)/dlx
         bxi = -2.*ai/dlx**2+ci
         cxi = (ai/dlx+0.5*bi)/dlx
      am(i)=dly**2*axi
      bm(i)=dly**2*bxi
      cm(i)=dly**2*cxi
  110 continue
c
c     set y direction
c
      do 120 j=1,nit
      dyj=1.0
      eyj=-2.0
      fyj=1.0
         an(j) = dyj
         bn(j) = eyj
         cn(j) = fyj
  120 continue
c
c     adjust edges in x direction unless periodic
c
      ax1 = am(1)
      cxm = cm(mit)
      go to (170,130,150,160,140),kswx
c
c     dirichlet-dirichlet in x direction
c
  130 am(1) = 0.0
      cm(mit) = 0.0
      go to 170
c
c     mixed-dirichlet in x direction
c
  140 am(1) = 0.0
      bm(1) = bm(1)+2.*alpha*dlx*ax1
      cm(1) = cm(1)+ax1
      cm(mit) = 0.0
      go to 170
c
c     dirichlet-mixed in x direction
c
  150 am(1) = 0.0
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*beta*dlx*cxm
      cm(mit) = 0.0
      go to 170
c
c     mixed - mixed in x direction
c
  160 continue
      am(1) = 0.0
      bm(1) = bm(1)+2.*dlx*alpha*ax1
      cm(1) = cm(1)+ax1
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*dlx*beta*cxm
      cm(mit) = 0.0
  170 continue
c
c     adjust in y direction unless periodic
c
      dy1 = an(1)
      fyn = cn(nit)
      gama=0.0
      xnu=0.0
      go to (220,180,200,210,190),kswy
c
c     dirichlet-dirichlet in y direction
c
  180 continue
      an(1) = 0.0
      cn(nit) = 0.0
      go to 220
c
c     mixed-dirichlet in y direction
c
  190 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      cn(nit) = 0.0
      go to 220
c
c     dirichlet-mixed in y direction
c
  200 an(1) = 0.0
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.*dly*xnu*fyn
      cn(nit) = 0.0
      go to 220
c
c     mixed - mixed direction in y direction
c
  210 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.0*dly*xnu*fyn
      cn(nit) = 0.0
  220 if (kswx .eq. 1) go to 270
c
c     adjust usol along x edge
c
      do 260 j=js,ns
         if (kswx.ne.2 .and. kswx.ne.3) go to 230
         usol(is,j) = usol(is,j)-ax1*usol(1,j)
         go to 240
  230    usol(is,j) = usol(is,j)+2.0*dlx*ax1*bda(j)
  240    if (kswx.ne.2 .and. kswx.ne.5) go to 250
         usol(ms,j) = usol(ms,j)-cxm*usol(k,j)
         go to 260
  250    usol(ms,j) = usol(ms,j)-2.0*dlx*cxm*bdb(j)
  260 continue
  270 if (kswy .eq. 1) go to 320
c
c     adjust usol along y edge
c
      do 310 i=is,ms
         if (kswy.ne.2 .and. kswy.ne.3) go to 280
         usol(i,js) = usol(i,js)-dy1*usol(i,1)
         go to 290
  280    usol(i,js) = usol(i,js)+2.0*dly*dy1*bdc(i)
  290    if (kswy.ne.2 .and. kswy.ne.5) go to 300
         usol(i,ns) = usol(i,ns)-fyn*usol(i,l)
         go to 310
  300    usol(i,ns) = usol(i,ns)-2.0*dly*fyn*bdd(i)
  310 continue
  320 continue
c
c     save adjusted edges in grhs if iorder=4
c
      if (iorder .ne. 4) go to 350
      do 330 j=js,ns
         grhs(is,j) = usol(is,j)
         grhs(ms,j) = usol(ms,j)
  330 continue
      do 340 i=is,ms
         grhs(i,js) = usol(i,js)
         grhs(i,ns) = usol(i,ns)
  340 continue
  350 continue
      iord = iorder
      pertrb = 0.0
c
c     check if operator is singular
c
      call chksn4(mbdcnd,nbdcnd,alpha,beta,cofx,singlr)
c
c     compute non-zero eigenvector in null space of transpose
c     if singular
c
      if (singlr) call tris4 (mit,am,bm,cm,dm,um,zm)
      if (singlr) call tris4 (nit,an,bn,cn,dn,un,zn)
c
c     adjust right hand side if necessary
c
  360 continue
      if (singlr) call ortho4 (usol,idmn,zn,zm,pertrb)
c
c     compute solution
c
c     save adjusted right hand side in grhs
      do 444 j=js,ns
      do 444 i=is,ms
      grhs(i,j)=usol(i,j)
  444 continue
      call genbun(np,nit,mp,mit,am,bm,cm,idmn,usol(is,js),ieror,w)
c     check if error detected in pois
c     this can only correspond to ierror=12
      if(ieror.eq.0) go to 224
c     set error flag if improper coefficients input to pois
      ierror=12
      return
  224 continue
      if (ierror .ne. 0) return
c
c     set periodic boundaries if necessary
c
      if (kswx .ne. 1) go to 380
      do 370 j=1,l
         usol(k,j) = usol(1,j)
  370 continue
  380 if (kswy .ne. 1) go to 400
      do 390 i=1,k
         usol(i,l) = usol(i,1)
  390 continue
  400 continue
c
c     minimize solution with respect to weighted least squares
c     norm if operator is singular
c
      if (singlr) call minso4 (usol,idmn,zn,zm,prtrb)
c
c     return if deferred corrections and a fourth order solution are
c     not flagged
c
      if (iord .eq. 2) return
      iord = 2
c
c     compute new right hand side for fourth order solution
c
      call defe4(cofx,idmn,usol,grhs)
      go to 360
      end

*deck @(#)zprint.f	5.1  11/6/94
      subroutine zprint(nz,ianz,iz,bl,alpha,beta,lbl,lalpha,lbeta,
     $                  bastyp,vname)
c***begin prologue     zprint
c***date written       850601  yymmdd
c***revision date      860701  (yymmdd)
c           1 july 1986  modified by pws at lanl
c            fixed format number for printing of dummy atom numbered 4
c             or more.
c
c***keywords           z-matrix
c***author             martin, richard (lanl)
c***source
c***purpose            prints the z-matrix.
c***description
c     call zprint(nz,ianz,iz,bl,alpha,beta,lbl,lalpha,lbeta,bastyp,vname)
c***references
c***routines called    fillel(util), zconv(l202)
c***end prologue       zprint
      implicit integer(a-z)
      real*8 bl(nz),alpha(nz),beta(nz)
      integer ianz(nz),iz(4,nz),lbl(nz),lalpha(nz),lbeta(nz)
      character el(106)*2
      character*(*) bastyp(nz), vname(nz)
      character*8 blnm,alpnm,betanm
      common/io/inp,iout
      data maxel/104/, el(1)/'x'/
      save maxel,el
c
 1010 format(1x,'z-matrix (angstroms and degrees):')
 1020 format(5x,'cd cent  el basis',6x,'n1',5x,'length',10x,'n2',4x,
     $       'alpha',10x,'n3',5x,'beta',11x,'j')
 2110 format(5x,i2,2x,i2,3x,a2,1x,a8)
 2120 format(5x,i2,      7x,a2,1x,a8)
 2210 format(5x,i2,2x,i2,3x,a2,1x,a8,3x,i2,2x,f9.6,1x,a8,1x)
 2220 format(5x,i2,      7x,a2,1x,a8,3x,i2,2x,f9.6,1x,a8,1x)
 2310 format(5x,i2,2x,i2,3x,a2,1x,a8,3x,i2,2x,f9.6,1x,a8,1x,
     $          i2,1x,f8.3,1x,a8,1x)
 2320 format(5x,i2,      7x,a2,1x,a8,3x,i2,2x,f9.6,1x,a8,1x,
     $          i2,1x,f8.3,1x,a8,1x)
 2410 format(5x,i2,2x,i2,3x,a2,1x,a8,3x,i2,2x,f9.6,1x,a8,1x,
     $          i2,1x,f8.3,1x,a8,1x,i2,1x,f8.3,1x,a8,1x,
     $          i2)
 2420 format(5x,i2,      7x,a2,1x,a8,3x,i2,2x,f9.6,1x,a8,1x,
     $          i2,1x,f8.3,1x,a8,1x,i2,1x,f8.3,1x,a8,1x,
     $          i2)
c
c     first converts from internal (bohr/radian) units to external
c     (angstrom/degree) units before printing.
c
c
c     print the heading.
      write(iout,1010)
      write(iout,1020)
c
c     first card of the z-matrix.
      if(nz.ge.1) then
         call fillel(0,maxel,el(2))
         call zconv('toangdeg',nz,bl,alpha,beta)
         icard=1
         idx=ianz(1)+2
         if(ianz(1).ge.0) then
            icent=1
            write(iout,2110) icard,icent,el(idx),bastyp(icard)
         else
            icent=0
            write(iout,2120) icard,el(idx),bastyp(icard)
         endif
      endif
c
c     second card of the z-matrix.
      if(nz.ge.2) then
         icard=2
         blnm=' '
         if(lbl(icard).ne.0) blnm=vname(iabs(lbl(icard)))(1:8)
         idx=ianz(2)+2
         if(ianz(2).ge.0) then
            icent=icent+1
            write(iout,2210) icard,icent,el(idx),bastyp(icard),iz(1,2),
     $                       bl(2),blnm
         else
            write(iout,2220) icard,      el(idx),bastyp(icard),iz(1,2),
     $                       bl(2),blnm
         endif
      endif
c
c     third card.
      if(nz.ge.3) then
         icard=3
         blnm=' '
         if(lbl(icard).ne.0) blnm=vname(iabs(lbl(icard)))(1:8)
         alpnm=' '
         if(lalpha(icard).ne.0) alpnm=vname(iabs(lalpha(icard)))(1:8)
         idx=ianz(3)+2
         if(ianz(3).ge.0) then
            icent=icent+1
            write(iout,2310) icard,icent,el(idx),bastyp(icard),iz(1,3),
     $                       bl(3),blnm,iz(2,3),alpha(3),alpnm
         else
            write(iout,2320) icard,      el(idx),bastyp(icard),iz(1,3),
     $                       bl(3),blnm,iz(2,3),alpha(3),alpnm
         endif
      endif
c
c     cards 4 through nz.
      if(nz.ge.4) then
         do 100 icard=4,nz
            blnm=' '
            if(lbl(icard).ne.0) blnm=vname(iabs(lbl(icard)))(1:8)
            alpnm=' '
            if(lalpha(icard).ne.0) alpnm=vname(iabs(lalpha(icard)))(1:8)
            betanm=' '
            if(lbeta(icard).ne.0) betanm=vname(iabs(lbeta(icard)))(1:8)
            idx=ianz(icard)+2
            if(ianz(icard).ge.0) then
               icent=icent+1
               write(iout,2410) icard,icent,el(idx),bastyp(icard),
     $                          iz(1,icard),bl(icard),blnm,
     $                          iz(2,icard),alpha(icard),alpnm,
     $                          iz(3,icard),beta(icard),betanm,
     $                          iz(4,icard)
            else
               write(iout,2420) icard,      el(idx),bastyp(icard),
     $                          iz(1,icard),bl(icard),blnm,
     $                          iz(2,icard),alpha(icard),alpnm,
     $                          iz(3,icard),beta(icard),betanm,
     $                          iz(4,icard)
            endif
  100    continue
      endif
c
c     convert back to internal units.
      if(nz.ge.1) call zconv('toborrad',nz,bl,alpha,beta)
c
c
      return
      end

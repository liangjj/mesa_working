*deck @(#)ptgrp.f	5.1  11/6/94
      subroutine ptgrp(maxap3,a,b,c,d,ian,atmchg,natoms,dump,
     $     pgrp,trvec,symflg,scr1,scr2,scr3,scr4,npop,nset)
      implicit real*8(a-h,o-z)
c
c     given the coordinates, c, and the atomic numbers or masses,
c     atmchg, of the natoms atoms in a molecule, determine the point
c     group and impose a standard orientation in cartesian space.  the
c     coordinates of the re-oriented molecule are returned in a and
c     the schonflies symbol for the point group is placed in pgrp.
c     b and d are scratch coordinate arrays while dump
c     is a debug print switch.
c
      character*(*) pgrp,symflg
      character*2 itoc
      integer pos
      logical dump
      common/io/inp,iout
      common/tol/   toler,tol2
      integer ian(natoms), npop(1), nset(1)
      real*8 a(maxap3,3), b(maxap3,3), c(natoms,3), d(maxap3,3)
      real*8 atmchg(1), trvec(3), prmom(3), praxes(3,3)
      real*8 t(3,3), v(3), scr1(1), scr2(1),
     $     scr3(1), scr4(1)
      data zero,one,two,four/0.0d+00,1.0d+00,2.0d+00,4.0d+00/
      save zero,one,two,four
c
 1000 format(' ptgrp-- translation vector:',3f12.6)
 1010 format(' ptgrp-- principal moments and axes of charge:',/
     $     '         moments:',3f12.7,/,
     $     '         axes   :',3f12.6/17x,3f12.6/17x,3f12.6)
 1020 format(' ptgrp-- the molecule is linear.')
 1030 format(' ptgrp-- the molecule is not linear.')
 1040 format(' ptgrp-- the molecule is an asymmetric top.')
 1050 format(' ptgrp-- the molecule is a symmetric top.')
 1060 format(' ptgrp-- the molecule is a spherical top.')
c
c     call rtrace(6hptgrp ,1)
      piovr4 = atan(one)
      pi     = four * piovr4
      halfpi = two  * piovr4
c
c     add three dummy atoms to trace the rotations of the molecule.
c
      numatm = natoms + 3
      do 60 iat=1,natoms
         a(iat,1) = c(iat,1)
         a(iat,2) = c(iat,2)
         a(iat,3) = c(iat,3)
 60   continue
      do 70 iat=1,3
         do 70 ixyz=1,3
 70         a(natoms+iat,ixyz) = zero
            a(natoms+1,1) = one
            a(natoms+2,2) = one
            a(natoms+3,3) = one
c
c     all symmetry elements must pass through the molecules charge
c     center.  translate the molecule so that this unique point is
c     at the origin of the fixed cartesian coordinate system.
c
            call qcentr(maxap3,natoms,a,atmchg,trvec)
            trvec(1) = - trvec(1)
            trvec(2) = - trvec(2)
            trvec(3) = - trvec(3)
            if (dump) write(iout,1000) (trvec(i),i=1,3)
            do 80 iat=1,natoms
               a(iat,1) = a(iat,1) + trvec(1)
               a(iat,2) = a(iat,2) + trvec(2)
               a(iat,3) = a(iat,3) + trvec(3)
 80         continue
c
c     calculate the principal moments and axes of charge.
c
            call secmom(maxap3,natoms,a,atmchg,prmom,praxes)
            if (dump) write(iout,1010) (prmom(i),i=1,3),
     $           ((praxes(j,i),i=1,3),j=1,3)
c
c     if the first moment is zero and the other two are equal,
c     the molecule is linear.
c
            if (abs(prmom(1)) .gt. tol2  .or.
     $           abs(prmom(3)-prmom(2)) .gt. tol2) goto 100
c
c     place the molecule on the z axis and distinguish between
c     d*h and c*v.
c
            if(dump) write(iout,1020)
            call putsym(maxap3,a,b,t,praxes(1,1),numatm,3)
            call oraxis(maxap3,a,b,natoms,atmchg,3)
            call reflct(maxap3,a,b,natoms,t,3)
            call equiv(maxap3,a,b,atmchg,natoms,itst)
            if(itst.eq.0) then
               pgrp='c*v'
               goto 2000
            else
               pgrp='d*h'
               goto 2000
            endif
c
c
c     classify the molecule as being a either a spherical top (itop=3),
c     a symmetric top (itop=2), or an asymmetric top(itop=1).  each
c     type of top will be handled seperately.
c
 100        continue
            if(dump) write(iout,1030)
            itop = 0
            tst1 = prmom(2) - prmom(3)
            tst2 = prmom(1) - prmom(3)
            tst3 = prmom(1) - prmom(2)
            if (abs(tst1) .lt. tol2) itop = itop + 1
            if (abs(tst2) .lt. tol2) itop = itop + 1
            if (abs(tst3) .lt. tol2) itop = itop + 1
            if (itop .ne. 3) itop = itop + 1
            goto(110,300,500), itop
c
c     *------------------------*
c     asymmetric top molecules
c     *------------------------*
c
c     these molecules can have no axes of order greater than 2.  thus
c     the possible point groups are:  d2h, d2, c2v, c2h, c2, ci, cs,
c     and c1.
c
c     align the principal axes with the cartesian axes.
c
 110        continue
            if(dump) write(iout,1040)
            call putsym(maxap3,a,b,t,praxes(1,3),numatm,3)
            call oraxis(maxap3,a,b,natoms,atmchg,3)
            call secmom(maxap3,natoms,a,atmchg,prmom,praxes)
            theta = halfpi
            if (abs(praxes(2,2)) .gt. toler)
     $           theta = - atan(praxes(1,2)/praxes(2,2))
            call rotate(maxap3,a,b,numatm,t,3,theta)
            call move(maxap3,b,a,numatm)
            call oraxis(maxap3,a,b,natoms,atmchg,2)
            call orptst(maxap3,a,natoms,ixyz)
            if (ixyz .ne. 0) call orplan(maxap3,a,b,atmchg,numatm,
     $           prmom,praxes,ixyz)
c
c     test z and y for c2.
c     if both are c2 then test for an inversion center.
c     if yes then d2h.
c     if no  then d2.
c     if only z is c2 go to 200.
c     if only y is c2, rotate y to z and go to 200.
c     if neither is c2, test x for c2.
c     if yes then rotate x to z and go to 200.
c     in no  then continue at 150.
c
            call rotate(maxap3,a,b,natoms,t,3,pi)
            call equiv(maxap3,a,b,atmchg,natoms,iztst)
            call rotate(maxap3,a,b,natoms,t,2,pi)
            call equiv(maxap3,a,b,atmchg,natoms,iytst)
            itst = 2*iztst + iytst + 1
            goto(140,130,200,120), itst
c
c     the molecule is either d2 or d2h.
c
 120        call invert(maxap3,a,b,natoms,t)
            call equiv(maxap3,a,b,atmchg,natoms,itst)
            pgrp='d2'
            if(itst.eq.0) then
               goto 2000
            else
               pgrp='d2h'
               call ord2h(maxap3,a,b,natoms,atmchg,ian)
               goto 2000
            endif
c
c     the y axis is c2 but the z axis is not.
c
 130        call rotate(maxap3,a,b,numatm,t,1,halfpi)
            call oraxis(maxap3,b,a,natoms,atmchg,3)
            call orplan(maxap3,b,a,atmchg,numatm,prmom,praxes,3)
            call move(maxap3,b,a,numatm)
            goto 200
c
c     neither y nor z axes are c2.  check x.
c
 140        call rotate(maxap3,a,b,natoms,t,1,pi)
            call equiv(maxap3,a,b,atmchg,natoms,itst)
            if (itst .eq. 0) goto 150
            call rotate(maxap3,a,b,numatm,t,2,halfpi)
            call oraxis(maxap3,b,a,natoms,atmchg,3)
            call orplan(maxap3,b,a,atmchg,numatm,prmom,praxes,3)
            call move(maxap3,b,a,numatm)
            goto 200
c
c     an asymmetric top molecule has no c2 axes.  the remaining
c     possibilities are cs, ci, and c1.  if cs, the symmetry plane
c     is made coincident with the xy plane.
c
 150        call invert(maxap3,a,b,natoms,t)
            call equiv(maxap3,a,b,atmchg,natoms,itst)
            if (itst .eq. 0) goto 160
            pgrp='ci'
            do 155 i = 1, 3
               do 155 ixyz = 1, 3
 155              a(natoms+i,ixyz) = zero
                  a(natoms+1,1) = one
                  a(natoms+2,2) = one
                  a(natoms+3,3) = one
                  goto 2000
c
 160              do 170 i=1,3
 170                 v(i) = zero
                     do 180 ixyz=1,3
                        call reflct(maxap3,a,b,natoms,t,ixyz)
                        call equiv(maxap3,a,b,atmchg,natoms,itst)
                        if (itst .eq. 0) goto 180
                        v(ixyz) = one
                        call putsym(maxap3,a,b,t,v,numatm,3)
                        call oraxis(maxap3,a,b,natoms,atmchg,3)
                        call orplan(maxap3,a,b,atmchg,numatm,prmom,
     $                       praxes,3)
                        goto 190
 180                 continue
c
                     pgrp='c1'
                     goto 2000
c
 190                 pgrp='cs'
                     call orcn(maxap3,a,b,d,d,atmchg,npop,nset,natoms,
     $                    dump,scr2,
     $                    scr1)
                     goto 2000
c
c     an asymmetric top molecule has but one c2 axis; the axis is
c     now coincident with the z axis.  the possible point groups are
c     c2h, c2v, and c2.
c
 200                 continue
                     call reflct(maxap3,a,b,natoms,t,3)
                     call equiv(maxap3,a,b,atmchg,natoms,itst)
                     if (itst .eq. 0) goto 210
                     pgrp='c2h'
                     call orcn(maxap3,a,b,d,d,atmchg,npop,nset,natoms,
     $                    dump,scr2,
     $                    scr1)
                     goto 2000
c
 210                 call reflct(maxap3,a,b,natoms,t,2)
                     call equiv(maxap3,a,b,atmchg,natoms,itst)
                     if (itst .eq. 0) goto 220
                     pgrp='c2v'
                     call orc2v(maxap3,a,b,natoms,atmchg)
                     goto 2000
c
 220                 pgrp='c2'
                     call orcn(maxap3,a,b,d,d,atmchg,npop,nset,natoms,
     $                    dump,scr2,
     $                    scr1)
                     goto 2000
c
c     *-----------------------*
c     symmetric top molecules
c     *-----------------------*
c
c     these molecules can belong to any axial point group, thus
c     only the cubic point groups (t, td, th, o, oh, i, ih) are
c     impossible.  however, execpt in rare cases the unique axis
c     is a rotation axis of order 3 or greater.  this routine is not
c     coded to identify the point group of these rare species.
c
c     align the unique axis with the z axis.
c
 300                 continue
                     if(dump) write(iout,1050)
                     if (abs(tst1) .lt. tol2) ixyz = 1
                     if (abs(tst2) .lt. tol2) ixyz = 2
                     if (abs(tst3) .lt. tol2) ixyz = 3
                     call putsym(maxap3,a,b,t,praxes(1,ixyz),numatm,3)
                     call oraxis(maxap3,a,b,natoms,atmchg,3)
c
c     test z for cn.
c     if n>1 then goto 330.
c     else quit.
c
                     call findcn(maxap3,natoms,a,b,d,atmchg,npop,
     $                    nset,3,norder)
                     if (norder .gt. 1) goto 330
                     symflg='axsymtop'
                     pgrp='c1'
                     goto 2000
c
c     the unique axis in a symmetric top molecule is a proper rotation
c     axis or order norder and is aligned with the cartesian z axis.
c
c     test z for s2n.
c     if no then continue.
c     else test for n dihedral planes.
c     if yes  then dnd.
c     if no   then s2n.
c     test for n c2 axes in the xy plane.
c     if no then continue at 400.
c     else test for a horizontal plane.
c     if yes  then dnh.
c     if no   then dn.
c
 330                 continue
                     theta = pi / norder
                     call rotate(maxap3,a,b,natoms,t,3,theta)
                     call reflct(maxap3,b,d,natoms,t,3)
                     call equiv(maxap3,a,d,atmchg,natoms,itst)
                     if (itst .eq. 0) goto 350
                     call findv(maxap3,a,b,d,natoms,npop,nset,
     $                    atmchg,itst)
                     if (itst .eq. 0) goto 340
                     pgrp(1:1)='d'
                     pgrp(2:)=itoc(norder)
                     pos=index(pgrp,' ')
                     pgrp(pos:pos)='d'
                     call ordn(maxap3,a,b,d,atmchg,npop,nset,natoms,
     $                    norder,dump)
                     goto 2000
c
 340                 pgrp(1:1)='s'
                     pgrp(2:)=itoc(2*norder)
                     call orcn(maxap3,a,b,d,d,atmchg,npop,nset,natoms,
     $                    dump,scr2,scr1)
                     goto 2000
c
 350                 call findc2(maxap3,a,b,d,npop,nset,atmchg,
     $                    natoms,itst)
                     if (itst .eq. 0) goto 400
                     call reflct(maxap3,a,b,natoms,t,3)
                     call equiv(maxap3,a,b,atmchg,natoms,itst)
                     if (itst .eq. 0) goto 360
                     pgrp(1:1)='d'
                     pgrp(2:)=itoc(norder)
                     pos=index(pgrp,' ')
                     pgrp(pos:pos)='h'
                     call ordn(maxap3,a,b,d,atmchg,npop,nset,
     $                    natoms,norder,dump)
                     if (norder .eq. 2) call ord2h(maxap3,a,b,
     $                    natoms,atmchg,ian)
                     goto 2000
c
 360                 pgrp(1:1)='d'
                     pgrp(2:)=itoc(norder)
                     call ordn(maxap3,a,b,d,atmchg,npop,nset,
     $                    natoms,norder,dump)
                     goto 2000
c
c     for a symmetric top molecule the possible point groups
c     have been limited to cnv, cnh, and cn.
c
 400                 continue
                     call findv(maxap3,a,b,d,natoms,npop,nset,
     $                    atmchg,itst)
                     if (itst .eq. 0) goto 410
                     pgrp(1:1) = 'c'
                     pgrp(2:)=itoc(norder)
                     pos=index(pgrp,' ')
                     pgrp(pos:pos) = 'v'
                     call ordn(maxap3,a,b,d,atmchg,npop,nset,
     $                    natoms,norder,dump)
                     if (norder .eq. 2) call orc2v(maxap3,a,b,natoms,
     $                    atmchg)
                     goto 2000
c
 410                 call reflct(maxap3,a,b,natoms,t,3)
                     call equiv(maxap3,a,b,atmchg,natoms,itst)
                     pgrp(1:1) = 'c'
                     pgrp(2:)=itoc(norder)
                     if(itst.ne.0) then
                        pos=index(pgrp,' ')
                        pgrp(pos:pos)='h'
                     endif
                     call orcn(maxap3,a,b,d,d,atmchg,npop,nset,natoms,
     $                    dump,scr2,scr1)
                     goto 2000
c
c     *-----------------------*
c     spherical top molecules
c     *-----------------------*
c
c     only the cubic point groups: t, td, th, o, oh, i, and, ih are
c     possible.  no provision is made in the subsequent code for the
c     the possibility that a spherical top molecule may belong to any
c     other point group.
c
c     find the highest order proper rotation axis and align it with
c     the z-axis.
c
 500                 continue
                     if (dump) write(iout,1060)
                     call sphere(maxap3,natoms,a,b,d,atmchg,nset,npop,
     $                    norder,dump,
     $                    scr1,scr2,scr3,scr4)
                     if (norder .ne. 0) goto 510
                     symflg='axsphtop'
                     pgrp='c1'
                     goto 2000
 510                 ihop = norder - 2
                     goto (520,560,580), ihop
c
c     a spherical top molecule has a two-fold axis aligned with
c     z and is t, td, or th.
c
c     test for a center of inversion.
c     if yes then th.
c     else test for a vertical plane.
c     if yes then td.
c     else t.
c
 520                 continue
                     call invert(maxap3,a,b,natoms,t)
                     call equiv(maxap3,a,b,atmchg,natoms,itst)
                     if (itst .eq. 0) goto 540
                     pgrp='th'
                     goto 550
c
 540                 call findv(maxap3,a,b,d,natoms,npop,nset,atmchg,
     $                    itst)
                     pgrp(1:1)='t'
                     if (itst .eq. 0) goto 550
                     pgrp(2:2)='d'
                     call rotate(maxap3,a,b,numatm,t,3,piovr4)
                     call move(maxap3,b,a,numatm)
 550                 call ordn(maxap3,a,b,d,atmchg,npop,nset,natoms,
     $                    2,dump)
                     goto 2000
c
c     a spherical top molecule has a four-fold axis aligned with z and
c     is either o or oh.
c
c     test for a center of inversion.
c     if yes then oh.
c     else o.
c
 560                 continue
                     call invert(maxap3,a,b,natoms,t)
                     call equiv(maxap3,a,b,atmchg,natoms,itst)
                     pgrp(1:1) = 'o'
                     if (itst .ne. 0) pgrp(2:2) = 'h'
                     call ordn(maxap3,a,b,d,atmchg,npop,nset,natoms,2,
     $                    dump)
                     goto 2000
c
c     a spherical top molecule has a five-fold axis aligned with z and
c     is either i or ih.
c
c     test for a center of inversion.
c     if yes then ih.
c     else i.
c
 580                 continue
                     call invert(maxap3,a,b,natoms,t)
                     call equiv(maxap3,a,b,atmchg,natoms,itst)
                     pgrp(1:1)='i'
                     if (itst .ne. 0) goto 600
                     call orcn(maxap3,a,b,d,d,atmchg,npop,nset,natoms,
     $                    dump,scr2,scr1)
                     goto 2000
 600                 pgrp(2:2) = 'h'
                     call ordn(maxap3,a,b,d,atmchg,npop,nset,natoms,2,
     $                    dump)
c
c     exit.
c     if requested, calculate ad print the moments of charge for the
c     reoriented molecule.
c
 2000                continue
                     call secmom(maxap3,natoms,a,atmchg,prmom,praxes)
                     if (dump) write(iout,1010)
     $                    (prmom(i),i=1,3),
     $                    ((praxes(j,i),i=1,3),j=1,3)
                     return
                     end

*deck @(#)ztoc.f	5.1  11/6/94
      subroutine ztoc(nz,ianz,iz,bl,alph,bet,ttest,numtet,
     $                natoms,ian,atmchg,c,cz,a,b,d,alpha,beta)
c***begin prologue     ztoc
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           z-matrix, coordinates
c***author             binkley, et. al., gaussian82
c***source             @(#)ztoc.f	5.1   11/6/94
c***purpose            computes cartesian coordinates give the z-matrix.
c***description
c     call ztoc(nz,ianz,iz,bl,alph,bet,ttest,numtet,
c               natoms,ian,atmchg,c,cz,a,b,d,alpha,beta)
c
c     routine to compute the cartesian coordinates, given the
c     z-matrix.  this routine returns coordinates both with,
c     and without the dummy atoms.
c
c     arguments:
c
c     nz     ... number of lines in the z-matrix.
c     ianz   ... the atomic numbers of the z-matrix centers.
c     iz     ... the integer components of the z-matrix.
c     bl     ... the bond-lengths from the z-matrix.
c     alph   ... the bond-angles from the z-matrix.  these angles
c                must be in radians.
c     bet    ... the dihedral angles from the z-matrix.  like
c                alph, these angles must also be in radians.
c     ttest  ... logical flag to enable testing for tetrahedral angles.
c                this feature is useful in obtaining exact tetrahedral
c                angles.  if any are found and this flag is set "true",
c                then exact values are used and a message is printed
c                indicating how many angles were changed.  the values
c                in alph and/or bet are updated.
c     natoms ... number of atoms (dummies removed), computed by this
c                routine.
c     ian    ... vector of length natoms, will receive atomic
c                numbers with dummies compressed out.
c     atmchg ... vector of length natoms, will receive atomic
c                charges with dummies compressed out.
c     c      ... coordinates, dummies compressed out.  stored
c                (x,y,z) for each atom.
c     cz     ... same as c, but with dummies still in.
c                the atomic number list, with dummies intact,
c                can be obtained from ianz.
c     a      ... scratch vector of length nz.
c     b      ... scratch vector of length nz.
c     d      ... scratch vector of length nz.
c     alpha  ... scratch vector of length nz.
c     beta   ... scratch vector of length nz.
c
c     this routine is dimension free, in the sense that any
c     restrictions are imposed by the calling routine.
c
c***references
c***routines called    vec(m202), vprod(m202), rzero(util),
c                      lnkerr(mdutil)
c***end prologue       ztoc
      implicit real*8(a-h,o-z)
      logical test,ttest,error,vecerr
      dimension ianz(nz),iz(4,nz),bl(nz),alpha(nz),beta(nz),ian(nz),
     $          atmchg(nz),c(nz),cz(nz),a(nz),b(nz),d(nz),alph(nz),
     $          bet(nz)
      dimension u1(3),u2(3),u3(3),u4(3),vj(3),vp(3),v3(3)
      common/io/inp,iout
      data zero,one,two,three/0.0d+00,1.0d+00,2.0d+00,3.0d+00/
      data four,f180/4.0d+00,1.80d+02/
      data tenm5,tenm6,tenm12/1.0d-05,1.0d-06,1.0d-12/
      data tetdat,toldat/109.471d+00,1.0d-03/
      save zero,one,two,three,four,f180
      save tenm5,tenm6,tenm12,tetdat,toldat
c
 1010 format(' invalid beta angle,type(z4),on z-matrix card:',i4)
 1020 format(' reference made to an undefined center on z-matrix card:',
     $        i4)
 1030 format(' multiple references to a center on the same card.'
     $        //' z-matrix card number:',i4)
 1040 format(' incipient floating point error detected on z-matrix card'
     $       //'number:',i4)
 1050 format (' angle alpha is outside valid range (0-180) on z-matrix'
     $        //'card number:',i4)
 1060 format(' bond length on z-matrix card number',i4,
     $       ' is not positive.')
 1070 format (' angle beta is outside valid range (0-180) on z-matrix'
     $        //'card number:',i4)
 2001 format(1x,i4,' tetrahedral angles replaced.')
c
c
c
      test(x) = abs(x-tettst).lt.tettol
c
c
c     check for nonsense in the connectivity.
c
      error=.false.
      if(nz.gt.1) then
         do 30 i=2,nz
            if(iz(1,i).ge.i.or.iz(1,i).le.0) then
               error=.true.
               write(iout,1020) i
            endif
c
            if(i.gt.2) then
               if(iz(2,i).ge.i.or.iz(2,i).le.0) then
                  error=.true.
                  write(iout,1020) i
               endif
               if(iz(2,i).eq.iz(1,i)) then
                  error=.true.
                  write(iout,1030) i
               endif
c
               if(i.gt.3) then
                  if(iz(3,i).ge.i.or.iz(3,i).le.0) then
                     error=.true.
                     write(iout,1020) i
                  endif
                  if(iz(3,i).eq.iz(1,i).or.iz(3,i).eq.iz(2,i)) then
                     error=.true.
                     write(iout,1030) i
                  endif
                  if(abs(iz(4,i)).gt.1) then
                     error=.true.
                     write(iout,1010) i
                  endif
               endif
            endif
   30    continue
      endif
      if(error) call lnkerr(' ')
c
c     compute local values of pi-related constants.
      call iosys('read real pi from rwf',1,pi,0,' ')
      torad=pi/f180
c
c     set up for laundering tetrahedral angles.
      if(ttest) then
         tetang=acos(-one/three)
         tettst=tetdat*torad
         tettol=toldat*torad
      endif
c
c     zero temporary coordinate array cz, move angles to local arrays,
c     test for alpha or beta out of range(0-180 degrees),  negative
c     bond lengths, and (optionally) test for tetrahedral angles.
      nz3=3*nz
      call rzero(cz,nz3)
c
      numtet=0
      if(nz.gt.1) then
         do 40 i=2,nz
            if(bl(i).le.zero) then
               error=.true.
               write(iout,1060) i
            endif
            if(i.gt.2) then
               alpha(i)=alph(i)
               if(alpha(i).lt.zero.or.alpha(i).gt.pi) then
                  error=.true.
                  write(iout,1050) i
               endif
               if(ttest.and.test(alpha(i))) then
                  alpha(i)=tetang
                  alph(i)=tetang
                  numtet=numtet+1
               endif
            endif
            if(i.gt.3) then
               beta(i)=bet(i)
               if(iz(4,i).ne.0
     $                    .and.(beta(i).lt.zero.or.beta(i).gt.pi)) then
                  error=.true.
                  write(iout,1070) i
               endif
               if(ttest.and.test(beta(i))) then
                  beta(i)=tetang
                  bet(i)=tetang
                  numtet=numtet+1
               endif
            endif
   40    continue
      endif
      if(numtet.ne.0) write (iout,2001) numtet
      if(error) call lnkerr(' ')
c
c     convert to cartesian coordinates
c     center 1 goes at the origin, center 2 along z, center 3 in xz.
c
c     z-coordinate, center 2.
      cz(6)=bl(2)
c
      if(nz.ge.3) then
c        x-coordinate, center 3.
         cz(7)=bl(3)*sin(alpha(3))
         if(iz(1,3).eq.1) then
c           z-coordinate, center3.
            cz(9)=bl(3)*cos(alpha(3))
         else
c           z-coordinate, center3, as a function of z-coordinate, center 2.
            cz(9)=cz(6)-bl(3)*cos(alpha(3))
         endif
      endif
c
      if(nz.ge.4) then
c        beware of linear molecules.
         do 50 i=4,nz
            ind3=(i-1)*3
            if(abs(cz(1+ind3-3)).ge.tenm5) goto 60
            cz(1+ind3)=bl(i)*sin(alpha(i))
            itemp=(iz(1,i)-1)*3
            jtemp=(iz(2,i)-1)*3
            cz(3+ind3)=cz(3+itemp)-bl(i)*cos(alpha(i))
     $                 *sign(one,cz(3+itemp)-cz(3+jtemp))
   50    continue
   60    k=i
         if(k.le.nz) then
c
c           we have handled the linear portions. now get to work.
            do 95 j=k,nz
               jnd3=(j-1)*3
               dcaj=cos(alpha(j))
               dsaj=sin(alpha(j))
               dcbj=cos(beta(j))
               dsbj=sin(beta(j))
               if(iz(4,j).eq.0) then
                  call vec(tenm6,vecerr,u1,cz,iz(2,j),iz(3,j))
                     if(vecerr) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  call vec(tenm6,vecerr,u2,cz,iz(1,j),iz(2,j))
                     if(vecerr) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  call vprod(vp,u1,u2)
                  arg=one-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))**2
                     if(arg.lt.zero) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  r=sqrt(arg)
                     if(r.lt.tenm6) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  do 70 i=1,3
   70                u3(i)=vp(i)/r
                  call vprod(u4,u3,u2)
                  do 80 i=1,3
                     vj(i)=bl(j)*
     $                     (-u2(i)*dcaj+u4(i)*dsaj*dcbj+u3(i)*dsaj*dsbj)
                     itemp=(iz(1,j)-1)*3
   80                cz(i+jnd3)=vj(i)+cz(i+itemp)
c
c
               else if((abs(iz(4,j))-1).eq.0) then
                  call vec(tenm6,vecerr,u1,cz,iz(1,j),iz(3,j))
                     if(vecerr) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  call vec(tenm6,vecerr,u2,cz,iz(2,j),iz(1,j))
                     if(vecerr) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  zeta=-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))
                  denom=one-zeta**2
                     if(abs(denom).lt.tenm6) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  a(j)=(-dcbj+zeta*dcaj)/denom
                  b(j)=(dcaj-zeta*dcbj)/denom
                  r=zero
                  gamma=pi/two
                  if (abs(zeta).ge.tenm6) then
                     if(zeta.lt.zero) r=pi
                        if(denom.lt.zero) then
                           error=.true.
                           write(iout,1040) j
                           call lnkerr(' ')
                        endif
                     gamma=atan(sqrt(denom)/zeta)+r
                  endif
                  d(j)=zero
                  if(abs(gamma+alpha(j)+beta(j)-two*pi)-tenm6.ge.zero)
     $               then
                     arg=(one+ a(j)*dcbj - b(j)*dcaj) / denom
                        if(arg.lt.zero) then
                           error=.true.
                           write(iout,1040) j
                           call lnkerr(' ')
                        endif
                     d(j)=float(iz(4,j))*sqrt(arg)
                  endif
                  call vprod(v3,u1,u2)
                  do 85 i=1,3
                     u3(i)=a(j)*u1(i)+b(j)*u2(i)+d(j)*v3(i)
                     vj(i)=bl(j)*u3(i)
                     itemp=(iz(1,j)-1)*3
   85                cz(i+jnd3)=vj(i)+cz(i+itemp)
c
c
               else
                  call vec(tenm6,vecerr,u1,cz,iz(1,j),iz(3,j))
                     if(vecerr) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  call vec(tenm6,vecerr,u2,cz,iz(2,j),iz(1,j))
                     if(vecerr) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  zeta=-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))
                  call vprod(v3,u1,u2)
                  v3mag=sqrt(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
                  denom=one - zeta**2
                     if(abs(denom).lt.tenm6) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  a(j)=v3mag*dcbj / denom
                  arg=(one-dcaj*dcaj-a(j)*dcbj*v3mag) / denom
                     if(arg.lt.zero) then
                        error=.true.
                        write(iout,1040) j
                        call lnkerr(' ')
                     endif
                  b(j)=sqrt(arg)
                  if((iz(4,j)-2).eq.0) then
                     d(j)=b(j)*zeta+dcaj
                  else
                     b(j)=-b(j)
                  endif
                  do 90 i=1,3
                     u3(i)=b(j)*u1(i)+d(j)*u2(i)+a(j)*v3(i)
                     vj(i)=bl(j)*u3(i)
                     itemp=(iz(j,1)-1)*3
   90                cz(i+jnd3)=vj(i)+cz(i+itemp)
               endif
   95       continue
         endif
      endif
c
c
c     eliminate dummy atoms. dummy atoms are characterized by
c     negative atomic numbers.  ghost atoms have zero atomic
c     numbers. ghost atoms are not eliminated.
      natoms=0
      iaind=0
      naind=0
      do 110 i=1,nz
         if(ianz(i).ge.0) then
            natoms=natoms+1
            ian(natoms)=ianz(i)
            atmchg(natoms)=float(ianz(i))
            c(1+naind)=cz(1+iaind)
            c(2+naind)=cz(2+iaind)
            c(3+naind)=cz(3+iaind)
            naind=naind+3
         endif
         iaind=iaind+3
  110 continue
c
c     'tidy' up coordinates.
      do 120 i=1,3*natoms
         if(abs(c(i)).le.tenm12) c(i)=zero
  120 continue
c
c
c
      return
      end

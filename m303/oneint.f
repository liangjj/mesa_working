*deck @(#)oneint.f	5.2  2/5/95
      subroutine oneint(c,ex,z,cont,s,ptprim,noprim,nocont,ptcont,
     #                  nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     #                  start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     #                  nx,ny,nz,minmom,ops,ds,nderiv,dhs)
c
c***module to form the overlap, kinetic-energy and potential-energy
c   one-electron integrals over generally contracted gaussian basis
c   sets.
c
c paul saxe               23 july 1984                  lanl
c
c byron lengsfield         5 august 1991                llnl
c                          modifications of nacme
c
      implicit integer (a-z)
c
      character*(*) ops
      real*8 c(3,nat),ex(nprim),z(maxcor),s(nnp),cont(ncont)
      real*8 zan(nat),ds(nnp,3,nat)
c..bhl
      real*8 dhs(*)
c..bhl
      real*8 h(21),wt(21),sdot,eonel,rjunk
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      logical dolp,logkey
c
      common/io/inp,iout
c
      data h   /0.0d+00
     2,         -.707106781186548d+00,  0.707106781186548d+00
     3,         -1.22474487139159d+00,0.0d+00,1.22474487139159d+00
     4,         -1.65068012388578d+00, -0.524647623275290d+00
     4,          0.524647623275290d+00, 1.65068012388578d+00
     5,   -2.02018287045609d+00,-0.958572464613819d+00,0.0d+00
     5,          0.958572464613819d+00, 2.02018287045609d+00
     6,         -2.350604973674d+00  , -1.335849074014d+00
     6,         -0.436077411928d+00  ,  0.436077411928d+00
     6,          1.335849074014d+00  ,  2.350604973674d+00/
      data wt  /1.77245385090552d+00
     2,         0.8862269254528d+00  ,  0.8862269254528d+00
     3,         0.2954089751509d+00  ,  1.181635900604d+00
     3,         0.2954089751509d+00
     4,         8.131283544725d-02   ,  8.049140900055d-01
     4,         8.049140900055d-01   ,  8.131283544725d-02
     5,         1.995324205905d-02   ,  3.936193231522d-01
     5,         9.453087204829d-01   ,  3.936193231522d-01
     5,         1.995324205905d-02
     6,         4.530009905509d-03   ,  1.570673203229d-01
     6,         7.246295952244d-01   ,  7.246295952244d-01
     6,         1.570673203229d-01   ,  4.530009905509d-03/
      data mxpts /21/
      save h,wt,mxpts
c
c
c     ----- zero the space to accumulate derivative overlap integrals
c
      if (nderiv.eq.1) then
        nbf2=nbasis*nbasis
        call rzero(ds,nnp*3*nat)
        call rzero(dhs,nbf2*3*nat)
      end if
c
c     ----- form the overlap integrals -----
c
      do 9 iatom=1,nat
         do 8 jatom=1,iatom
            do 7 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 7
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 6 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 6
c
                  imax=maxmom(itype)
                  jmax=maxmom(jtype)
                  nprimi=noprim(iatom,itype)
                  nprimj=noprim(jatom,jtype)
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c     ----- allocate core for temporary vectors, etc. -----
c
                  a=1
                  b=a+npint*(nderiv*2+1)
                  alpha=b+npint
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=1
                  xyz=max(expon+npint,prmint+lenblk*npint*(nderiv*2+2))
                  top1=xyz+3*npint*(imax+1)*(jmax+1)*(nderiv+1)
c
                  conint=prmint+npint*lenblk*(nderiv*2+2)
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  top=max(top1,top2)
                  call getscm(top,z,junk,'overlap integrals',0)
c
c     ----- form the two-dimensional integrals -----
c
                  call sints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       ex,z(alpha),c,z(ainv),z(xyza),
     #                       z(expon),z(xyz),npint,nprim,nat,
     #                       nbtype,h,wt,mxpts,z(a),z(b),nderiv)
c
c     ----- form the primitive integrals -----
c
                  call fmoned(z(prmint),z(xyz),npint,lenblk,imax,jmax,
     #                   mintyp(itype),mintyp(itype)+nocart(itype)-1,
     #                   mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     #                   nx,ny,nz,nderiv)
c
c     ----- transform to contracted functions -----
c
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c     ----- transfer integrals to total array -----
c
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     #                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
c     ----- transform derivative integrals, and sum into array -----
c
c..bhl                   if (nderiv.eq.1.and.iatom.ne.jatom) then
c..bhl                   call transd
c
                  if (nderiv.eq.1) then
                     call tranhs(z(prmint),z(conint),nprimi,nprimj,
     #                      nconti,ncontj,itype,jtype,
     #                      cont(ptcont(iatom,itype)),
     #                      cont(ptcont(jatom,jtype)),
     #                      z(tmp1),len1,lenblk,minmom(itype),
     #                      maxmom(itype),minmom(jtype),
     #                      maxmom(jtype),nocart,iatom,jatom,
     #                      ds,nnp,nat,start,nbtype,nobf,dhs,nbasis)
                  end if
c
    6          continue
    7       continue
    8    continue
    9 continue
c
c     ----- print the integrals -----
c
      if (logkey(ops,'print=gradient=s',.false.,' ')) then
         write (iout,10)
   10    format ('1',//,' the overlap integrals')
         call print(s,nnp,nbasis,iout)
         do 13 iatom=1,nat
            do 12 coord=1,3
               write (iout,11) iatom,coord
   11          format ('1',//,'   derivative overlap, atom',i3,
     #                 ';  coordinate ',i1)
               call print(ds(1,coord,iatom),nnp,nbasis,iout)
   12       continue
   13    continue
      end if
c
c     ----- put overlap integrals on the read-write file -----
c
      call iosys('write real "ao derivative overlap integrals" '//
     $     'on rdints',nnp*3*nat,ds,0,' ')
c..bhl.unicos
      call iosys('write real "ao derivative overlap integrals" '//
     $     'on rwf',nnp*3*nat,ds,0,' ')
c..bhl aug 5, 1991
c     write(iout,*)' storing half derivative overlap integrals '
      call iosys('write real "ao half_deriv overlap integrals" '//
     $     'on rwf',nbf2*3*nat,dhs,0,' ')
c..bhl aug 5, 1991
c..bhl.unicos
c
c     ----- form the kinetic-energy integrals -----
c
      call rzero(ds,nnp*3*nat)
c
      do 19 iatom=1,nat
         do 18 jatom=1,iatom
            do 17 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 17
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 16 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 16
c
                  imax=maxmom(itype)
                  jmax=maxmom(jtype)
                  nprimi=noprim(iatom,itype)
                  nprimj=noprim(jatom,jtype)
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c     ----- allocate core for temporary vectors, etc. -----
c
                  a=1
                  b=a+npint*(nderiv*2+1)
                  bp2=b+npint
                  bm2=bp2+npint
                  alpha=bm2+npint
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=1
                  xyz=max(expon+npint,prmint+lenblk*npint*(2*nderiv+2))
                  t1=xyz+3*npint*(imax+1)*(jmax+1)*2*(nderiv+1)
                  top1=t1+npint
c
                  conint=prmint+npint*lenblk*(2*nderiv+2)
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  top=max(top1,top2)
                  call getscm(top,z,junk,'kinetic integrals',0)
c
c     ----- form the two-dimensional integrals -----
c
                  call tints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       ex,z(alpha),c,z(ainv),z(xyza),
     #                       z(expon),z(xyz),npint,nprim,nat,
     #                       nbtype,h,wt,mxpts,z(a),z(b),z(bp2),z(bm2),
     #                       nderiv)
c
c     ----- form the primitive integrals -----
c
                  call fmt(z(prmint),z(xyz),npint,lenblk,
     #                     imax,jmax,z(t1),
     #                   mintyp(itype),mintyp(itype)+nocart(itype)-1,
     #                   mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     #                   nx,ny,nz,nderiv)
c
c     ----- transform to contracted functions -----
c
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c     ----- transfer integrals to total array -----
c
                 call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     #                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
c     ----- transform derivative integrals, and sum into array -----
c
                  if (nderiv.eq.1.and.iatom.ne.jatom) then
                     call transd(z(prmint),z(conint),nprimi,nprimj,
     #                           nconti,ncontj,itype,jtype,
     #                           cont(ptcont(iatom,itype)),
     #                           cont(ptcont(jatom,jtype)),
     #                           z(tmp1),len1,lenblk,minmom(itype),
     #                           maxmom(itype),minmom(jtype),
     #                           maxmom(jtype),nocart,iatom,jatom,
     #                           ds,nnp,nat,start,nbtype,nobf)
                  end if
c
   16          continue
   17       continue
   18    continue
   19 continue
c
c     ----- print the integrals -----
c
      if (logkey(ops,'print=gradient=t',.false.,' ')) then
         write (iout,20)
   20    format (/,' the kinetic-energy integrals')
         call print(s,nnp,nbasis,iout)
         do 23 iatom=1,nat
            do 22 coord=1,3
               write (iout,21) iatom,coord
   21          format ('1',//,'   derivative kinetic, atom',i3,
     #                 ';  coordinate ',i1)
               call print(ds(1,coord,iatom),nnp,nbasis,iout)
   22       continue
   23    continue
      end if
c
c     ----- form the potential-energy integrals -----
c
      do 29 iatom=1,nat
         do 28 jatom=1,iatom
            do 27 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 27
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 26 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 26
c
                  imax=maxmom(itype)
                  jmax=maxmom(jtype)
                  nroots=(imax+jmax+nderiv)/2+1
                  nprimi=noprim(iatom,itype)
                  nprimj=noprim(jatom,jtype)
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c     ----- allocate core for temporary vectors, etc. -----
c
                  a=1
                  b=a+npint*(nderiv*2+1)
                  aiaj=b+npint*(nderiv*2+1)
                  xyz0=aiaj+npint
                  xyz1=xyz0+npint*3
                  rysrt=xyz1+npint*3
                  ryswt=rysrt+npint*nroots
                  alpha=ryswt+npint*nroots
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=expon+npint
                  xyz=max(expon+npint,prmint+lenblk*npint*7)
                  t1=xyz+3*npint*(imax+1)*(jmax+1)*3
                  top1=t1+npint
c
                  conint=prmint+npint*lenblk*7
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  top=max(top1,top2)
                  call getscm(top,z,junk,'potential integrals',0)
c
c     ----- form the primitive integrals -----
c
                  call vints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       ex,z(alpha),c,z(ainv),z(xyza),
     #                       z(expon),z(xyz),npint,nprim,nat,
     #                       nbtype,h,wt,mxpts,z(a),z(b),z(aiaj),
     #                       z(xyz0),z(rysrt),z(ryswt),nroots,
     #                       z(prmint),lenblk,zan,z(xyz1),z(t1),
     #                 mintyp(itype),mintyp(itype)+nocart(itype)-1,
     #                 mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     #                 nx,ny,nz,nderiv,z(conint),
     #                           nconti,ncontj,itype,jtype,
     #                           cont(ptcont(iatom,itype)),
     #                           cont(ptcont(jatom,jtype)),
     #                           z(tmp1),len1,minmom(itype),
     #                           minmom(jtype),nocart,
     #                           ds,nnp,start,nobf)
c
c     ----- transform to contracted functions -----
c
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c     ----- transfer integrals to total array -----
c
                  call put1el(s,z(conint),start,iatom,jatom,itype,
     #                        jtype,nconti,ncontj,nnp,lenblk,nat,
     #                        nbtype,nobf)
c
   26          continue
   27       continue
   28    continue
   29 continue
c
c     --- print the potential integrals ---
      if (logkey(ops,'print=gradient=v',.false.,' ')) then
         write (iout,30)
   30    format (/,' the potential energy integrals')
         call print(s,nnp,nbasis,iout)
         do 33 iatom=1,nat
            do 32 coord=1,3
               write (iout,31) iatom,coord
   31          format ('1',//,'   derivative t+v, atom',
     #                 i3,';  coordinate ',i1)
               call print(ds(1,coord,iatom),nnp,nbasis,iout)
   32       continue
   33    continue
      end if
c
c     are there any ecp centers?
c
      dolp=.false.
      do 35 iatom=1,nat
         do 34 lptype=nbtype+1,ntypes
            if(noprim(iatom,lptype).gt.0) dolp=.true.
   34    continue
   35 continue
      if(dolp) then
c
c        ----- form the effective core potential integrals -----
c
         call ldata
         call ztab
         do 39 iatom=1,nat
            do 38 jatom=1,iatom
               do 37 itype=1,nbtype
                  if (noprim(iatom,itype).le.0) go to 37
                  if (iatom.ne.jatom) then
                     jtypmx=nbtype
                  else
                     jtypmx=itype
                  end if
                  do 36 jtype=1,jtypmx
                     if (noprim(jatom,jtype).le.0) go to 36
c
c                    find the basis function type that corresponds
c                    to the derivative shell block.
                     dimax=maxmom(itype)+nderiv
                     djmax=maxmom(jtype)+nderiv
                     dimin=minmom(itype)-nderiv
                     djmin=minmom(jtype)-nderiv
                     dimin=max(dimin,0)
                     djmin=max(djmin,0)
                     ditype=0
                     djtype=0
                     do 51 ii=1,nbtype
                        if(dimax.eq.maxmom(ii).and.dimin.eq.minmom(ii))
     $                       ditype=ii
                        if(djmax.eq.maxmom(ii).and.djmin.eq.minmom(ii))
     $                      djtype=ii
   51                continue
                     if(nderiv.ne.0.and.ditype.eq.0.or.djtype.eq.0)
     $                  call lnkerr('derivatives for this basis set '
     $                              //'unavailable.')
                     imax=maxmom(itype)
                     jmax=maxmom(jtype)
                     nprimi=noprim(iatom,itype)
                     nprimj=noprim(jatom,jtype)
                     nconti=nocont(iatom,itype)
                     ncontj=nocont(jatom,jtype)
                     npint=nprimi*nprimj
                     lenblk=nocart(itype)*nocart(jtype)
                     if(nderiv.ne.0) then
                        dlen=nocart(ditype)*nocart(djtype)
                     else
                        dlen=lenblk
                     endif
c
c        ----- allocate core for temporary vectors, etc. -----
c
                     prmint=1
                     ntpse=prmint+(npint*7*lenblk)
                     nlp=ntpse+10
                     zlp=nlp+100
                     clp=zlp+100
                     qq=clp+100
                     scr=qq+(11*9*9*npint)
                     top1=scr+npint*dlen
c
                     conint=top1
                     tmp1=conint+nconti*ncontj*lenblk
                     len1=nconti*nprimj
                     top2=tmp1+len1
c
                     top=max(top1,top2)
                     call getscm(top,z,junk,'ecp integrals',0)
c
c        ----- form the primitive integrals -----
c
                     call lpints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       nprim,ncont,nat,nbtype,ntypes,npint,lenblk,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       c,ex,noprim,ptprim,ptcont,cont,maxmom,
     #                       mintyp(itype),
     #                       mintyp(itype)+nocart(itype)-1,
     #                       mintyp(jtype),
     #                       mintyp(jtype)+nocart(jtype)-1,
     #                       nx,ny,nz,z(prmint),
     #                       z(ntpse),z(nlp),z(zlp),z(clp),z(qq),
     #                       z(scr),dlen,nderiv,mintyp(ditype),
     #                       mintyp(ditype)+nocart(ditype)-1,
     #                       mintyp(djtype),
     #                       mintyp(djtype)+nocart(djtype)-1,
     #                           z(conint),nconti,ncontj,itype,jtype,
     #                           cont(ptcont(iatom,itype)),
     #                           cont(ptcont(jatom,jtype)),
     #                           z(tmp1),len1,minmom(itype),
     #                           minmom(jtype),nocart,
     #                           ds,nnp,start,nobf)
c
c     ----- transform to contracted functions -----
c
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c        ----- transfer integrals to total array -----
c
                     call put1el(s,z(conint),start,iatom,jatom,itype,
     #                 jtype,nconti,ncontj,nnp,lenblk,nat,nbtype,
     #                 nobf)
c
   36             continue
   37          continue
   38       continue
 39      continue
      endif
c
c     ----- print the integrals -----
c
c     these are actually the entire one-electron part, not just
c     the ecp-contribution
      if (logkey(ops,'print=gradient=ecp',.false.,' ')) then
         write (iout,40)
   40    format (/,' the ecp integrals')
         call print(s,nnp,nbasis,iout)
         do 43 iatom=1,nat
            do 42 coord=1,3
               write (iout,41) iatom,coord
   41          format ('1',//,'   derivative t+v+ecp, atom',i3,
     #                 ';  coordinate ',i1)
               call print(ds(1,coord,iatom),nnp,nbasis,iout)
   42       continue
   43    continue
      end if
c
c     ----- put one-electron integrals on the integral file ----
c
      call iosys('write real "ao derivative one-electron integrals" '//
     $     'on rdints',nnp*3*nat,ds,0,' ')
c
c        ----- end timing of this section -----
c
  100 continue
c
      return
      end

*deck @(#)oneint.f	5.2  4/17/95
      subroutine oneint(c,ex,z,cont,s,ptprim,noprim,nocont,ptcont,
     #                  nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     #                  start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     #                  nx,ny,nz,minmom,ops,nderiv,grad,tempg,d,
     #                  f,ndmat,d2e,ld2e,nd2e,nd1e,hf,ci,mcscf,dft)
c
c***begin prologue     oneint
c***date written       840723   (yymmdd)
c***revision date      871121   (yymmdd)
c
c 21 november 1987     pws at lanl
c    adding sections for ci and mcscf derivatives.
c
c 14 november 1987     pws at lanl
c    adding second derivatives.
c
c 21 august 1986       pws at lanl
c    changing to handle open-shell cases by using generalized fock
c    f coefficients, which are half the occupancy of orbitals.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)oneint.f	5.2   4/17/95
c
c***purpose
c
c***description
c
c   module to form the overlap, kinetic-energy and potential-energy
c   one-electron integrals over generally contracted gaussian basis
c   sets.
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit integer (a-z)
c
      character*(*) ops
      real*8 d2e(nd2e)
      real*8 ld2e(nd2e)
      real*8 c(3,nat),ex(nprim),z(maxcor),s(nnp),cont(ncont)
      real*8 zan(nat),grad(3,nat),tempg(3,nat),d(nnp)
      real*8 h(21),wt(21),sdot,eonel,rjunk,f(ndmat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      logical dolp,logkey
      logical hf
      logical ci
      logical mcscf
      logical dft
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
c     ----- start timing if enabled -----
c
c
c     ----- zero the space to accumulate derivative overlap integrals
c
      if (nderiv.eq.2) then
         call rzero(ld2e,nd2e)
         call rzero(tempg,3*nat)
         nint=10
      else if (nderiv.eq.1) then
         call rzero(tempg,3*nat)
         nint=4
      end if
c
c     ----- fetch the lagrangian matrix -----
c
      if (mcscf) then
         call iosys('read real "mcscf ao lagrangian" from rwf',
     $        nnp,d,0,' ')
      else if (ci) then
         go to 700
      else if (dft) then
         call iosys('read real "dft nrg weighted density matrix" '//
     $        'from rwf',nnp,d,0,' ')
      else 
         call iosys('read real "scf ao lagrangian" from rwf',
     $        nnp,d,0,' ')
         
      end if
c
c     ----- double off-diagonals to account for lij and lji -----
c
      ij=0
      do 201 i=1,nbasis
         do 200 j=1,i-1
            ij=ij+1
            d(ij)=d(ij)*2.0d+00
  200    continue
         ij=ij+1
  201 continue
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
                  xyz=max(expon+npint,prmint+lenblk*npint*nint)
                  top1=xyz+3*npint*(imax+1)*(jmax+1)*(nderiv+1)
c
                  conint=prmint+npint*lenblk*nint
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  if (top1.gt.maxcor.or.top2.gt.maxcor) then
                     call lnkerr('m702: oneint..not enough core for s')
                  end if
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
     #                   nx,ny,nz,nderiv,nint)
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
                  if (nderiv.ge.1.and.iatom.ne.jatom) then
                     call transd(z(prmint),z(conint),nprimi,nprimj,
     #                           nconti,ncontj,itype,jtype,
     #                           cont(ptcont(iatom,itype)),
     #                           cont(ptcont(jatom,jtype)),
     #                           z(tmp1),len1,lenblk,minmom(itype),
     #                           maxmom(itype),minmom(jtype),
     #                           maxmom(jtype),nocart,iatom,jatom,
     #                           d,nnp,nat,start,nbtype,nobf,tempg)
                  end if
c
c     ----- transform second derivative integrals, and sum into array -----
c
                  if (nderiv.ge.2.and.iatom.ne.jatom) then
                     call transf(z(prmint),z(conint),nprimi,nprimj,
     #                           nconti,ncontj,itype,jtype,
     #                           cont(ptcont(iatom,itype)),
     #                           cont(ptcont(jatom,jtype)),
     #                           z(tmp1),len1,lenblk,minmom(itype),
     #                           maxmom(itype),minmom(jtype),
     #                           maxmom(jtype),nocart,iatom,jatom,
     #                           d,nnp,nat,start,nbtype,nobf,nint,
     #                           ld2e,nd2e)
                  end if
c
c
    6          continue
    7       continue
    8    continue
    9 continue
c
c     ----- form the gradient contribution -----
c
      do 203 iatom=1,nat
         do 202 coord=1,3
            tempg(coord,iatom)=-2.0d+00*tempg(coord,iatom)
            grad(coord,iatom)=grad(coord,iatom)+tempg(coord,iatom)
  202    continue
  203 continue
c
c     ----- and the second-derivative contribution -----
c
      if (nderiv.ge.2) then
         do 303 i=1,nd1e
            do 302 j=1,i
               ij=i*(i-1)/2+j
               ld2e(ij)=-2.0d+00*ld2e(ij)
               d2e(ij)=d2e(ij)+ld2e(ij)
  302       continue
  303    continue
      end if
c
      if (logkey(ops,'print=gradient=overlap',.false.,' ')) then
         write (iout,204)
  204    format ('1     the overlap contribution to the scf ',
     #           'gradients:')
         call matout(tempg,3,nat,3,nat,iout)
c
         if (nderiv.eq.2) then
            write (iout,304)
  304       format (/,'     the overlap contribution to the scf ',
     #              'force constants')
            call print(ld2e,nd2e,nd1e,iout)
         end if
      end if
c
 700  continue
c
c     ----- form the kinetic-energy integrals -----
c
      eonel=0.0d+00
c
c     ----- add the various density matrices together, weighted by f
c
      if (mcscf) then
         call iosys('read real "mcscf ao 1pdm" from rwf',nnp,d,0,' ')
      else if (ci) then
         call iosys('read real "ci ao 1pdm" from rwf',nnp,d,0,' ')
      else
         call rzero(d,nnp)
         call iosys('rewind "hf density matrix" on rwf',0,0,0,' ')
         do 1000 dmat=1,ndmat
            call iosys('read real "hf density matrix" from rwf '//
     #           'without rewinding',nnp,s,0,' ')
            do 999 i=1,nnp
               d(i)=d(i)+2.0d+00*f(dmat)*s(i)
 999        continue
 1000    continue
      end if
c
c     ----- double off-diagonals to account for dij and dji -----
c
      ij=0
      do 206 i=1,nbasis
         do 205 j=1,i-1
            ij=ij+1
            d(ij)=d(ij)*2.0d+00
  205    continue
         ij=ij+1
  206 continue
c
c     ----- zero space for the gradients and second derivatives -----
c
      if (nderiv.eq.2) then
         call rzero(ld2e,nd2e)
         call rzero(tempg,3*nat)
         nint=10
      else if (nderiv.eq.1) then
         call rzero(tempg,3*nat)
         nint=4
      end if
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
                  xyz=max(expon+npint,prmint+lenblk*npint*nint)
                  t1=xyz+3*npint*(imax+1)*(jmax+1)*2*(nderiv+1)
                  top1=t1+npint
c
                  conint=prmint+npint*lenblk*nint
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  if (top1.gt.maxcor.or.top2.gt.maxcor) then
                     call lnkerr('m702: oneint..not enough core for t')
                  end if
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
     #                   nx,ny,nz,nderiv,nint)
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
                  if (nderiv.ge.1.and.iatom.ne.jatom) then
                     call transd(z(prmint),z(conint),nprimi,nprimj,
     #                    nconti,ncontj,itype,jtype,
     #                    cont(ptcont(iatom,itype)),
     #                    cont(ptcont(jatom,jtype)),
     #                    z(tmp1),len1,lenblk,minmom(itype),
     #                    maxmom(itype),minmom(jtype),
     #                    maxmom(jtype),nocart,iatom,jatom,
     #                    d,nnp,nat,start,nbtype,nobf,tempg)
                  end if
c
c     ----- transform second derivative integrals, and sum into array -----
c
                  if (nderiv.ge.2.and.iatom.ne.jatom) then
                     call transf(z(prmint),z(conint),nprimi,nprimj,
     #                           nconti,ncontj,itype,jtype,
     #                           cont(ptcont(iatom,itype)),
     #                           cont(ptcont(jatom,jtype)),
     #                           z(tmp1),len1,lenblk,minmom(itype),
     #                           maxmom(itype),minmom(jtype),
     #                           maxmom(jtype),nocart,iatom,jatom,
     #                           d,nnp,nat,start,nbtype,nobf,nint,
     #                           ld2e,nd2e)
                  end if
c
   16          continue
   17       continue
   18    continue
   19 continue
c
c     ----- form the gradient contribution -----
c
      eonel=eonel+sdot(nnp,d,1,s,1)
      do 208 iatom=1,nat
         do 207 coord=1,3
            grad(coord,iatom)=grad(coord,iatom)+tempg(coord,iatom)
  207    continue
  208 continue
c
c     ----- and the second-derivative contribution -----
c
      if (nderiv.ge.2) then
         do 308 i=1,nd1e
            do 307 j=1,i
               ij=i*(i-1)/2+j
               ld2e(ij)=ld2e(ij)
               d2e(ij)=d2e(ij)+ld2e(ij)
 307        continue
 308     continue
      end if
c
c
      if (logkey(ops,'print=gradient=kinetic',.false.,' ')) then
         write (iout,209)
  209    format (//,'     the kinetic-energy contribution to the scf ',
     #           'gradients:',//)
         call matout(tempg,3,nat,3,nat,iout)
c
         if (nderiv.eq.2) then
            write (iout,309)
 309        format (/,'     the kinetic-energy contribution to the ',
     $           'scf force constants')
            call print(ld2e,nd2e,nd1e,iout)
         end if
      end if
c
c     ----- form the potential-energy integrals -----
c
c
c     ----- zero space for the gradients and second derivatives -----
c
      if (nderiv.eq.2) then
         call rzero(ld2e,nd2e)
         call rzero(tempg,3*nat)
         nxyz=6
         nint=28
      else if (nderiv.eq.1) then
         call rzero(tempg,3*nat)
         nint=7
         nxyz=3
      end if
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
                  xyz=max(expon+npint,prmint+lenblk*npint*nint)
                  t1=xyz+3*npint*(imax+1)*(jmax+1)*nxyz
                  top1=t1+npint
c
                  conint=prmint+npint*lenblk*nint
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
             if (top1.gt.maxcor.or.top2.gt.maxcor) then
                write(iout,*) 'top1,top2,maxcor',top1,top2,maxcor
                call lnkerr('m702: oneint..not enough core for v')
             end if
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
     #                           d,nnp,start,nobf,tempg,nxyz,ld2e,nd2e)
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
c     ----- form the gradient contribution -----
c
      eonel=eonel+sdot(nnp,d,1,s,1)
      do 211 iatom=1,nat
         do 210 coord=1,3
            grad(coord,iatom)=grad(coord,iatom)+tempg(coord,iatom)
  210    continue
  211 continue
c
c     ----- and the second-derivative contribution -----
c
      if (nderiv.ge.2) then
         do 311 i=1,nd1e
            do 310 j=1,i
               ij=i*(i-1)/2+j
               d2e(ij)=d2e(ij)+ld2e(ij)
 310        continue
 311     continue
      end if
c
      if (logkey(ops,'print=gradient=potential',.false.,' ')) then
         write (iout,212)
  212    format (//,'     the potential-energy contribution to the ',
     #           'scf gradients:',//)
         call matout(tempg,3,nat,3,nat,iout)
c
         if (nderiv.eq.2) then
            write (iout,312)
 312        format (/,'     the potential-energy contribution to the ',
     $           'scf force constants')
            call print(ld2e,nd2e,nd1e,iout)
         end if
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
         call rzero(tempg,3*nat)
         if (nderiv .eq. 2) then
            call rzero(ld2e,nd2e)
         endif
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
     $                      ditype=ii
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
                     if (nderiv.eq.2) then
                        ntpse=prmint+(npint*28*lenblk)
                     else
                        ntpse=prmint+(npint*7*lenblk)
                     endif
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
                     if (top1.gt.maxcor.or.top2.gt.maxcor) then
                        write(iout,*) 'top1,top2,maxcor',
     $                                 top1,top2,maxcor
                        call lnkerr('m702: oneint..not enough core'
     #                              //' for lp')
                     end if
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
     #                           d,nnp,start,nobf,tempg,ld2e,nd2e)
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
     #                        jtype,nconti,ncontj,nnp,lenblk,nat,nbtype,
     #                        nobf)
c
   36             continue
   37          continue
   38       continue
   39    continue
c
c     ----- form the gradient contribution -----
c
         eonel=eonel+sdot(nnp,d,1,s,1)
         do 221 iatom=1,nat
            do 220 coord=1,3
                  grad(coord,iatom)=grad(coord,iatom)+tempg(coord,iatom)
  220       continue
  221    continue
c
c     ----- and the second-derivative contribution -----
c
      if (nderiv.ge.2) then
         do 321 i=1,nd1e
            do 320 j=1,i
               ij=i*(i-1)/2+j
               d2e(ij)=d2e(ij)+ld2e(ij)
 320        continue
 321     continue
      end if
c
c        print the derivative integrals.
c
         if (logkey(ops,'print=gradient=ecp',.false.,' ')) then
            write (iout,223)
 223        format (//,'     the effective core potential contribution',
     #              ' to the scf gradients:',//)
            call matout(tempg,3,nat,3,nat,iout)
            if (nderiv.eq.2) then
               write (iout,323)
 323           format (/,'     the ecp contribution to the ',
     $              'scf force constants')
               call print(ld2e,nd2e,nd1e,iout)
            end if
         end if
c
      endif
c
      if (logkey(ops,'print=gradient=one-electron',.false.,' ')) then
         write (iout,213)
  213    format ('1',//,'     the one-electron and nuclear repulsion ',
     #           'contribution to the scf gradients:',//)
         call matout(grad,3,nat,3,nat,iout)
      end if
c
c     ----- check up on the one-electron energies -----
c
      if (logkey(ops,'hf=quadratic',.false.,' ').or.
     $     logkey(ops,'mcscf',.false.,' ').or.
     $     logkey(ops,'ci',.false.,' ')) then
c
c        ----- mcscf does not cleanly separate one- and two-electron
c              energies (core fock business) so procrastinate
c
         call iosys('write real "mcscf m702 1e energy" to rwf',
     $              1,eonel,0,' ')
         if(logkey(ops,'m702=check-energy',.false.,' ')) then
            write(iout,9214) eonel
 9214       format(5x,'checking the 1e- energy:'
     $           /,10x,'one-electron energy:',4x,f14.8)
         end if
      else
c
c        --- these will be tested later down the line.
         call iosys('read real "hf 1e energy" from rwf',1,rjunk,0,' ')
         call iosys('write real "hf m702 1e energy" to rwf',
     $               1,eonel,0,' ')
c         if (abs((eonel-rjunk)/rjunk).gt.1.0d-06) then
c           write (iout,214) eonel,rjunk
c214        format (/,' calculated one-electron energy:',f14.8,
c    #           /,'   one-electron energy from scf:',f14.8)
c           call lnkerr('error in one-electron energy !!!')
c        end if
      end if
c
c
      return
      end

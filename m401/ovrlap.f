*deck @(#)ovrlap.f	5.1  11/6/94
      subroutine ovrlap(ci,cj,exi,exj,z,conti,contj,s,ptprmi,ptprmj,
     #                  noprmi,noprmj,nocnti,nocntj,ptcnti,ptcntj,
     #                  nati,natj,nprmi,nprmj,maxcor,ntypes,nbtype,
     #                  ncnti,ncntj,strti,strtj,nbasi,nbasj,nocart,
     #                  nobf,maxmom,minmom,mintyp,nx,ny,nz,lenxyz)
c***begin prologue     ovrlap
c***date written       230784  yymmdd
c***revision date      850801  yymmdd
c***keywords           overlap, integrals
c***author             saxe, paul, and martin, richard (lanl)
c***source             @(#)ovrlap.f	5.1   11/6/94
c***purpose            forms overlap integral matrix between two generally
c                      contracted gaussian basis sets.
c***description
c     call ovrlap(ci,cj,exi,exj,z,conti,contj,s,ptprmi,ptprmj,
c       noprmi,noprmj,nocnti,nocntj,ptcnti,ptcntj,
c       nati,natj,nprmi,nprmj,maxcor,ntypes,nbtype,
c       ncnti,ncntj,strti,strtj,nbasi,nbasj,nocart,
c       nobf,maxmom,minmom,mintyp,nx,ny,nz,lenxyz)
c***references         (none)
c***routines called    lnkerr(mdutil), stwod(util), fmonel(util), trans1(util),
c                      putrec(m401)
c***end prologue       ovrlap
      implicit integer (a-z)
c
      real*8 z(maxcor),s(nbasi,nbasj)
      real*8 h(21),wt(21)
      real*8 ci(3,nati),exi(nprmi),conti(ncnti)
      real*8 cj(3,natj),exj(nprmj),contj(ncntj)
      integer ptprmi(nati,ntypes),noprmi(nati,ntypes)
      integer nocnti(nati,ntypes)
      integer ptcnti(nati,ntypes),strti(nati,ntypes)
      integer ptprmj(natj,ntypes),noprmj(natj,ntypes)
      integer nocntj(natj,ntypes)
      integer ptcntj(natj,ntypes),strtj(natj,ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nocart(ntypes),nx(lenxyz),ny(lenxyz),nz(lenxyz)
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
c     ----- form the overlap integrals -----
c
      do 9 iatom=1,nati
         do 8 jatom=1,natj
            do 7 itype=1,nbtype
               if (noprmi(iatom,itype).le.0) go to 7
               do 6 jtype=1,nbtype
                  if (noprmj(jatom,jtype).le.0) go to 6
c
                  imax=maxmom(itype)
                  jmax=maxmom(jtype)
                  nprimi=noprmi(iatom,itype)
                  nprimj=noprmj(jatom,jtype)
                  nconti=nocnti(iatom,itype)
                  ncontj=nocntj(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c     ----- allocate core for temporary vectors, etc. -----
c
                  a=1
                  b=a+npint
                  alpha=b+npint
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=1
                  xyz=max(expon+npint,prmint+lenblk*npint)
                  top1=xyz+3*npint*(imax+1)*(jmax+1)
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  if (top1.gt.maxcor.or.top2.gt.maxcor) then
                     call lnkerr('m401: ovrlap..not enough core for s')
                  end if
c
c     ----- form the two-dimensional integrals -----
c
                  call stwod(imax,jmax,nprimi,nprimj,
     #                       exi(ptprmi(iatom,itype)),
     #                       exj(ptprmj(jatom,jtype)),
     #                       z(alpha),ci(1,iatom),cj(1,jatom),
     #                       z(ainv),z(xyza),
     #                       z(expon),z(xyz),npint,
     #                       nbtype,h,wt,mxpts,z(a),z(b))
c
c     ----- form the primitive integrals -----
c
                  call fmonel(z(prmint),z(xyz),npint,lenblk,imax,jmax,
     #                   mintyp(itype),mintyp(itype)+nocart(itype)-1,
     #                   mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     #                   nx,ny,nz)
c
c     ----- transform to contracted functions -----
c
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,conti(ptcnti(iatom,itype)),
     #                        contj(ptcntj(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c     ----- transfer integrals to total array -----
c
                  call putrec(s,z(conint),strti,strtj,iatom,jatom,itype,
     #                        jtype,nconti,ncontj,nbasi,nbasj,lenblk,
     #                        nati,natj,nbtype,nobf)
c
    6          continue
    7       continue
    8    continue
    9 continue
c
c     ----- end timing of this section -----
c
c
      return
      end

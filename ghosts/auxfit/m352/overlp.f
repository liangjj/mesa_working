*deck @(#)overlp.f	1.1  4/28/92
      subroutine overlp(c,ex,z,iz,cont,s,ptprim,noprim,nocont,ptcont,
     #                  nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     #                  start,nbasis,nocart,nobf,maxmom,mintyp,
     #                  nx,ny,nz,minmom,bflabl,ops)
c
c***begin prologue     overlp.f
c***date written       840723  
c***revision date      920427
c      april 27,1992   rlm at lanl
c         modifying oneint.f to do just overlap integrals.
c***keywords           
c***author             saxe, paul and martin, richard (lanl) 
c***source             @(#)overlp.f	1.1   4/28/92
c***purpose            
c***description
c      module to form the overlap integrals over generally contracted
c      gaussian basis sets.
c
c***references
c
c***routines called
c
c***end prologue       overlp.f
      implicit integer (a-z)
c
c     ----- input arrays(unmodified) -----
      character*(*) ops, bflabl(*)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
c
c     ----- scratch working space -----
      real*8 z(maxcor),s(nnp)
      integer iz(*)
c
c     ----- local scratch arrays -----
      real*8 h(21),wt(21)
c
c     ----- input flags -----
c
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
   10 format(5x,'overlap integrals:')
c
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
                     call lnkerr('m352: overlp..not enough core for s')
                  end if
c
c     ----- form the two-dimensional integrals -----
c
                  call sints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       ex,z(alpha),c,z(ainv),z(xyza),
     #                       z(expon),z(xyz),npint,nprim,nat,
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
    6          continue
    7       continue
    8    continue
    9 continue
c
c
      return
      end

*deck @(#)v0.f	1.1  4/25/95
      subroutine v0(c,ex,z,cont,ptprim,noprim,nocont,ptcont,
     $     nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     $     start,nbf,zan,nocart,nobf,maxmom,mintyp,
     $     nx,ny,nz,minmom,ops,v,vnuc,dolp,bflabl,
     $     ngrid,cgrid)
c***begin prologue     v0.f
c***date written       940304    (yymmdd)  
c***revision date      11/6/94
c
c***keywords           
c***author             martin, r.l. and saxe, p.w.
c***source             @(#)v0.f	1.1   4/25/95
c***purpose            forms the electrostatic potential on a grid
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       v0.f
      implicit none
c     --- input variables -----
      integer nat,nbf,nprim,ncont,maxcor,ntypes,nbtype,nnp,ngrid
c     --- input arrays (unmodified) ---
      character*(*) ops
      character*(*) bflabl(*)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 zan(nat)
      real*8 cgrid(3,ngrid)
c     --- input arrays (scratch) ---
      real*8 z(maxcor)
c     --- output arrays ---
      real*8 v(nnp,ngrid),vnuc(ngrid)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout,mxpts
      integer iatom,jatom,itype,jtypmx,jtype
      integer imax,jmax,nprimi,nprimj,nconti,ncontj,lenblk
      integer a,aiaj,b,xyz
      integer alpha,ainv,xyza,expon,xyz0,rysrt,ryswt,xyz1,t1
      integer prmint,conint,nroots
      integer npint,top1,top2,tmp1,len1
      integer at,catom
      real*8 h(21),wt(21)
      real*8 cdif(3),rsq
      real*8 sdot
      real*8 zero,one,thrsh
c
      logical dolp,logkey
c
      parameter (zero=0.0d+00,one=1.0d+00,thrsh=1.0d-08)
      data mxpts /21/
      save h,wt,mxpts
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
 1000 format(1x,'electrostatic potential integrals.')
 1010 format(/5x,'electrostatic potential, grid point',i4)
c
c     --- v will contain the potential matrices
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
                  nroots=(imax+jmax)/2+1
                  nprimi=noprim(iatom,itype)
                  nprimj=noprim(jatom,jtype)
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c                 --- allocate core for temporary vectors, etc. ---
                  a=1
                  b=a+npint
                  aiaj=b+npint
                  xyz0=aiaj+npint
                  xyz1=xyz0+npint*3
                  rysrt=xyz1+npint*3
                  ryswt=rysrt+npint*nroots
                  alpha=ryswt+npint*nroots
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=expon+npint
                  xyz=max(expon+npint,prmint+lenblk*npint)
                  t1=xyz+3*npint*(imax+1)*(jmax+1)*2
                  top1=t1+npint
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
               if (top1.gt.maxcor.or.top2.gt.maxcor) then
                  write(iout,*) 'need',max(top1,top2),' have',maxcor
                  call lnkerr('m619: v0..not enough core')
               end if
c
c              --- form the primitive integrals ---
c              loop over the points we wish to evaluate
               do 25 catom=1,ngrid
                  call v0atr(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       ex,z(alpha),c,z(ainv),z(xyza),
     #                       z(expon),z(xyz),npint,nprim,nat,
     #                       h,wt,mxpts,z(a),z(b),z(aiaj),
     #                       z(xyz0),z(rysrt),z(ryswt),nroots,
     #                       z(prmint),lenblk,zan,z(xyz1),z(t1),
     #                 mintyp(itype),mintyp(itype)+nocart(itype)-1,
     #                 mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     #                 nx,ny,nz,cgrid(1,catom))
c
c                 --- transform to contracted functions ---
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c                 --- transfer integrals to total array ---
                  call put1el(v(1,catom),z(conint),start,iatom,jatom,
     #                        itype,jtype,nconti,ncontj,nnp,lenblk,nat,
     #                        nbtype,nobf)
   25          continue
   26          continue
   27       continue
   28    continue
   29 continue
c
c     --- print the integrals if requested.
      if(logkey(ops,'print=properties=v0',.false.,' ')) then
         write(iout,1000)
         do 33 catom=1,ngrid
            write (iout,1010) catom
            call trtosq(z,v(1,catom),nbf,nnp)
            call wlmat(z,nbf,nbf,bflabl,bflabl)
   33    continue
      endif
c
c     --- compute nuclear contribution and combine
      do 36 catom=1,ngrid
         vnuc(catom)=zero
         do 35 at=1,nat
            call vsub(cdif,c(1,at),cgrid(1,catom),3)
            rsq=sdot(3,cdif,1,cdif,1)
            if(rsq.ge.thrsh) then
               vnuc(catom)=vnuc(catom)+zan(at)/sqrt(rsq)
            endif
   35     continue
   36 continue
      if(dolp) then
c
c        the ecp properties are computed by the simple expedient of 
c        using the reduced nuclear charge.  some research indicates
c        that for "mo better" properties the valence orbitals should be 
c        explicitly orthogonalized to the core orbitals in order 
c        to reintroduce the proper nodal patterns.
      endif
c
c
      return
      end

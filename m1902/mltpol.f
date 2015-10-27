*deck @(#)mltpol.f	5.1  11/6/94
      subroutine mltpol(c,ex,cont,s,ptprim,noprim,nocont,ptcont,
     $                  nat,nprim,ntypes,nbtype,nnp,ncont,
     $                  start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     $                  nx,ny,nz,minmom,lmult,nmat,ops,refops,bflabl)
c***begin prologue     mltpol.f
c***date written       860826  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)mltpol.f	5.1   11/6/94
c***purpose            form electric multipole integrals over contracted
c                      basis sets.
c***description
c   the formation of the multipole integrals is done by manipulation of
c   overlap integrals.   
c
c***references
c
c***routines called
c
c***end prologue       mltpol.f
      implicit none
c     --- input variables -----
      integer nat,nprim,maxcor,ntypes,nbtype,nnp
      integer nbasis,ncont,nmat,lmult
      character*(*) ops
      character*(*) refops
      character*(*) bflabl(*)
      real*8 zan(nat)
c     --- input arrays (unmodified) ---
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
c     --- input arrays (scratch) ---
      real*8 z
c     --- output arrays ---
      real*8 s(nnp,nmat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxpts
      integer first, need, ngot, wpadti, idum
      integer iatom,jatom,itype,jtype,jtypmx
      integer imax,jmax,nprimi,nprimj,nconti,ncontj,npint,lenblk
      integer a,b,alpha,ainv,xyza,expon,prmint
      integer xyz,top1,conint,tmp1,len1,top2
      integer xyzprp,prptyp,prpmom,minprp,prp,powx,powy,powz
      integer at
      real*8 ctest(3)
      real*8 h(21),wt(21)
      logical logkey
      character*6 mult
      character*1 itoc
      character*4 funcnm
      real*8 cdif(3),nucprp
      real*8 jcx,jcy,jcz
      real*8 zero,one
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      data mxpts /21/
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
      save h,wt,mxpts
c
      common/io/inp,iout
      pointer ( pz,z(1) )
      first=0
      maxcor=0
c
 1000 format(1x,'electric multipole integrals.')
 1010 format(/5x,a6,':')
c
c     --- form the overlap integrals.
c         the multipoles will be referred to the origin.
      call rzero(ctest,3)
      do 9 iatom=1,nat
         do 8 jatom=1,iatom
            jcx=c(1,jatom)-ctest(1)
            jcy=c(2,jatom)-ctest(2)
            jcz=c(3,jatom)-ctest(3)
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
c                 --- allocate core for temporary vectors, etc. ---
                  a=1
                  b=a+npint
                  alpha=b+npint
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  xyz=expon+npint
                  prmint=xyz+3*npint*(imax+1)*(jmax+lmult+1)
                  xyzprp=prmint+lenblk*npint
                  top1=xyzprp+3*npint*(imax+1)*(jmax+1)
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
                  need=wpadti( max(top1,top2,nbasis*nbasis) )
                  if(first.eq.0) then
c
c                    this is the first pass through.  get the memory as
c                    nothing has been allocated.  this is all scratch.
c
                     call getmem(need,pz,ngot,'scratch',0)
                     first=1
                     maxcor=ngot
                  else
c                    
c                    this is not the first pass through.  if we need
c                    more memory then we have allocated on the last pass,
c                    we deallocate what we have and then allocate what we need.
c   
                     if(maxcor.lt.need) then
                        call getmem(-ngot,pz,idum,'scratch',idum)
                        call getmem(need,pz,ngot,'scratch',0)
                        maxcor=ngot
                     endif
                  endif
c
c
c                 --- form the two-dimensional integrals ---
c                     have sints form multipole integrals with respect
c                     to center j.
                  call sints(iatom,jatom,imax,jmax+lmult,nprimi,nprimj,
     $                       ptprim(iatom,itype),ptprim(jatom,jtype),
     $                       ex,z(alpha),c,z(ainv),z(xyza),
     $                       z(expon),z(xyz),npint,nprim,nat,
     $                       nbtype,h,wt,mxpts,z(a),z(b))
c
c                 --- convert the overlap integrals
c                     into multipole integrals.
                  prptyp=0
                  do 5 prpmom=0,lmult
                     minprp=mintyp(prpmom+1)
                     do 4 prp=minprp,minprp+nocart(prpmom+1)-1
                        powx=nx(prp)
                        powy=ny(prp)
                        powz=nz(prp)
                        prptyp=prptyp+1
c
c                       --- form the multipole one-dimensional integrals ---
                        call fmmult(z(xyzprp),z(xyz),npint,
     $                              imax,jmax,
     $                              lmult,powx,powy,powz,
     $                              jcx,jcy,jcz)
c
c                       --- form the primitive integrals ---
                        call fmonel(z(prmint),z(xyzprp),npint,lenblk,
     $                              imax,jmax,
     $                              mintyp(itype),
     $                              mintyp(itype)+nocart(itype)-1,
     $                              mintyp(jtype),
     $                              mintyp(jtype)+nocart(jtype)-1,
     $                              nx,ny,nz)
c
c                       --- transform to contracted functions ---
                        call trans1(z(prmint),z(conint),nprimi,nprimj,
     $                              nconti,ncontj,
     $                              cont(ptcont(iatom,itype)),
     $                              cont(ptcont(jatom,jtype)),z(tmp1),
     $                              len1,lenblk,minmom(itype),
     $                              maxmom(itype),minmom(jtype),
     $                              maxmom(jtype),nocart)
c
c                       --- transfer integrals to total array ---
                        call put1el(s(1,prptyp),z(conint),start,iatom,
     $                              jatom,itype,jtype,nconti,ncontj,nnp,
     $                              lenblk,nat,nbtype,nobf)
    4                continue
    5             continue
    6          continue
    7       continue
    8    continue
    9 continue
c
c     --- print the integrals.
      prptyp=0
      do 30 prpmom=0,lmult
         minprp=mintyp(prpmom+1)
         do 20 prp=minprp,minprp+nocart(prpmom+1)-1
            powx=nx(prp)
            powy=ny(prp)
            powz=nz(prp)
            prptyp=prptyp+1
            mult='e'//itoc(prpmom)//funcnm(powx,powy,powz)
c           if(logkey(ops,'print=properties='//mult(1:2),.false.,
c     $        refops)) then
               if(prpmom.eq.0) write(iout,1000)
               write(iout,1010) mult
               call trtosq(z,s(1,prptyp),nbasis,nnp)
               call wlmat(z,nbasis,nbasis,bflabl,bflabl)
c           endif
   20    continue
   30 continue
c
c     --- compute the nuclear contribution and put the multipole 
c         integrals on the read-write file.
      prptyp=0
      do 50 prpmom=0,lmult
         minprp=mintyp(prpmom+1)
         do 40 prp=minprp,minprp+nocart(prpmom+1)-1
            prptyp=prptyp+1
            powx=nx(prp)
            powy=ny(prp)
            powz=nz(prp)
            nucprp=zero
            do 35 at=1,nat
               if(powx.eq.0) then
                  cdif(1)=one
               else
                  cdif(1)=c(1,at)-ctest(1)
               endif
               if(powy.eq.0) then
                  cdif(2)=one
               else
                  cdif(2)=c(2,at)-ctest(2)
               endif
               if(powz.eq.0) then
                  cdif(3)=one
               else
                  cdif(3)=c(3,at)-ctest(3)
               endif
               nucprp=nucprp+zan(at)
     $                *(cdif(1)**powx)*(cdif(2)**powy)*(cdif(3)**powz)
   35       continue
            mult='e'//itoc(prpmom)//funcnm(powx,powy,powz)
            call iosys('create real '//mult//' on rwf',
     $                  4+nnp,0,0,' ')
            call iosys('write real '//mult//' on rwf after rewinding',
     $                  3,ctest,0,' ')
            call iosys('write real '//mult//' on rwf without rewinding',
     $                  1,nucprp,0,' ')
            call iosys('write real '//mult//' on rwf without rewinding',
     $                  nnp,s(1,prptyp),0,' ')
   40    continue
   50 continue
c
c
c     get rid of all the scratch.  we no longer need it.
c
      call getmem(-ngot,pz,idum,'scratch',idum)
      return
      end

*deck @(#)del.f	5.1  11/6/94
      subroutine del(c,ex,cont,s,ptprim,noprim,nocont,ptcont,
     $                  nat,nprim,ntypes,nbtype,nnp,ncont,
     $                  start,nbasis,nocart,nobf,maxmom,mintyp,
     $                  nx,ny,nz,minmom,ops,refops,bflabl)
c***begin prologue     del.f
c***date written       860906  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)del.f	5.1   11/6/94
c***purpose            forms del and del**4 integrals 
c                      over contracted basis. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       del.f
      implicit none
c     --- input variables -----
      integer nat,nprim,maxcor,ntypes,nbtype,nnp
      integer nbasis,ncont
      character*(*) ops
      character*(*) refops
      character*(*) bflabl(*)
c     --- input arrays (unmodified) ---
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
c     --- input arrays (scratch) ---
      real*8 z
c     --- output arrays ---
      real*8 s(nnp,3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxpts
      integer first, need, ngot, wpadti, idum
      integer iatom,jatom,itype,jtype,jtypmx
      integer imax,jmax,nprimi,nprimj,nconti,ncontj,npint,lenblk
      integer a,b,ap2,am2,bp1,bm1,bp2,bm2,alpha,ainv,xyza,expon,prmint
      integer xyz,t1,top1,conint,tmp1,len1,top2,coord
      real*8 zero,ctest(3)
      real*8 h(21),wt(21)
      logical logkey
      character*5 mult,prpnam(3)
c
      parameter (zero=0.0d+00)
c
      data mxpts /21/
      data prpnam/'delx','dely','delz'/
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
      save prpnam,h,wt,mxpts
c
      common/io/inp,iout
      pointer (pz,z(1))
c
 1000 format(1x,'del integrals. ')
 1010 format(/5x,a5)
 1020 format(1x,'del**4 integrals.')
      first=0
      maxcor=0
c
c     --- form the del integrals.
      if(logkey(ops,'properties=del',.false.,refops)) then
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
c                    --- allocate core for temporary vectors, etc. ---
                     a=1
                     b=a+npint
                     bp1=b+npint
                     bm1=bp1+npint
                     alpha=bm1+npint
                     ainv=alpha+2*npint
                     xyza=ainv+npint
                     expon=xyza+3*npint
                     prmint=1
                     xyz=max(expon+npint,prmint+lenblk*npint)
                     t1=xyz+3*npint*(imax+1)*(jmax+1)*2
                     top1=t1+npint
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
C
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
c                    --- form the two-dimensional integrals ---
                     call pone(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                          ptprim(iatom,itype),ptprim(jatom,jtype),
     $                          ex,z(alpha),c,z(ainv),z(xyza),
     $                          z(expon),z(xyz),npint,nprim,nat,
     $                          nbtype,h,wt,mxpts,z(a),z(b),
     $                          z(bp1),z(bm1))
c
                     do 15 coord=1,3
c
c                       --- form the primitive integrals -----
                        call fmdel(z(prmint),z(xyz),npint,lenblk,
     $                              imax,jmax,
     $                              mintyp(itype),
     $                              mintyp(itype)+nocart(itype)-1,
     $                              mintyp(jtype),
     $                              mintyp(jtype)+nocart(jtype)-1,
     $                              nx,ny,nz,coord)
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
                        call put1el(s(1,coord),z(conint),start,iatom,
     $                              jatom,itype,jtype,nconti,ncontj,nnp,
     $                              lenblk,nat,nbtype,nobf)
   15                continue
   16             continue
   17          continue
   18       continue
   19    continue
c
c        --- print the integrals.
         if(logkey(ops,'print=properties=del',.false.,refops)) then
            write(iout,1000)
            do 20 coord=1,3
               write(iout,1010) prpnam(coord)
               call trtosq(z,s(1,coord),nbasis,nnp)
               call wlmat(z,nbasis,nbasis,bflabl,bflabl)
   20       continue
         endif
c
c        --- and store on the read-write file.
         call rzero(ctest,3)
         do 22 coord=1,3
            mult=prpnam(coord)
            call iosys('create real '//mult//' on rwf',4+nnp,0,0,' ')
            call iosys('write real '//mult//' on rwf after rewinding',
     $                 3,ctest,0,' ')
            call iosys('write real '//mult//' on rwf without rewinding',
     $                 1,zero,0,' ')
            call iosys('write real '//mult//' on rwf without rewinding',
     $                 nnp,s(1,coord),0,' ')
   22    continue
      endif
c     get rid of all the scratch.  we no longer need it.
c
      call getmem(-ngot,pz,idum,'scratch',idum)
c
c     --- form the del**4 integrals.

      if(logkey(ops,'properties=mv',.false.,refops)) then
         first=0
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
                     nprimi=noprim(iatom,itype)
                     nprimj=noprim(jatom,jtype)
                     nconti=nocont(iatom,itype)
                     ncontj=nocont(jatom,jtype)
                     npint=nprimi*nprimj
                     lenblk=nocart(itype)*nocart(jtype)
c
c                    --- allocate core for temporary vectors, etc. ---
                     a=1
                     b=a+npint
                     ap2=b+npint
                     am2=ap2+npint
                     bp2=am2+npint
                     bm2=bp2+npint
                     alpha=bm2+npint
                     ainv=alpha+2*npint
                     xyza=ainv+npint
                     expon=xyza+3*npint
                     prmint=1
                     xyz=max(expon+npint,prmint+lenblk*npint)
                     t1=xyz+3*npint*(imax+2)*(jmax+2)*2
                     top1=t1+npint
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
C
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
c                    --- form the two-dimensional integrals ---
                     call pfour(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                          ptprim(iatom,itype),ptprim(jatom,jtype),
     $                          ex,z(alpha),c,z(ainv),z(xyza),
     $                          z(expon),z(xyz),npint,nprim,nat,
     $                          nbtype,h,wt,mxpts,z(a),z(b),
     $                          z(ap2),z(am2),z(bp2),z(bm2))
c
c                    --- form the primitive integrals ---
                     call fmt(z(prmint),z(xyz),npint,lenblk,
     $                        imax,jmax,z(t1),
     $                      mintyp(itype),mintyp(itype)+nocart(itype)-1,
     $                      mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     $                      nx,ny,nz)
c
c                    --- transform to contracted functions ---
                     call trans1(z(prmint),z(conint),nprimi,nprimj,
     $                          nconti,ncontj,cont(ptcont(iatom,itype)),
     $                          cont(ptcont(jatom,jtype)),z(tmp1),len1,
     $                          lenblk,minmom(itype),maxmom(itype),
     $                          minmom(jtype),maxmom(jtype),nocart)
c
c                    --- transfer integrals to total array ---
                     call put1el(s(1,1),z(conint),start,iatom,jatom,
     $                           itype,jtype,
     $                           nconti,ncontj,nnp,lenblk,nat,nbtype,
     $                           nobf)
   26             continue
   27          continue
   28       continue
   29    continue
c
c        --- print the integrals.
         if(logkey(ops,'print=properties=mv',.false.,refops)) then
            write(iout,1020)
            call trtosq(z,s(1,1),nbasis,nnp)
            call wlmat(z,nbasis,nbasis,bflabl,bflabl)
         endif
c
c        --- put the mass-velocity integrals on the read-write file.
         call rzero(ctest,3)
         call iosys('create real "del4 integrals" on rwf',
     $               4+nnp,0,0,' ')
         call iosys('write real "del4 integrals" on rwf',3,ctest,0,' ')
         call iosys('write real "del4 integrals" on rwf '//
     $              'without rewinding',1,zero,0,' ')
         call iosys('write real "del4 integrals" on rwf '//
     $              'without rewinding',nnp,s(1,1),0,' ')
      endif
c
c
c
c     get rid of all the scratch.  we no longer need it.
c
      call getmem(-ngot,pz,idum,'scratch',idum)
      return
      end

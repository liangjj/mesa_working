*deck @(#)efld.f	5.1  11/6/94
      subroutine efld(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     nat,nprim,ntypes,nbtype,nnp,ncont,
     $     start,nbf,zan,nocart,nobf,maxmom,mintyp,
     $     nx,ny,nz,minmom,ops,refops,ds,nderiv,dolp,bflabl,
     $     ngrid,cgrid,lens)
c***begin prologue     efld.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard and saxe,paul (lanl)
c***source             @(#)efld.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       efld.f
      implicit none
c     --- input variables -----
      character*(*) ops
      character*(*) refops
      character*(*) bflabl(*)
      integer nat,nprim,maxcor,ntypes,nbtype,nnp,nbf
      integer ncont,nderiv,ngrid,lens
c     --- input arrays (unmodified) ---
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 zan(nat)
      real*8 cgrid(3,ngrid)
c     --- input arrays (scratch) ---
      real*8 z
c     --- output arrays ---
      real*8 ds(nnp,0:lens-1,ngrid)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      character*4 mult,prpnam(0:9)
      integer inp,iout
      integer mxpts
      integer first, need, ngot, wpadti, idum
      integer i,at
      integer catom,nxyz,nint,iatom,jatom,itype,jtypmx,jtype
      integer imax,jmax,nroots,nprimi,nprimj,nconti,ncontj
      integer npint,lenblk
      integer a,b,aiaj,xyz0,xyz1,rysrt,ryswt,alpha,ainv
      integer xyza,expon,prmint,xyz,t1,top1
      integer conint,tmp1,len1,top2
      integer coord
      real*8 h(21),wt(21)
      real*8 ctest(3),cdif(3),nucprp,rsq,rsqi,ri,r3i,r5i,sdot
      real*8 zero,one,thrsh
      real*8 junk,three
      logical dolp,logkey
c
      parameter (zero=0.0d+00,one=1.0d+00,thrsh=1.0d-08)
      parameter (three=3.0d+00)
c
      common/io/inp,iout
      pointer (pz,z(1))
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
      data prpnam/'v0  ','v1x ','v1y ','v1z ','v2xx','v2yy','v2zz',
     $            'v2xy','v2xz','v2yz'/
c
      save prpnam,h,wt,mxpts
c
 1000 format(/5x,'electrostatic potential:grid point',f12.3)
 1010 format(/5x,'electric field:x, grid point',f12.3)
 1011 format(/5x,'electric field:y, grid point',f12.3)
 1012 format(/5x,'electric field:z, grid point',f12.3)
 1020 format(/5x,'electric field gradient:xx, grid point',f12.3)
 1021 format(/5x,'electric field gradient:yy, grid point',f12.3)
 1022 format(/5x,'electric field gradient:zz, grid point',f12.3)
 1023 format(/5x,'electric field gradient:xy, grid point',f12.3)
 1024 format(/5x,'electric field gradient:xz, grid point',f12.3)
 1025 format(/5x,'electric field gradient:yz, grid point',f12.3)
      first=0
      maxcor=0
c
c     --- zero the space to accumulate the property integrals.
c         ds will contain the potential, then the electric field integrals,
c         then the field gradients.
      call rzero(ds,nnp*lens*ngrid)
c
c     --- loop over the points at which we want the property.
      do 30 catom=1,ngrid
         call vmove(ctest,cgrid(1,catom),3)
c        --- zero space for the gradients and second derivatives
         if (nderiv.eq.2) then
            nxyz=6
            nint=28
         else if (nderiv.eq.1) then
            nxyz=3
            nint=7
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
c                    --- allocate core for temporary vectors, etc. -----
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
                     need=wpadti( max(top1,top2,nbf*nbf) )
c
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
c                    --- form the primitive integrals -----
c
                     call eints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                          ptprim(iatom,itype),ptprim(jatom,jtype),
     $                          ex,z(alpha),c,z(ainv),z(xyza),
     $                          z(expon),z(xyz),npint,nprim,nat,
     $                          nbtype,h,wt,mxpts,z(a),z(b),z(aiaj),
     $                          z(xyz0),z(rysrt),z(ryswt),nroots,
     $                          z(prmint),lenblk,zan,z(xyz1),z(t1),
     $                    mintyp(itype),mintyp(itype)+nocart(itype)-1,
     $                    mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     $                          nx,ny,nz,nderiv,z(conint),
     $                          nconti,ncontj,itype,jtype,
     $                          cont(ptcont(iatom,itype)),
     $                          cont(ptcont(jatom,jtype)),
     $                          z(tmp1),len1,minmom(itype),
     $                          minmom(jtype),nocart,
     $                          ds(1,0,catom),nnp,start,nobf,ctest,
     $                          nxyz,lens)
c
c                     --- transform the potential integrals 
c                         to contracted functions -----
c                         the field and field gradient, if requested,
c                         are transformed in eints.
                     call trans1(z(prmint),z(conint),nprimi,nprimj,
     $                           nconti,ncontj,
     $                           cont(ptcont(iatom,itype)),
     $                           cont(ptcont(jatom,jtype)),
     $                           z(tmp1),len1,
     $                           lenblk,minmom(itype),maxmom(itype),
     $                           minmom(jtype),maxmom(jtype),nocart)
c
c                    --- transfer integrals to total array -----
                     call put1el(ds(1,0,catom),z(conint),start,
     $                           iatom,jatom,itype,jtype,nconti,ncontj,
     $                           nnp,lenblk,nat,nbtype,nobf)
c
   26             continue
   27          continue
   28       continue
   29    continue
   30 continue
c
c
      if(logkey(ops,'print=properties=v0',.false.,refops)) then
         do 33 catom=1,ngrid
            write (iout,1000) (cgrid(i,catom),i=1,3)
            call trtosq(z,ds(1,0,catom),nbf,nnp)
            call wlmat(z,nbf,nbf,bflabl,bflabl)
   33    continue
      endif
      if(logkey(ops,'print=properties=v1',.false.,refops)) then
         do 35 catom=1,ngrid
            do 34 coord=1,3
               if(coord.eq.1) then
                  write (iout,1010) (cgrid(i,catom),i=1,3)
               else if(coord.eq.2) then
                  write (iout,1011) (cgrid(i,catom),i=1,3)
               else if(coord.eq.3) then
                  write (iout,1012) (cgrid(i,catom),i=1,3)
               endif
               call trtosq(z,ds(1,coord,catom),nbf,nnp)
               call wlmat(z,nbf,nbf,bflabl,bflabl)
   34       continue
   35    continue
      endif
      if(logkey(ops,'print=properties=v2',.false.,refops)) then
         do 43 catom=1,ngrid
            do 42 coord=4,9
               if(coord.eq.4) then
                  write (iout,1020) (cgrid(i,catom),i=1,3)
               else if(coord.eq.5) then
                  write (iout,1021) (cgrid(i,catom),i=1,3)
               else if(coord.eq.6) then
                  write (iout,1022) (cgrid(i,catom),i=1,3)
               else if(coord.eq.7) then
                  write (iout,1023) (cgrid(i,catom),i=1,3)
               else if(coord.eq.8) then
                  write (iout,1024) (cgrid(i,catom),i=1,3)
               else if(coord.eq.9) then
                  write (iout,1025) (cgrid(i,catom),i=1,3)
               endif
               call trtosq(z,ds(1,coord,catom),nbf,nnp)
               call wlmat(z,nbf,nbf,bflabl,bflabl)
   42       continue
   43    continue
      endif
c
c     --- compute nuclear contribution and write them to the rwf.
      do 57 coord=0,lens-1
         mult=prpnam(coord)
         call iosys('create real '//mult//' on rwf',ngrid*(4+nnp),
     $              0,0,' ')
         do 56 catom=1,ngrid
            nucprp=zero
            do 55 at=1,nat
               call vsub(cdif,c(1,at),cgrid(1,catom),3)
               rsq=sdot(3,cdif,1,cdif,1)
               if(rsq.ge.thrsh) then
                  rsqi=one/rsq
                  ri=sqrt(rsqi)
                  r3i=ri*rsqi
                  r5i=ri*rsqi*rsqi
                  if(coord.eq.0) then
c                    --- potential
                     nucprp=nucprp+zan(at)*ri
                  else if(coord.lt.4) then
c                    --- field
                     nucprp=nucprp-zan(at)*cdif(coord)*r3i
                  else
c                    field gradient
                     if(coord.eq.4) then
                        junk=cdif(1)*cdif(1)
                        nucprp=nucprp-(zan(at)*(three*junk-rsq)*r5i)
                     else if(coord.eq.5) then
                        junk=cdif(2)*cdif(2)
                        nucprp=nucprp-(zan(at)*(three*junk-rsq)*r5i)
                     else if(coord.eq.6) then
                        junk=cdif(3)*cdif(3)
                        nucprp=nucprp-(zan(at)*(three*junk-rsq)*r5i)
                     else if(coord.eq.7) then
                        junk=cdif(1)*cdif(2)
                        nucprp=nucprp-(zan(at)*three*junk*r5i)
                     else if(coord.eq.8) then
                        junk=cdif(1)*cdif(3)
                        nucprp=nucprp-(zan(at)*three*junk*r5i)
                     else if(coord.eq.9) then
                        junk=cdif(2)*cdif(3)
                        nucprp=nucprp-(zan(at)*three*junk*r5i)
                     endif
                  endif
               endif
   55       continue
            call iosys('write real '//mult//' on rwf without rewinding',
     $                 3,cgrid(1,catom),0,' ')
            call iosys('write real '//mult//' on rwf without rewinding',
     $                 1,nucprp,0,' ')
            call iosys('write real '//mult//' on rwf without rewinding',
     $                 nnp,ds(1,coord,catom),0,' ')
   56    continue
   57 continue
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
c
c     get rid of all the scratch.  we no longer need it.
c
      call getmem(-ngot,pz,idum,'scratch',idum)
      return
      end

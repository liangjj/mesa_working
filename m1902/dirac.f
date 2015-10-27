*deck @(#)dirac.f	5.1  11/6/94
      subroutine dirac(c,ex,cont,s,ptprim,noprim,nocont,ptcont,
     $                  nat,nprim,ntypes,nbtype,nnp,ncont,
     $                  start,nbasis,nocart,nobf,maxmom,mintyp,
     $                  nx,ny,nz,minmom,ops,refops,bflabl,ngrid,
     $                  cgrid)
c***begin prologue     dirac.f
c***date written       860908  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)dirac.f	5.1   11/6/94
c***purpose            forms dirac delta functin integrals over contracted
c                      basis. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       dirac.f
      implicit none
c     --- input variables -----
      integer nat,nprim,maxcor,ntypes,nbtype,nnp
      integer nbasis,ncont,ngrid
      real*8 cgrid(3,ngrid)
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
      integer first, need, ngot, wpadti, idum
      integer iatom,jatom,itype,jtype,jtypmx,catom
      integer imax,jmax,nprimi,nprimj,nconti,ncontj,npint,lenblk
      integer alpha,expon,prmint
      integer xyz,top1,conint,tmp1,len1,top2
      real*8 zero,ctest(3)
      logical logkey
c
      parameter (zero=0.0d+00)
c
      common/io/inp,iout
      pointer (pz,z(1))
c
 1000 format(1x,'fermi contact integrals.')
 1010 format(/5x,'grid point',i3,':')
      first=0
      maxcor=0
c
c     --- form the fermi contact integrals.
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
c                 --- allocate core for temporary vectors, etc. ---
                  alpha=1
                  expon=alpha+2*npint
                  xyz=expon+npint
                  prmint=xyz+3*npint*(imax+1)*(jmax+1)
                  top1=prmint+3*npint*(imax+1)*(jmax+1)
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
                  need=wpadti( max(top1,top2,nbasis*nbasis) )
                  if(first.eq.0) then
c
c                 this is the first pass through.  get the memory as
c                 nothing has been allocated.  this is all scratch.
c
                     call getmem(need,pz,ngot,'scratch',0)
                     first=1
                     maxcor=ngot
                  else
c                 
c                 this is not the first pass through.  if we need
c                  more memory then we have allocated on the last pass,
c                  we deallocate what we have and then allocate what we need.
c 
                     if(maxcor.lt.need) then
                        call getmem(-ngot,pz,idum,'scratch',idum)
                        call getmem(need,pz,ngot,'scratch',0)
                        maxcor=ngot
                     endif
                  endif
c
c
c                 --- loop over the points we wish to test.
                  do 5 catom=1,ngrid
                     ctest(1)=cgrid(1,catom)
                     ctest(2)=cgrid(2,catom)
                     ctest(3)=cgrid(3,catom)
c
c                    --- form the primitive integrals.
                     call fermi(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                          ptprim(iatom,itype),ptprim(jatom,jtype),
     $                          ex,z(alpha),c,
     $                          z(expon),z(xyz),npint,nprim,nat,
     $                          nbtype,ctest)
c
c                    --- form the primitive integrals ---
                     call fmonel(z(prmint),z(xyz),npint,lenblk,
     $                           imax,jmax,
     $                           mintyp(itype),
     $                           mintyp(itype)+nocart(itype)-1,
     $                           mintyp(jtype),
     $                           mintyp(jtype)+nocart(jtype)-1,
     $                           nx,ny,nz)
c
c                    --- transform to contracted functions ---
                     call trans1(z(prmint),z(conint),nprimi,nprimj,
     $                           nconti,ncontj,
     $                           cont(ptcont(iatom,itype)),
     $                           cont(ptcont(jatom,jtype)),z(tmp1),
     $                           len1,lenblk,minmom(itype),
     $                           maxmom(itype),minmom(jtype),
     $                           maxmom(jtype),nocart)
c
c                    --- transfer integrals to total array ---
                     call put1el(s(1,catom),z(conint),start,iatom,
     $                           jatom,itype,jtype,nconti,ncontj,nnp,
     $                           lenblk,nat,nbtype,nobf)
    5             continue
    6          continue
    7       continue
    8    continue
    9 continue
c
c     --- print the integrals.
      if(logkey(ops,'print=properties=fermi',.false.,refops)) then
         write(iout,1000)
         do 30 catom=1,ngrid
            write(iout,1010) catom
            call trtosq(z,s(1,catom),nbasis,nnp)
            call wlmat(z,nbasis,nbasis,bflabl,bflabl)
   30    continue
      endif
c
c     --- put the fermi contact integrals on the read-write file.
      call iosys('create real dirac_delta on rwf',
     $            ngrid*(4+nnp),0,0,' ')
      do 40 catom=1,ngrid
         call iosys('write real dirac_delta on rwf without rewinding ',
     $              3,c(1,catom),0,' ')
         call iosys('write real dirac_delta on rwf without rewinding ',
     $              1,zero,0,' ')
         call iosys('write real dirac_delta on rwf without rewinding ',
     $              nnp,s(1,catom),0,' ')
   40 continue
c
c
c     get rid of all the scratch.  we no longer need it.
c
      call getmem(-ngot,pz,idum,'scratch',idum)
      return
      end

*deck @(#)bigdij.f	5.3 11/28/95
      subroutine bigdij(noprim,nbtype,nocont,ncont,nocart,
     $                  nat,start,nbf,nnp,d,nnshl,dijmax,
     $                  z,left,ndmat)
c***begin prologue     bigdij.f
c***date written       950327   (yymmdd)
c***revision date      11/28/95
c***keywords           direct, j-matrix, density
c***author             r.l. martin(lanl)
c***source             @(#)bigdij.f	5.3 11/28/95
c***purpose            finds largest density matrix element in the set
c                      of integrals from each shell block pair.
c***description
c***references         
c
c***routines called    
c
c***end prologue       bigdij.f
c
      implicit none
c     --- input variables -----
      integer nat,nbtype,ncont,ndmat
      integer nbf,nnp,nnshl,left
c     --- input arrays (unmodified) ---
      integer noprim(nat,nbtype)
      integer nocont(nat,nbtype)
      integer nocart(nbtype)
      integer start(nat,nbtype)
      real*8 d(nnp,ndmat)
c     --- input arrays (scratch) ---
      real*8 z(left)
c     --- output arrays ---
      real*8 dijmax(nnshl)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout,stderr
      integer symcen(2),angmom(2)
      integer iatom,jatom
      integer itype,jtype
      integer jtmax
      integer nfi,nfj
      integer nconti,ncontj,nprimi,nprimj
      integer istrt,jstrt
      integer isamax,ij,ijshl,i,j
c
      logical debug
c
c
      parameter (debug=.false.)
c
      common /io/inp,iout
c
c     --- find the largest element of the density matrix in each
c         ij shell block pair.
c     write(iout,*) 'density',(d(i),i=1,nnp)
      ijshl=0
      do 8000 iatom=1,nat
         symcen(1)=iatom
         do 7000 itype=1,nbtype
            nprimi=noprim(iatom,itype)
            if (nprimi.le.0) go to 7000
            nfi=nocart(itype)
            nconti=nocont(iatom,itype)
            angmom(1)=itype
            istrt=start(symcen(1),angmom(1))
            do 6000 jatom=1,iatom
               symcen(2)=jatom
               if (jatom.eq.iatom) then
                  jtmax=itype
               else
                  jtmax=nbtype
               end if
               do 5000 jtype=1,jtmax
                  nprimj=noprim(jatom,jtype)
                  if (nprimj.le.0) go to 5000
                  nfj=nocart(jtype)
                  ncontj=nocont(jatom,jtype)
                  angmom(2)=jtype
                  jstrt=start(symcen(2),angmom(2))
                  ijshl=ijshl+1
c
c
                  if (left .lt. nconti*ncontj*nfi*nfj*ndmat) then
                     write(iout,*) 'not enough core in bigdij'
                     call lnkerr('bigdij -- not enough core')
                  endif
                  call get1dm(d,nnp,ndmat,z,nconti,ncontj,nfi,nfj,
     $                        istrt,jstrt)
c                 write(iout,*) 'ijshl',
c     $                 ((z(i+(j-1)*nconti*ncontj*nfi*nfj),i=1,
c     $                 nconti*ncontj*nfi*nfj),j=1,ndmat)
                  ij=isamax(nconti*ncontj*nfi*nfj,z,1)
                  dijmax(ijshl)=abs(z(ij))
                  if (ndmat.eq.2) then
                     ij=isamax(nconti*ncontj*nfi*nfj,
     $                    z(1+nconti*ncontj*nfi*nfj),1)
                     dijmax(ijshl)=max(dijmax(ijshl),abs(z(ij)))
                  endif
                  write(iout,*) 'dijmax(',ijshl,')=',dijmax(ijshl)
c
 5000          continue
 6000       continue
 7000    continue
 8000 continue
c
c
      return
      end

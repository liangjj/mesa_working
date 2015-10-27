*deck @(#)minbas.f	5.1  11/6/94
      subroutine minbas(natoms,ex,cf,ptprim,nprim,ptcont,ncont,
     $                 start,pstart,nctype,atomno,atomz,
     $                 ncart,nobf,maxmom,minmom,mintyp,
     $                 nx,ny,nz,nat,nex,ncf,nbf,ntypes,nbtype,lenxyz,
     $                 ediag,corflg)
c***begin prologue     minbas
c***date written       850601  yymmdd
c***revision date      860118  860118
c***keywords           minimal basis, huckel, initial guess
c***author             martin, richard (lanl)
c***source             @(#)minbas.f	5.1   11/6/94
c***purpose            constructs a minimum basis set for a molecule.
c                      the routine also fills vectors containing the
c                      diagonal huckel energies and  orbital types.
c***description
c     call minbas(natoms,ex,cf,ptprim,nprim,ptcont,ncont,
c                 start,pstart,nctype,atomno,atomz,
c                 ncart,nobf,maxmom,minmom,mintyp,
c                 nx,ny,nz,nat,nex,ncf,nbf,ntypes,nbtype,lenxyz,
c                 ediag,corflg)
c
c***revisions
c          870118  rlm  fixed dimensioning problem with maxcf.
c          851015  pws  added pstart to dummy arguments and pstart and
c                       npf to call to fmstrt.
c***references         (none)
c***routines called    replic(util), sto3g(m412), togenc(util), fmstrt(util),
c                      ehuckl(m412), normal(util), traktm(mdutil),
c                      lnkerr(mdutil)
c***end prologue       minbas
      implicit integer(a-z)
      parameter (maxsh=19,maxtyp=25,maxpr=3*19,maxcf=3*105)
      parameter (lensto=3)
      real*8 ex(nex),cf(ncf),atomz(nat),ediag(nbf)
      real*8 aex(maxpr),acf(maxcf)
      integer ptprim(nat,ntypes),nprim(nat,ntypes),ncont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),ncart(ntypes)
      integer atomno(nat),corflg(nbf),pstart(nat,ntypes)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nctype(ntypes),nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer aptco(maxsh),aptpr(maxsh),anprim(maxsh),ancont(maxsh)
      integer typeno(maxsh)
      character*8 typnam(maxtyp)
c
      common/io/inp,iout
c
      data typnam /'s','p','d','f','g','h','sp','pd','spd','df','pdf',
     #             'spdf','fg','dfg','pdfg','spdfg',
     #             'us','up','ud','uf','ug','s-u','p-u','d-u','f-u'/
      save typnam
c
 1000 format(1x,'expected dimensions from sizmb differ from actual needs
     $ in m412(minbas)')
 1010 format(' expected:   nbasj,   nex,   ncf',/12x,3i6)
 1020 format('    found:  nbasis, numex,   numcf',/12x,3i6)
c
c
c     start timing.
c
c     loop through the atoms, getting each basis set.
      numex=0
      numcf=0
      numbf=0
      do 100 atom=1,natoms
         if(atomno(atom).eq.0) goto 100
c
c        have we already picked up this basis?
         do 90 i=1,atom-1
            if(atomno(atom).eq.atomno(i)) then
               call replic(i,atom,ptprim,ptcont,nprim,ncont,atomz,
     $                     natoms,ntypes)
               goto 100
            endif
   90    continue
         call sto3g(atomno(atom),atomz(atom),typeno,
     $               anprim,ancont,aptpr,aptco,aex,
     $               acf,typnam,ncart,
     $               nctype,nshell,maxsh,maxpr,maxcf,ntypes)
c
c       form arrays used in integral programs, creating
c       general contraction scheme as well.
        call togenc(numex,numcf,ntypes,natoms,ptprim,ptcont,
     $              nprim,ncont,ex,cf,maxpr,maxcf,nshell,
     $              typeno,anprim,ancont,aptpr,aptco,
     $              aex,acf,atom,iout,nctype)
c
c
  100 continue
c
c     form 'start', an array pointing to the first basis function
c     for each atom, shell type combination.
      call fmstrt(start,nprim,ncont,natoms,nbtype,nobf,
     $            nbasis,ncbf,ncart,npf,pstart)
c
c     check that the basis size is what was expected.
      if(nbasis.ne.nbf.or.numex.gt.nex.or.numcf.gt.ncf) then
         write(iout,1000)
         write(iout,1010) nbf,nex,ncf
         write(iout,1020) nbasis,numex,numcf
         call lnkerr(' ')
      endif
c
c     retrieve the diagonal huckel elements, and set the core flag.
      call ehuckl(start,nprim,ncont,natoms,nbtype,nbasis,
     $            ncart,minmom,maxmom,atomno,ediag,corflg)
c
c     normalize the basis.
      call normal(natoms,nbtype,ptprim,ptcont,nprim,
     $            ncont,ex,numex,cf,numcf,iout,
     $            nctype,minmom,maxmom)
c
c     stop timing.
c
c
      return
      end

*deck @(#)pm1991.f	5.1  11/6/94
      subroutine pm1991(z,a)
c***begin prologue    pm1991.f
c***date written      860722   (yymmdd)
c***revision date     11/6/94
c
c     28 november 1989 rlm at lanl
c     fixed core allocation bug associated with printing natural orbitals
c     in conjunction with a properties run.
c
c     11 august 1987   pws at lanl
c     fixed core allocation and other junk to handle ions on the end
c     of the list of atoms. nactiv is the number of 'real' atoms,
c     natoms, the total including ions.
c
c     10 august 1987   pws at lanl
c     correcting cases where number of roots requested is greater
c     than the number of configurations.
c
c     10 october 1986  rlm at lanl
c     added evaluation of one-electron properties.
c***keywords          mulliken,population,wave-function analysis,properties
c***author            saxe, paul,    (lanl)
c***source            $W%  11/6/94
c***purpose           to calculate and print the mulliken population
c                     and evaluate the one-electron properties.
c***description       m1991 will form and print mulliken populations
c       given the appropriate density and overlap matrices for either
c       hartree-fock, mcscf or ci wavefunctions.
c
c
c       currently, to get hf populations, one must have run the scf code
c       in this run. hf may be combined with ci or mcscf/casscf as long as
c       both the scf and ci/mcscf have been run; however, it is currently
c       not possible to have ci and mcscf/casscf in the same run since
c       the natural orbitals over-write one another.
c
c***references
c
c***routines called
c***end prologue       pm1991.f
      implicit none
c     --- input variables -----
      integer a(*)
      real*8 z(*)
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxnbf,maxatm
      parameter (maxnbf=2000,maxatm=2000)
c
      integer inp,iout
      integer ntodo,nshell,nbf,natoms,nbtype,nae,nbe
      integer nactiv,i,j,ij,bfstrt,top,s,temp,nnp
      integer mxcore,mintyp,nx,ny,nz,lenxyz,ncont,nobf
      integer c,eig,vset,mo,vector,d,totpop,atomp,atomc,orbpop,bfend
      integer iadtwp,wpadti
      character*16 names(maxatm)
      character*16 bflabl(maxnbf)
      character*16 basnam(maxatm)
      character*8 caltyp
      character*4096 ops
      character*4096 refops
      character chrkey*32
      character*32 namvec
      character*32 occ
      character*16 namchk
      logical logkey
c
      common /io/ inp,iout
c
   20 format (' m1991:')
   24 format(/5x,'wavefunction analysis for ',a32)
c
c     --- load the reference options string.
      refops=' '
c
c     --- recover the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- open the checkpoint file ---
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',1000000,0,0,namchk)
c
c     --- default other options ---
      caltyp='rhf'
      ntodo=1
      nshell=1
c
c     --- get the constants for dividing up core ---
      call iosys('read integer "number of basis functions" from chk',
     $           1,nbf,0,' ')
      call iosys('read integer "number of atoms" from chk',
     $           1,natoms,0,' ')
      call iosys('read integer "number basis types" from chk',
     $           1,nbtype,0,' ')
      call iosys('read integer "number of alpha electrons" from chk',
     $           1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from chk',
     $           1,nbe,0,' ')
      call iosys('length of "power of x" on chk',lenxyz,0,0,' ')
      nnp=nbf*(nbf+1)/2
c
c     --- determine the number of "real" atoms, i.e. those with a
c     basis set.  this is somewhat sleazy, and should be changed someday.
c     it was put here by r.l. martin,june87, in order to bypass a great
c     deal of printing in population matrices associated with point charges.
      call iosys('read character "atomic basis name" from chk',
     $           -1,0,0,basnam)
      nactiv=0
      do 10 i=1,natoms
         if(basnam(i).eq.'nobasis') goto 10
         nactiv=nactiv+1
   10 continue
c
c     --- allocate core ---
      bfstrt=1
      bfend=bfstrt+nactiv
      ncont=bfend+nactiv
      nobf=ncont+natoms*nbtype
      mintyp=nobf+nbtype
      nx=mintyp+nbtype
      ny=nx+lenxyz
      nz=ny+lenxyz
      vset=nz+lenxyz
      s=iadtwp(vset+nbf)
      d=s+nbf**2
      c=d+nbf**2
      eig=c+nbf**2
      totpop=eig+nbf
      atomp=totpop+nbf**2*ntodo
      orbpop=atomp+nactiv**2*ntodo
      atomc=orbpop+nbf*ntodo
      temp=atomc+nactiv*ntodo
      top=wpadti(temp+nnp)
c
      call getscm(top,a,mxcore,'m1991',0)
c
c     --- read in basis set information ---
      call iosys('read integer "pointer to cartesians" from chk',
     $            nbtype,a(mintyp),0,' ')
      call iosys('read integer "power of x" from chk',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from chk',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from chk',-1,a(nz),0,' ')
      call iosys('read integer "number of contraction coefficients" '//
     $            'from chk', natoms*nbtype,a(ncont),0,' ')
      call iosys('read integer "number of pure functions" from chk',
     $            nbtype,a(nobf),0,' ')
      call iosys('read character "basis function labels" from chk',
     $          -1,0,0,bflabl)
c
c     --- read in the overlap integrals ---
      call iosys('read real "overlap integrals" from chk',nnp,
     $            z(temp),0,' ')
      call trtosq(z(s),z(temp),nbf,nnp)
c
c     --- read in the z-matrix names ---
      call iosys('read character "z-names w/o dummies" from chk',
     $            -1,0,0,names)
c
c     --- read in the vectors,eigenvalues, and occupations ---
c         which vectors are we to examine?
      namvec=chrkey(ops,'vectors',' ', ' ')
      if (namvec(1:10).eq.'"no vector') then
         occ='"no occ'//namvec(11:32)
      else if (namvec(1:14).eq.'"mcscf vector"') then
         occ='"mcscf orbital energies"'
      else
         occ='"orbital energies"'
      end if
      call iosys('read real '//namvec//' from chk',
     $           nbf**2,z(c),0,' ')
      call iosys('read real '//occ//' from chk',
     $           nbf,z(eig),0,' ')
c     --- which orbitals are we to analyze?
      call izero(a(vset),nbf)
      call intarr(ops,'orbitals',a(vset),nbf,' ')
c
      write(iout,20)
      write(iout,24) namvec
c
c     --- loop over the orbitals to be examined ---
      do 100 mo=1,nbf
         if(a(vset+mo-1).ne.0) then
            vector=a(vset+mo-1)
            write(iout,6) vector
 6          format (/5x,'orbital:',i5)
c           --- move the vector into the temporary array.
            call scopy(nbf,z(c+(vector-1)*nbf),1,z(temp),1)
c           --- generate a density matrix.
            ij=0
            do 90 i=1,nbf
               do 80 j=1,nbf
                  z(d+ij)=z(temp+i-1)*z(temp+j-1)
                  ij=ij+1
   80          continue
   90       continue
c
c           --- print the vector if wished ---
            if (logkey(ops,'print=vector',.false.,' ')) then
               call wvec(z(temp),z(eig+vector-1),nbf,1,bflabl,' ')
            end if
c
c           --- form the population matrices ---
            call mulpop(nbf,natoms,z(s),z(d),z(totpop),z(atomp),
     $               z(atomc),z(orbpop),a(ncont),a(nobf),nbtype,
     $               a(bfstrt),a(bfend),names,ntodo,a(mintyp),
     $               a(nx),a(ny),a(nz),caltyp,ops,bflabl,refops,
     $               nactiv)
         end if
  100 continue
c
c     --- close the checkpoint file
      call iosys('close chk',namchk,0,0,' ')
c
c     --- terminate this code and start the next ---
      call chainx(0)
c
c
      return
      end

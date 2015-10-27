*deck @(#)pm611.f	5.3 5/12/95
      subroutine pm611(z,a)
c***begin prologue     pm611.f
c***date written       930611  (yymmdd)
c***revision date      5/12/95
c   july 1, 1994       rlm at lanl
c      modifying to do ylm component analysis of the density.
c   march 13, 1994     rlm at lanl
c      fixing bug associated with number of lebedev points for l=13
c***keywords           m611, link 611, density-functional, dft, quadrature
c***author             martin, richard (lanl)
c***source             @(#)pm611.f	5.3 5/12/95
c***purpose            evaluates and analyzes the density on a grid
c***description
c
c
c***references         
c                      A.D. Becke, J. Chem. Phys. 88, 2547 (1988)
c
c***routines called
c***end prologue       pm611.f
c
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a(*)
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer nbf,nnp,multip,natoms,nae,nbe
      integer nprim,ncont,ntypes,nbtype
      integer ptprim,noprim,ptcont,nocont,start
      integer nocart,nobf,minmom,maxmom,mintyp
      integer nx,ny,nz
      integer cont,ex,top
      integer inp,iout
      integer ndmat
      integer ian,lenxyz
      integer c
      integer nradial,order,lmax
      integer d,dtot
      integer maxnbf,maxatm,mxcont,maxl,bigl
      integer canget
      integer iadtwp,wpadti,itobyt
      integer imaxl,charge,size
      integer intkey
      integer itch
      integer atom
      integer minesz
      real*8 toang
      
      parameter (maxnbf=2000,maxatm=2000)
c
      character*8 prtflg
      character*4096 ops
      character*16 names(maxatm)
c
      real*8 fpkey
      real*8 zero,two
      parameter (zero=0.0d+00,two=2.0d+00)
      real*8 dencut,dmcut,defcut,toosmall
      parameter (defcut=1.0d-16,toosmall=1.0d-50)
      real*8 cchg,ochg
c
      logical timeit,debug
      logical prnt,adjust,logkey
      parameter (timeit=.true.)
      parameter (debug=.false.)
c
      common /io/     inp,iout
c
      data prnt/.true./
      save prnt
c
 1000 format(1x,'m611:')
 1010 format(5x,'memory available(bytes)',12x,i9)
 1016 format(5x,'grid integration information:')
 1017 format(8x,'density-matrix cutoff:           ',1pe8.1)
 1018 format(8x,'density cutoff:                  ',1pe8.1)
 1019 format(8x,'voronoi cells are size adjusted')
 1021 format(8x,'maximum l:                       ',i3,
     $      /8x,'nradial:                         ',i3,
     $      /8x,'spline-order:                    ',i3)
 1230 format(5x,'becke charges   and  "radius":')
 1231 format(8x,a8,2x,f8.5,f8.5)
 1232 format(8x,'total   ',2x,f8.5)
 1240 format(16x,'      closed      open     total  "radius"')
 1241 format(8x,a8,2x,4f10.5)
 1242 format(8x,'total   ',2x,3f10.5)
c
c     --- get max core available in integers
      call getscm(0,z,canget,'m611:',0)
c
c     ----- recover the options string -----
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
c
c     ----- check on options and set defaults -----
c     --- quadrature options ---
      nradial=intkey(ops,'charges=nradial',51,' ')
      lmax=intkey(ops,'charges=lmax',11,' ')
      order=intkey(ops,'charges=spline-order',3,' ')
      order=order+1
      dmcut=fpkey(ops,'scf=denmat-cutoff',defcut,' ')
      dencut=fpkey(ops,'scf=density-cutoff',toosmall,' ')
      minesz=intkey(ops,'scf=minesz',100,' ')
      adjust=logkey(ops,'scf=adjustcell',.true.,' ')
c
c     ----- has printing been turned off externally? -----
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
      if(prnt) then
         write(iout,1000)
         write(iout,1010) canget*itobyt(1)
         write(iout,1016)
         write(iout,1021) lmax,nradial-1,order-1
         write(iout,1017) dmcut
         write(iout,1018) dencut
         if(adjust) then
            write(iout,1019) 
         endif
      endif
c
c     ----- get some basic information -----
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
      call iosys('read integer "spin multiplicity" from rwf',
     $           1,multip,0,' ')
      call iosys('read integer "number of alpha electrons" from rwf',
     $           1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $           1,nbe,0,' ')
c
      if (nbf.gt.maxnbf) then
         call lnkerr('character core not long enough')
      end if
c
c     ----- get the number of density matrices ---
      call iosys('read integer "number of hf density matrices" '//
     $           'from rwf',1,ndmat,0,' ')
c
c     --- retrieve information regarding basis set 
c         geometry and atomic numbers. first allocate some core.
      ian=1
      c=iadtwp(ian+natoms)
      top=wpadti(c+3*natoms)
      call iosys('read integer "atomic numbers" from rwf',
     $            -1,a(ian),0,' ')
      call iosys('read real coordinates from rwf',-1,z(c),0,' ')
      call basis(natoms,nbf,nprim,ncont,ntypes,nbtype,lenxyz,
     $           mxcont,maxl,
     $           ptprim,noprim,ptcont,nocont,start,nocart,nobf,
     $           minmom,maxmom,mintyp,nx,ny,nz,cont,ex,
     $           top,z,a)
c
c     --- allocate core for density matrix and properties
      bigl=max(maxl,2)
      imaxl=top
      d=iadtwp(imaxl+natoms)
      dtot=d+ndmat*nnp
      charge=dtot+nnp
      size=charge+ndmat*natoms
      top=wpadti(size+natoms)
      itch=iadtwp(top)
c
c     --- read the density matrices and let's go
c         first one is the closed shell piece, second one open.
      call iosys('read real "hf density matrix" from rwf',
     $            ndmat*nnp,z(d),0,' ')
c
c     --- form the total density
c      call rzero(z(dtot),nnp)
c      call saxpy(nnp,two,z(d),1,z(dtot),1)
c      if(ndmat.eq.2) then
c         call vadd(z(dtot),z(dtot),z(d+nnp),nnp)
c      endif
c
c
      call charges(z(itch),z(d),nbf,nnp,ndmat,natoms,nradial,lmax,
     $             order,dmcut,dencut,adjust,z(itch),a(ian),
     $             z(c),z(ex),z(cont),a(ptprim),a(noprim),a(nocont),
     $             a(ptcont),mxcont,nprim,ntypes,nbtype,ncont,a(start),
     $             a(nocart),a(nobf),a(maxmom),a(minmom),a(mintyp),
     $             a(nx),a(ny),a(nz),lenxyz,z(charge),z(size),a(imaxl),
     $             bigl,minesz,ops)
c
c     --- report the charge analysis
      call iosys('read character "z-names w/o dummies" from rwf',
     $           -1,0,0,names)
      cchg=zero
      ochg=zero
      if(ndmat.eq.1) then
         write(iout,1230)
         do 10 atom=1,natoms
            write(iout,1231) names(atom),z(charge+atom-1),
     $                       z(size+atom-1)*toang
            cchg=cchg+z(charge+atom-1)
   10    continue
         write(iout,1232) cchg
      else
         write(iout,1230)
         write(iout,1240)
         do 20 atom=1,natoms
            write(iout,1241) names(atom),z(charge+atom-1),
     $                       z(charge+natoms+atom-1),
     $                       z(charge+atom-1)+z(charge+natoms+atom-1),
     $                       z(size+atom-1)*toang
            cchg=cchg+z(charge+atom-1)
            ochg=ochg+z(charge+natoms+atom-1)
   20    continue
         write(iout,1242) cchg,ochg,cchg+ochg
      endif
c
c     --- form total charges and send to rwf ---
      if(ndmat.eq.2) then
         call vadd(z(charge),z(charge),z(charge+natoms),natoms)
      endif
      call iosys('write real "becke charges" to rwf',
     $            natoms,z(charge),0,' ')
      call iosys('write real "becke sizes" to rwf',
     $            natoms,z(size),0,' ')
c
c     --- and exit gracefully ---
      call chainx(0)
c
c
      return
      end

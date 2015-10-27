*deck @(#)genbas.f	5.3  5/12/95
      subroutine genbas(atsymb,atomno,atomz,basnam,
     $                  atname,natoms,namdat,inp,iout,bflabl,ops)
c***begin prologue     genbas.f
c***date written       850101  
c***revision date      5/12/95      
c   june 2, 1993       rlm at lanl
c      adding a new basis function type: pdfgh, used by ecp
c      routines for doing second derivatives of f functions.
c
c***keywords           
c***author             saxe, paul (lanl 
c***source             @(#)genbas.f	5.3   5/12/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       genbas.f
      implicit none
c     --- input variables -----
      integer inp,iout,natoms
c     --- input arrays (unmodified) ---
      character*(*) atsymb(0:*),basnam(natoms),atname(natoms),ops
      character*(*) bflabl(*),namdat
      real*8 atomz(natoms)
c     --- input arrays (scratch) ---
      integer a
      real*8 z, t
      integer ngot(2), idum, top
c     --- output arrays ---
      integer atomno(natoms)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ntypes,nbtype
      integer amxsh,amxpr,amxcf,mxpr,mxcf,nindex
      integer shlblks 
c
      parameter (ntypes=26,nbtype=17)
      parameter (amxsh=100,amxpr=300,amxcf=500)
c
c     --- note that nindex is used in some of the integral routines for
c     scratch space. if this parameter is changed here in order to 
c     add a new basis function type, you should search for its
c     definition in other routines (particularly in m312,m323,m712
c     and library/intgrl).
c
      parameter (mxpr=5000,mxcf=15000,nindex=56)
c
      integer nctype(ntypes),ncart(ntypes),nobf(ntypes)
      integer maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(nindex),ny(nindex),nz(nindex)
      integer typeno,ancont,anprim,aptco,aptpr,aex
      integer iadtwp,acf,wpadti,wptoin,nprim,momatm
      integer ncont,pstart,start,ex,cf
      integer i,j,iostat,numex,numcf,atom,nshell
      integer charge,intkey,ne,multip,nae,nbe,bf,type,cont
      integer bfptr,powx,powy,powz,nbasis,ncbf,npf,maxgrp,maxprm
      integer ptprim,ptcont,numdat,nunp
      integer mxcont,maxl,minl,maxblk,d1len,d2len
      integer iatom,itype,lptype,d1type,d2type,dminl,d2minl
      integer nactiv,dmaxl,d2maxl,ind
      integer mxcore
      character*8 typnam(ntypes),shname(amxsh)
      character funcnm*8
      real*8 zero,half
      logical dolp,logkey
      pointer (p,z(1)), (p,a(1))
      pointer (pt,t(1))
c
      data typnam /'s','p','d','f','g','h','sp','pd','spd','df','pdf',
     $             'spdf','fg','dfg','pdfg','spdfg','pdfgh',
     $             'us','up','ud','uf','ug','s-u','p-u','d-u','f-u'/
      data nctype / 1,  1,  1,  1,  1,  1,  2,   2,   3,    2,   3,
     $              4,     2,   3,    4,     5,      5,
     $              2,   2,   2,   2,   2,   2,    2,    2,    2/
      data ncart  / 1,  3,  6, 10, 15, 21,  4,   9,  10,    16,  19,
     $              20,    25,  31,   34,    35,     55,
     $              0,   0,   0,   0,   0,   0,    0,    0,    0/
c     note that this is identical to ncart, except for the composite
c     shells.  we use it in the code as if it were the same as nocart,
c     but it should probably be used if projecting out the s component
c     of the 6 cartesian d functions, etc.
      data nobf   / 1,  3,  6, 10, 15, 21,  4,   8,   9,    12,  15,
     $              16,    16,  21,   24,    25,     35,
     $              0,   0,   0,   0,   0,   0,    0,    0,    0/
      data minmom / 0,  1,  2,  3,  4,  5,  0,   1,   0,    2,   1,
     $              0,     3,   2,    1,     0,      1,
     $              0,   0,   0,   0,   0,   0,    0,    0,    0/
      data maxmom / 0,  1,  2,  3,  4,  5,  1,   2,   2,    3,   3,
     $              3,     4,   4,    4,     4,      5,
     $              0,   1,   2,   3,   4,   0,    1,    2,    3/
      data mintyp / 1,  2,  5, 11, 21, 36,  1,   2,   1,    5,   2,
     $              1,     11,  5,    2,     1,      2,
     $              0,   0,   0,   0,   0,   0,    0,    0,    0/
cc      data nx /0,  1,0,0,  2,1,0,1,0,0,  3,2,1,0,2,1,0,1,0,0,
cc     $         4,3,2,1,0,3,2,1,0,2,1,0,1,0,0/
cc      data ny /0,  0,1,0,  0,1,2,0,1,0,  0,1,2,3,0,1,2,0,1,0,
cc     $         0,1,2,3,4,0,1,2,3,0,1,2,0,1,0/
cc      data nz /0,  0,0,1,  0,0,0,1,1,2,  0,0,0,0,1,1,1,2,2,3,
cc     $         0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/
      data nx /0,  1,0,0,  2,0,0,1,1,0,  3,0,0,1,2,2,1,0,0,1,
     $         4,3,2,1,0,3,2,1,0,2,1,0,1,0,0,
     $         5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0/
      data ny /0,  0,1,0,  0,2,0,1,0,1,  0,3,0,2,1,0,0,1,2,1,
     $         0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,
     $         0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0/
      data nz /0,  0,0,1,  0,0,2,0,1,1,  0,0,3,0,0,1,2,2,1,1,
     $         0,0,0,0,0,1,1,1,1,2,2,2,3,3,4,
     $         0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5/
c
      data zero/0.0d+00/, half/0.5d+00/
      save typnam,nctype,ncart,nobf,minmom,maxmom,mintyp,nx,ny,nz
c
 1000 format(' molecular parameters:')
 1010 format(5x,'charge:                ',i5,
     $      /5x,'multiplicity:          ',i5,
     $      /5x,'number alpha electrons:',i5,
     $      /5x,'number beta electrons: ',i5)
 1020 format(1x,'the number of electrons,',i4,', is incompatible with'
     $         //' the multiplicity,',i2,'.')
c
c     --- read in atomic numbers, charges and basis sets ---
      call iosys('read integer "atomic numbers" from rwf',
     $           -1,atomno,0,' ')
      call iosys('read real "nuclear charges" from rwf',-1,atomz,0,' ')
      call iosys('read character "atomic basis name" from rwf',
     $           -1,0,0,basnam)
      call iosys('read character "z-names w/o dummies" from rwf',
     $           -1,0,0,atname)
c
c     --- allocate core ---
      typeno=1
      anprim=typeno+amxsh
      ancont=anprim+amxsh
      aptpr=ancont+amxsh
      aptco=aptpr+amxsh
      aex=iadtwp(aptco+amxsh)
      acf=aex+amxpr
      ptprim=wpadti(acf+amxcf)
      ptcont=ptprim+ntypes*natoms
      nprim=ptcont+ntypes*natoms
      ncont=nprim+ntypes*natoms
      momatm=ncont+ntypes*natoms
      start=momatm+natoms
      pstart=start+ntypes*natoms
      ex=iadtwp(pstart+ntypes*natoms)
      cf=ex+mxpr
      top=wpadti(cf+mxcf)
      call iosys('read integer mxcore from rwf',1,mxcore,0,' ')
c
      if (top.gt.mxcore) then
         write (iout,1) top,mxcore
         call lnkerr('m102 core size problem')
      else
         call getmem(top,p,ngot(1),'genbas',0)
      end if
c
c     --- check for illegal case of atoms with same basis and
c         different charge
      do 81 i=1,natoms
         do 80 j=1,i-1
            if(basnam(i).eq.'nobasis'.or.basnam(j).eq.'nobasis') then
            else
               if (atomno(i).eq.atomno(j)
     $             .and.basnam(i).eq.basnam(j)
     $             .and.atomz(i).ne.atomz(j))
     $                call lnkerr('illegal atoms '//
     $                         'with same basis but differing charge')
            endif
   80   continue
   81 continue
c
c     --- open up the data file, giving it a unit number ---
      numdat=1
      open (unit=numdat,file=namdat,access='sequential',
     $      form='formatted',iostat=iostat,status='old')
c
c     --- if the data file does not exist, we'll sink or swim
c         on the input
      if (iostat.ne.0) numdat=inp
c
c     --- loop through the atoms, getting each basis set ---
      numex=0
      numcf=0
      do 100 atom=1,natoms
         if(basnam(atom).eq.'nobasis') goto 100
c
c        --- check if already done this atom and basis ---
         do 90 i=1,atom-1
            if (atomno(atom).eq.atomno(i).and.
     $          basnam(atom).eq.basnam(i)) go to 95
   90    continue
         call getbas(atsymb(atomno(atom)),atomno(atom),atomz(atom),
     $               basnam(atom),a(typeno),
     $               a(anprim),a(ancont),a(aptpr),a(aptco),z(aex),
     $               z(acf),shname,nshell,amxsh,amxpr,amxcf,typnam,
     $               nctype,ntypes,numdat,inp,iout)
c
c        --- form arrays used in integral programs, creating
c            general contraction scheme also.
         call togenc(numex,numcf,ntypes,natoms,a(ptprim),a(ptcont),
     $               a(nprim),a(ncont),z(ex),z(cf),mxpr,mxcf,nshell,
     $               a(typeno),a(anprim),a(ancont),a(aptpr),a(aptco),
     $               z(aex),z(acf),atom,iout,nctype)
         go to 100
c
   95    continue
c
         call replic(i,atom,a(ptprim),a(ptcont),a(nprim),a(ncont),
     $               atomz,natoms,ntypes)
c
  100 continue
c
c     --- we now have the atomic charges, possibly adjusted for 
c     effective core potentials.  read the net electronic charge 
c     (default=0)
      charge=intkey(ops,'q',0,' ')
c
c     --- compute the number of electrons from the nuclear charges and
c         the net electronic charge.  avoid being fooled by ghosts, 
c        dummies,or point charges
      ne=0
      do 120 i=1,natoms
         if(abs(atomz(i)).gt.zero.and.basnam(i).ne.'nobasis')
     $                       ne=ne+int(atomz(i)+half)
  120 continue
      ne=ne-charge
c
c     --- read the spin multiplicity; default to the low-spin state..
      multip=intkey(ops,'2s+1',mod(ne,2)+1,' ')
c
c     --- check that the molecular charge and multiplicity are consistent.
c         a system with an even number of electrons cannot have an even
c         multiplicity.
      if(mod(ne,2).eq.mod(multip,2)) then
         write(iout,1020) ne,multip
         call lnkerr(' ')
      endif
c
c     --- store the charge and multiplicity.
      call iosys('write integer charge to rwf',1,charge,0,' ')
      call iosys('write integer "spin multiplicity" to rwf',
     $            1,multip,0,' ')
c
c     --- compute and store the number of alpha electrons,
c         the number of beta electrons, and the total.
      nunp=multip-1
      nae=(ne+nunp)/2
      nbe=(ne-nunp)/2
      write(iout,1000)
      write(iout,1010) charge,multip,nae,nbe
      call iosys('write integer "number of electrons" to rwf',
     $            1,ne,0,' ')
      call iosys('write integer "number of alpha electrons" to rwf',
     $            1,nae,0,' ')
      call iosys('write integer "number of beta electrons" to rwf',
     $            1,nbe,0,' ')
c
c     --- form an array containing basis function labels.
      bf=1
      do 160 atom=1,natoms
         bflabl(bf)(1:8)=atname(atom)
         do 150 type=1,nbtype
            do 140 cont=1,a(ncont+(type-1)*natoms+atom-1)
               bfptr=mintyp(type)-1
               do 130 i=1,nobf(type)
                  bfptr=bfptr+1
                  powx=nx(bfptr)
                  powy=ny(bfptr)
                  powz=nz(bfptr)
                  bflabl(bf)(9:)=funcnm(powx,powy,powz)
                  if(bflabl(bf)(9:9).eq.' ') bflabl(bf)(9:)='s'
                  bf=bf+1
  130          continue
  140       continue
  150    continue
  160 continue
      bf=bf-1
c
c
c     --- form 'start', an array pointing to the first basis
c         function for each atom, shell type combination ---
      call fmstrt(a(start),a(nprim),a(ncont),natoms,nbtype,nobf,
     $            nbasis,ncbf,ncart,npf,a(pstart))
c
c     --- find the maximum angular momentum on each atom ---
      call moment(a(nprim),a(momatm),maxmom,natoms,nbtype,maxgrp)
c
      call iosys('write integer "maximum momentum group" to rwf',
     $            1,maxgrp,0,' ')
      call iosys('write integer "atom maximum momentum" to rwf',
     $             natoms,a(momatm),0,' ')
c
c     --- save the number of basis functions and cartesian bf. ---
      write(iout,109)
  109 format(1x,'basis dimensions:')
      write (iout,110) npf,ncbf,nbasis
  110 format(5x,'$ primitive cartesians ',i5,
     $     /,5x,'$ contracted cartesians',i5,
     $     /,5x,'$ basis functions      ',i5)
      call iosys('write integer "number of basis functions" to rwf',
     $            1,nbasis,0,' ')
      call iosys('write integer "number of cartesian basis '//
     $           'functions" to rwf',1,ncbf,0,' ')
      call iosys('write integer "number of primitive functions" '//
     $           'to rwf',1,npf,0,' ')
c
c     --- normalize the basis ---
      call normal(natoms,nbtype,a(ptprim),a(ptcont),a(nprim),
     $                  a(ncont),z(ex),numex,z(cf),numcf,iout,
     $                  nctype,minmom,maxmom)
c
c     --- print out basis information if desired ---
      if(logkey(ops,'print=basis',.false.,' ')) then
         call basprt(natoms,atsymb,atomz,basnam,typnam,nbtype,a(ptprim),
     $            a(ptcont),a(nprim),a(ncont),z(ex),z(cf),numex,numcf,
     $            a(start),iout,atomno,nctype,minmom,maxmom)
      end if
c
c     --- print out ecp information if desired ---
      if(logkey(ops,'print=ecp',.false.,' ')) then
         call ecpprt(natoms,atsymb,atomz,basnam,typnam,ntypes,nbtype,
     $            a(ptprim),a(ptcont),a(nprim),a(ncont),z(ex),z(cf),
     $            numex,numcf,a(start),iout,atomno,nctype,minmom,maxmom)
      end if
c
c     --- find the largest parameters governing shell blocks.
c         this information is used in the integral links to allocate core.
      maxprm=0
      mxcont=0
      maxl=0
      minl=99
      maxblk=0
      d1len=0
      d2len=0

      shlblks=0
      do 20 iatom=1,natoms
         do 10 itype=1,nbtype
            ind=(itype-1)*natoms +iatom-1
            if(a(nprim+ind).gt.0) then
               shlblks=shlblks+1
               maxprm=max(maxprm,a(nprim+ind))
               mxcont=max(mxcont,a(ncont+ind))
               if(maxmom(itype).gt.maxl) then
                  maxl=maxmom(itype)
                  minl=minmom(itype)
               endif
               maxblk=max(maxblk,ncart(itype)**2)
            endif
   10    continue
   20 continue
c
c     --- perhaps we're going to do ecp's.
      dolp=.false.
      do 40 iatom=1,natoms
         do 30 lptype=nbtype+1,ntypes
            ind=(lptype-1)*natoms +iatom-1
            if(a(nprim+ind).gt.0) dolp=.true.
   30    continue
   40 continue
      if(dolp) then
         d1type=0
         d2type=0
         dmaxl=maxl+1
         dminl=max(0,minl-1)
         d2maxl=maxl+2
         d2minl=max(0,minl-2)
         do 50 itype=1,nbtype
            if(dmaxl.eq.maxmom(itype).and.dminl.eq.minmom(itype))
     $      d1type=itype
            if(d2maxl.eq.maxmom(itype).and.d2minl.eq.minmom(itype))
     $      d2type=itype
   50    continue
         if(d1type.ne.0) then
            d1len=ncart(d1type)*ncart(d1type)
            d2len=ncart(d2type)*ncart(d2type)
         endif
      endif
c
      call iosys('write integer "number of shell blocks" to rwf',
     $            1,shlblks,0,' ')
      call iosys('write integer maxprm to rwf',1,maxprm,0,' ')
      call iosys('write integer maxcont to rwf',1,mxcont,0,' ')
      call iosys('write integer maxl to rwf',1,maxl,0,' ')
      call iosys('write integer maxblk to rwf',1,maxblk,0,' ')
      call iosys('write integer d1maxblk to rwf',1,d1len,0,' ')
      call iosys('write integer d2maxblk to rwf',1,d2len,0,' ')
      call iosys('write integer dolp to rwf',1,dolp,0,' ')
c
c     --- and write the basis set information to the rwf ---
      call iosys('write character "atomic basis name" to rwf',
     $            16*natoms,0,0,basnam)
      call iosys('write real "nuclear charges" to rwf',
     $            natoms,atomz,0,' ')
      call iosys('write integer "pointer to primitives" to rwf',
     $            natoms*ntypes,a(ptprim),0,' ')
      call iosys('write integer "pointer to contraction '//
     $           'coefficients" to rwf',natoms*ntypes,a(ptcont),0,' ')
      call iosys('write integer "number of primitives" to rwf',
     $            natoms*ntypes,a(nprim),0,' ')
      call iosys('write integer "number of contraction coefficients"'//
     $           ' to rwf',natoms*ntypes,a(ncont),0,' ')
      call iosys('write integer "pointer to first function" to rwf',
     $            natoms*ntypes,a(start),0,' ')
      call iosys('write integer "pointer to first primitive" to rwf',
     $            natoms*ntypes,a(pstart),0,' ')
      call iosys('write real exponents to rwf',numex,z(ex),0,' ')
      call iosys('write real "contraction coefficients" to rwf',
     $            numcf,z(cf),0,' ')
      call iosys('write integer "number basis types" to rwf',
     $            1,nbtype,0,' ')
c
c     --- calculate the number of atoms with basis functions
      call iosys('read character "atomic basis name" from rwf',
     $            -1,0,0,basnam)
      nactiv=0
      do 1069 i=1,natoms
         if(basnam(i).eq.'nobasis') goto 1069
         nactiv=nactiv+1
 1069 continue
      call iosys('write integer '//
     $           '"number of atoms with basis functions"'//
     $           ' to rwf',1,nactiv,0,' ')
c
c     --- write basis function ordering info --
      call iosys('write character "basis function labels" to rwf',
     $            16*bf,0,0,bflabl)
      call iosys('write integer "power of x" to rwf',nindex,nx,0,' ')
      call iosys('write integer "power of y" to rwf',nindex,ny,0,' ')
      call iosys('write integer "power of z" to rwf',nindex,nz,0,' ')
      call iosys('write integer nctype to rwf',ntypes,nctype,0,' ')
      call iosys('write integer "number of cartesians" to rwf',
     $            ntypes,ncart,0,' ')
      call iosys('write integer "number of pure functions" to rwf',
     $            ntypes,nobf,0,' ')
      call iosys('write integer "maximum momentum" to rwf',
     $            ntypes,maxmom,0,' ')
      call iosys('write integer "minimum momentum" to rwf',
     $            ntypes,minmom,0,' ')
      call iosys('write integer "pointer to cartesians" to rwf',
     $            ntypes,mintyp,0,' ')
c
c     --- calculate the primitive to contracted transformation matrix
      top=wptoin(npf*nbasis)
      call iosys('read integer mxcore from rwf',1,mxcore,0,' ')
      if (top.gt.mxcore) then
         write (iout,1) top,mxcore
         call lnkerr('m102 core size problem')
      else
         call getmem(top,pt,ngot(2),'genbas',0)
      end if
c
      call maket(t,npf,nbasis,a(ptcont),natoms,nbtype,a(nprim),
     $           a(ncont),z(cf),numcf,ncart,minmom,maxmom,a(start),
     $           a(pstart))
c
      call iosys('write real t(prim,cont) to rwf',
     $            npf*nbasis,t,0,' ')
      call getmem(-ngot(1),p,idum,'genbas',idum)
      call getmem(-ngot(2),pt,idum,'genbas',idum)
c
c
      return
    1    format (//,' $$$$$ m102 need ',i8,' words of core, have ',
     $           i8,//)
      end

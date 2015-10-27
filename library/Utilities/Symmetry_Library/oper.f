*deck @(#)oper.f	5.1 11/6/94
      subroutine oper(maxap3,grp,subgrp,axis,natoms,maxop,nop,cstd,
     $                ctmp,prtsym,dump,redsym,maxmom,ncart,nbf,
     $                atmchg)
c
c***begin prologue     oper
c***date written       871004   (yymmdd)
c***revision date      910610   (yymmdd)
c
c      10 june   1991  rlm at lanl
c                      passing atmchg array along to permut.
c                      zeroing bftran array.
c***keywords           symmetry operations
c***author             saxe, paul (lanl)
c***source             @(#)oper.f	5.1 11/6/94
c
c***purpose            do all the symmetry setup of transformation
c                      matrices, etc.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       oper
c
c
      implicit integer (a-z)
c
      parameter (maxatm=2000)
      parameter (maxlam=5,maxrep=14)
c
      character*8 lirrep(maxrep), gengrp
      character*12 labop(120)
      character*(*) grp,subgrp,axis
      integer a1,a2,a3,a4,ngot(4)
      real*8 atmchg(natoms)
      integer numso(maxrep)
      integer ncart(0:maxmom)
      integer lambda(maxrep)
      logical dump,redsym,prtsym
      real*8 z1,z2,z3,z4
      real*8 cstd(maxap3,3),ctmp(3,natoms)
c
      common/io/inp,iout
      pointer (p1,z1(1)), (p1,a1(1))
      pointer (p2,z2(1)), (p2,a2(1))
      pointer (p3,z3(1)), (p3,a3(1))
      pointer (p4,z4(1)), (p4,a4(1))
c
c
      nocont=1
      bfstrt=nocont+natoms*(maxmom+1)
      need=bfstrt+natoms*(maxmom+1)
c
c      call getscm(need,a,junk,'oper',0)
      call getmem(need,p1,ngot(1),'oper',0)
c
      call iosys('read integer "number of contraction coefficients" '//
     $     'from rwf',natoms*(maxmom+1),a1(nocont),0,' ')
      call iosys('read integer "pointer to first function" '//
     $     'from rwf',natoms*(maxmom+1),a1(bfstrt),0,' ')
ctemp
c
c     this is a fix for situations where the symmetry program has been 
c     tricked by passing it a number of atoms less than the actual.
c     in particular, in cases where point charges are at the end of
c     the list, the sheer number of them causes these routines to be
c     prohibitively expensive.  we must temporarily modify natoms to the 
c     correct number, rearrange the nocont and bfstrt arrays, and then
c     proceed.
      call iosys('read integer "number of atoms" from rwf',
     $            1,tmpatm,0,' ') 
      if (tmpatm.ne.natoms) then
         tcont=1
         tstrt=tcont+tmpatm*(maxmom+1)
         tneed=tstrt+tmpatm*(maxmom+1)
         call getmem(tneed,p2,ngot(2),'oper',0)
         call iosys('read integer "number of contraction '//
     $              'coefficients" from rwf',tmpatm*(maxmom+1),
     $               a2(tcont),0,' ')
         call iosys('read integer "pointer to first function" '//
     $              'from rwf',tmpatm*(maxmom+1),a2(tstrt),0,' ')
c        rearrange
         k=0
         l=0
         do 120 j=0,maxmom
            do 110 i=1,tmpatm
               k=k+1
               if(i.le.natoms) then
                  l=l+1
                  a1(nocont+l-1)=a2(tcont+k-1) 
                  a1(bfstrt+l-1)=a2(tstrt+k-1)
               end if
  110       continue
  120    continue
         call getmem(-ngot(2),p2,idum,'oper',idum)
      end if
cend
c
c
c     copy the coordinates in the standard orientation to an internal array.
c
      do 90 i=1,natoms
         do 85 j=1,3
            ctmp(j,i)=cstd(i,j)
 85      continue
 90   continue
c
c     ----- determine number of generators, etc. for this group -----
c
      call group(grp,ngen,nirrep,lirrep,gengrp,naxis,axis,nop,lambda,
     #           labop,prtsym,redsym,subgrp)
c
      call iosys('write integer "symmetry: number of generators" '//
     $     'to rwf',1,ngen,0,' ')
      call iosys('write integer "number of irreducible '//
     $     'representations" to rwf',1,nirrep,0,' ')
      call iosys('write integer "number of symmetry operations" to rwf',
     $     1,nop,0,' ')
c
c     dump debugging information?
c
      if(dump) then
         write (iout,1) ngen,nirrep,gengrp,axis,naxis,nop,
     $        (lirrep(i),i=1,nirrep)
 1       format (5x,' number of generators: ',i5,/,
     #        5x,' number of irreps    : ',i5,/,
     #        5x,' generic group       : ',a5,/,
     #        5x,' principal axis      : ',a5,/,
     #        5x,' order of main axis  : ',i5,/,
     #        5x,' number of operations: ',i5,/,
     #        5x,' labels of irreps    : ',14a5)
      endif
c
c     ----- size of matrices for function transformations -----
c
      lnbftr=0
      do 10 angmom=0,maxmom
         lnbftr=lnbftr+ncart(angmom)**2*nop
 10   continue
c
c     ----- size of the generator matrices -----
c
      lengen=1
      do 5 irrep=1,nirrep
         lengen=lengen+lambda(irrep)**2
    5 continue
      lengen=lengen-1
c
      maxfnc=ncart(maxmom)
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
c
c
c     ----- core allocation -----
c
      momatm=1
      nx=momatm+natoms
      ny=nx+lenxyz
      nz=ny+lenxyz
      ptbftr=nz+lenxyz
      atprmt=ptbftr+(maxmom+1)
      genpt=atprmt+natoms*nop
      symat=genpt+nirrep
      labels=symat+natoms
c
      t=iadtwp(labels+nbf)
      bftran=t+3*3*nop
      t1=bftran+lnbftr
      gen=t1+3*natoms
      gamma=gen+lengen*ngen
      char=gamma+lengen*nop
      t2=char+nirrep*nop
      bfchar=t2+ncart(maxmom)*maxlam
      bftrac=bfchar+(maxmom+1)*nop
      p=bftrac+(maxmom+1)*nop
      coeffs=p+(ncart(maxmom)*natoms)**2*maxlam
      temp=coeffs+nbf*nbf
      need=wpadti(temp+nbf)
c
c      call getscm(need,a,junk,'oper',0)
      call getmem(need,p2,ngot(2),'oper',0)
c
      call iosys('read integer "atom maximum momentum" from rwf',
     $            natoms,a2(momatm),0,' ')
      call iosys('read integer "power of x" from rwf',
     $            lenxyz,a2(nx),0,' ')
      call iosys('read integer "power of y" from rwf',
     $            lenxyz,a2(ny),0,' ')
      call iosys('read integer "power of z" from rwf',
     $            lenxyz,a2(nz),0,' ')
c
c     ---- calculate the transformation matrices of the coordinates ---
c
      call tmat(z2(t),nop,gengrp,naxis,axis)
      if(dump) then
         write(iout,*) 'coordinate transformation matrices'
         off=0
         do 130 op=1,nop
            write(iout,*) 'oper:',op
            call matout(z2(t+off),3,3,3,3,iout)
            off=off+9
  130    continue
      endif
c
      call iosys('write real "symmetry: coordinate transformations" '//
     $     'to rwf',3*3*nop,z2(t),0,' ')
c
c     ----- and the transformation of the basis functions (p, d, f, ...
c
      call rzero(z2(bftran),lnbftr)
      call trmat(z2(bftran),z2(t),a2(ptbftr),maxmom,nop,ncart,a2(nx),
     $           a2(ny),a2(nz))
c
      call iosys('write real "symmetry: basis transformations" to rwf',
     $     lnbftr,z2(bftran),0,' ')
      call iosys('write integer "symmetry: transformation pointers "'//
     $     'to rwf',maxmom+1,a2(ptbftr),0,' ')
c
c     ----- and the permutation array for the atoms -----
c
      call permut(ctmp,natoms,z2(t),nop,z2(t1),a2(atprmt),atmchg)
c
      call iosys('write integer "symmetry: atom permutations" to rwf',
     $     natoms*nop,a2(atprmt),0,' ')
c
c     ----- work out the symmetry related sets of atoms -----
c
      call atmset(a2(symat),a2(atprmt),natoms,nop,ns,mcu)
c
      call iosys('write integer "symmetry: atom sets" to rwf',
     $     natoms,a2(symat),0,' ')
c
c     ----- allocate a bit more core -----
c
      nsf=ncart(maxmom)*mcu
c
      ptsc=1
      nsymat=ptsc+(maxmom+1)*ns
      relatm=nsymat+ns
      need=relatm+mcu*ns
c
c      call getscm(need,a,junk,'oper',0)
      call getmem(need,p3,ngot(3),'oper',0)
c
c     ----- and fill in the number of symmetry related atoms, and the
c           list of related atoms
c
      call relate(a2(symat),natoms,a3(nsymat),ns,a3(relatm),mcu)
c
      call iosys('write integer "number of symmetry distinct atoms" '//
     $     'to rwf',1,ns,0,' ')
      call iosys('write integer "maximum number of symmetry related '//
     $     'atoms" to rwf',1,mcu,0,' ')
      call iosys('write integer "number of symmetry related atoms" '//
     $     'to rwf',ns,a3(nsymat),0,' ')
      call iosys('write integer "symmetry related atoms" to rwf',
     $     mcu*ns,a3(relatm),0,' ')
c
c     ----- and generate pointers to the symmetry contraction sets -----
c
      call symset(a3(ptsc),ns,maxmom,a3(nsymat),ncart,lnsc,a3(relatm),
     $     mcu,a2(momatm),natoms,naords,maxsao)
c
      call iosys('write integer "symmetry contraction pointers" '//
     $     'to rwf',(maxmom+1)*ns,a3(ptsc),0,' ')
c
c     ----- generate the representation matrices -----
c
      call genrtr(gengrp,naxis,z2(gen),nirrep,ngen,z2(gamma),nop,lambda,
     #            a2(genpt),lengen,z2(char),z2(t),lirrep)
c
c     ----- dump debug information? -----
c
      if(dump) then
         call sprint(a2(genpt),lambda,lirrep,z2(gen),ngen,nirrep,
     $        lengen,z2(gamma),nop,z2(char),labop,z2(t),z2(bftran),
     $        a2(ptbftr),maxmom,a2(atprmt),natoms,ncart)
      endif
c
c     ----- work out the symmetry-adapted linear combinations -----
c           note that the salc transformation matrix which comes back
c           is orthonormal on the unit metric.
c
      aords=1
      sc=iadtwp(aords+maxsao*naords)
      need=wpadti(sc+lnsc)
c
c      call getscm(need,a,junk,'oper',0)
      call getmem(need,p4,ngot(4),'oper',0)
      call salc(a2(atprmt),a2(symat),z2(bfchar),z2(bftrac),lambda,
     $          z2(char),a2(genpt),z2(gamma),z2(bftran),a2(ptbftr),
     #          z2(p),z2(coeffs),maxfnc,natoms,nbf,nop,maxmom,
     $          nirrep,lengen,ncart,z4(sc),lnsc,z2(t2),a2(labels),
     #          maxlam,dump,ncart,a3(nsymat),ns,a3(relatm),mcu,
     $          a2(momatm),a3(ptsc),a4(aords),maxsao,naords,
     $          numso,a1(nocont),a1(bfstrt),z2(temp))
c
      if(prtsym) then
         write (iout,234) (lirrep(i),i=1,nirrep)
 234     format(5x,'number of salcs :',3x,10a5)
 235     format(5x,'                 ',10i5)
         write (iout,235) (numso(i),i=1,nirrep)
      endif
      if (dump) then
         write(iout,236)
  236    format(5x,'salc transformation matrix:')
         call matout(z2(coeffs),nbf,nbf,nbf,nbf,iout)
      end if
c
      call iosys('write character "group symbol" to rwf',0,0,0,subgrp)
      call iosys('write integer "number of symmetry orbitals" to rwf',
     $     nirrep,numso,0,' ')
      call iosys('write character "labels of irreducible '//
     $     'representations" to rwf',nirrep*len(lirrep(1)),0,0,lirrep)
      call iosys('write integer "degeneracies of irreducible '//
     $     'representations" to rwf',nirrep,lambda,0,' ')
      call iosys('write real "salc transformation matrix" to rwf',
     $     nbf*nbf,z2(coeffs),0,' ')
      call iosys('write integer "number of ao reduction sets" to rwf',
     $     1,naords,0,' ')
      call iosys('write integer "maximum number of related aos" to rwf',
     $     1,maxsao,0,' ')
      call iosys('write integer "ao reduction sets" to rwf',
     $     maxsao*naords,a4(aords),0,' ')
      call iosys('write real "symmetry contraction matrices" to rwf',
     $     lnsc,z4(sc),0,' ')
c
      call getmem(-ngot(1),p1,idum,'oper',idum)
      call getmem(-ngot(2),p2,idum,'oper',idum)
      call getmem(-ngot(3),p3,idum,'oper',idum)
      call getmem(-ngot(4),p4,idum,'oper',idum)
c
c
      return
      end

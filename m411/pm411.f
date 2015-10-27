*deck @(#)pm411.f	5.1  11/6/94
      subroutine pm411(z,a)
c***begin prologue     m411
c***date written       910528   (yymmdd)
c***revision date               (yymmdd)
c***author             schneider, barry    (lanl)
c***source             @(#)pm411.f	5.1   11/6/94
c***purpose            prepares an initial guess for the scf link
c***                   checking for linear dependence.
c***description        the strategy is to use the lowest set of orbitals
c***                   from the guess and to compliment these with a
c***                   linearly independent set constructed from the
c***                   overlap matrix or the salc's. 
c
c***references
c
c***routines called
c
c***end prologue       m411
      implicit integer(a-z)
c
c
      integer a(*)
      real*8 z(1)
      parameter (maxnbf=2000, maxrep=14)
      integer numso(maxrep), lambda(maxrep), numout(maxrep)
      integer index(maxnbf)
      real*8 maxerr, fpkey, tol
      character ops*4096
      character*8 lirrep(maxrep)
      character*16 bflabl(maxnbf)
      character *3 itoc
      logical prnt, logkey, chklin
      logical sym
c
c
      common/io/inp,iout
c
 
c
c
 1000 format (//,5x,'m411:check for linear dependence')
 1001 format (/,5x,'calculation performed using symmetry')
 1002 format (/,5x,'no. of irreps',1x,i3,
     1        //,5x,'no. input orbitals in each irrep',(/,10(i4,1x)))
 1003 format(/,5x,'irrep',1x,a8)
 1004 format (/,5x,'no. frozen orbitals of this irrep',1x,i3)
 1005 format (/,5x,'no. linearly independent orbitals of this irrep',
     1              1x,i3)
c
c
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnt=logkey(ops,'print=m411',.false.,' ')
      chklin=logkey(ops,'m411=check-linear-dependence',.false.,' ')
c
      freeze=0
      if(logkey(ops,'m411=freeze',.false.,' ')) then
         freeze=intkey(ops,'m411=freeze',0,' ')
         call intarr(ops,'m411=index',index,freeze,' ')
      endif      
c
      maxerr=fpkey(ops,'maxerr',1.d-06,' ')
      tol=fpkey(ops,'eigtol',1.d-06,' ')
c
c
      sym=.false.
      if (logkey(ops,'scf=symmetry',.false.,' ')) then
          sym=.true.
      endif
c
c     signal our presence.
         write(iout,1000)
c
c     retrieve the basis set information.
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbasis,0,' ')
      nnp=nbasis*(nbasis+1)/2
      nbsq=nbasis*nbasis
c
c     allocate core
      words=9*nbsq+2*nbasis+nnp
      if (sym) then
          words=words+nbsq+iadtwp(2*nbasis)
      endif
      call iosys ('read integer maxsiz from rwf',1,maxsiz,0,' ')
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m411',0)
c
c     allocate core
      c=1
      eigval=c+nbsq
      eigvec=eigval+nbasis
      s=eigvec+nbsq
      tri=s+nbsq
      eigtmp=tri+nnp      
      eigvt=eigtmp+nbasis
      t1=eigvt+nbsq
      t2=t1+nbsq
      ovab=t2+nbsq
      t3=ovab+nbsq
      t4=t3+nbsq
      if (sym) then
          salc=t4+nbsq            
          vecsym=wpadti(salc+nbsq)
          tmpsym=vecsym+nbasis
      endif
c
c
      if (sym) then
          write (iout,1001)
          call iosys('read integer "number of irreducible '//
     $               'representations" from rwf',1,nirrep,0,' ')
          call iosys('read integer "number of symmetry orbitals" '//
     1                'from rwf',nirrep,numso,0,' ')
          call iosys('read integer "degeneracies of irreducible'//
     $               ' representations" from rwf',nirrep,lambda,0,' ')
          call iosys('read character "labels of irreducible '//
     $               'representations" from rwf',nirrep*len(lirrep(1)),
     $                0,0,lirrep)
          call iosys('read real "salc transformation matrix" from rwf',
     1               -1,z(salc),0,' ')
          write (iout,1002) nirrep, (numso(i),i=1,nirrep)
      endif
c
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z(tri),0,' ')
c
c     ----- square up overlap and put in s -----
c
      call trtosq(z(s),z(tri),nbasis,nnp)       
c
c    ----- orthonormalize basis functions or salcs -----
c
c     read the guess from the rwf.
c
      call iosys('read real "guess vector" from rwf',-1,z(c),0,' ')
      call iosys('read real "guess orbital energies" from rwf',-1,
     $           z(eigval),0,' ')
      if (sym) then
          call iosys ('read integer "guess vector symmetries" from '//
     1                'rwf',nbasis,a(vecsym),0,' ')
      endif
c
c     ----- copy frozen vectors into first slots of z(eigvec) -----
c     ----- then do same with eigval and vecsym -----
      if (freeze.ne.0) then
          call cpfrez(z(c),z(eigval),a(vecsym),z(eigvec),z(eigtmp),
     1                a(tmpsym),index,nbasis,freeze,sym)
      endif
 
c     ----- at this point the frozen orbitals are in the first freeze -----
c     ----- positions in eigvec and the corresponding eigenvalues -----
c     ----- in eigval and symmetries in vecsym -----
c     ----- the remaining orbitals, eigenvalues and symmetry numbers -----
c     ----- appear in c, eigval and vecsym respectively -----
      ioff=0
      joff=freeze*nbasis  
      koff=freeze   
c
      if (sym) then
c
          total=0
          ntot=0
          do 10 irrep=1,nirrep
             if (numso(irrep).ne.0) then
                 write (iout,1003) lirrep(irrep)
c
c     ----- schmidt orthogonalize the salcs to the frozen functions -----
c
c     ----- copy the frozen orbitals of this symmetry into t2 -----
c     ----- and overwrite the index array with their positions ----- 
                 if (freeze.ne.0) then
                     count=0
                     iloc=vecsym
                     jloc=eigvec
                     do 20 nfreze=1,freeze
                        if ( a(iloc).eq.irrep ) then
                             count=count+1
                             index(count)=nfreze
                             call scopy(nbasis,z(jloc),1,z(t2),1)
                        endif        
                        iloc=iloc+1
                        jloc=jloc+nbasis
   20                continue  
                    write (iout,1004) count
                    if (count.ne.0) then
                        call schmab(z(s),z(t2),z(salc+ioff),z(ovab),
     1                              z(t1),nbasis,count,numso(irrep),
     2                              prnt,chklin)
                    endif
                 endif
c	
c     ----- compute schmidt orthogonalized salc*salc overlap matrix -----
c 
c
c     -----  s * salc and put in t1 -----
c
                 call aeqbc(z(t1),1,nbasis,
     1                      z(s),1,nbasis,
     2                      z(salc+ioff),1,nbasis,
     3                      nbasis,nbasis,numso(irrep))
c
c      ----- multiply result above by transpose of salc and put in t2 -----
c
                 call aeqbc(z(t2),1,numso(irrep),
     1                      z(salc+ioff),nbasis,1,
     2                      z(t1),1,nbasis,
     3                      numso(irrep),nbasis,numso(irrep))
 
c
                 nnpi=numso(irrep)*(numso(irrep)+1)/2
c
c     ----- diagonalize salc salc overlap in t2 -----
c     ----- will not request s**-1/2 to be computed so tri is scratch -----
c
                 call linvec(z(t2),z(tri),z(eigvt),z(eigtmp),z(t1),
     1                       z(t3),numso(irrep),numout(irrep),nnpi,tol,
     2                       'no s**-1/2','none',0)
                 write (iout,1005) numout(irrep)
                 ntest=count+numout(irrep)
                 if (ntest.gt.numso(irrep)) then
                     call lnkerr('too many vectors in symmetry-'//
     1                           itoc(irrep))
                 endif
                 ntot=ntot+ntest
c
c     ----- express vectors in original ao basis -----
c
                 call aeqbc(z(t1),1,nbasis,
     1                      z(salc+ioff),1,nbasis,
     2                      z(eigvt),1,numso(irrep),
     3                      nbasis,numso(irrep),numout(irrep))
c
c     ----- put in proper slots in eigvec -----  
c
                 call scopy(nbasis*numout(irrep),z(t1),1,
     $                      z(eigvec+joff),1)
                 call scopy(numout(irrep),z(eigtmp),1,z(eigval+koff),1)
                 call putsym(a(vecsym+koff),numout(irrep),irrep)
c
c     ----- overwrite the salcs of this symmetry with the frozen -----
c     ----- orbitals and the new orthogonal vectors of this symmetry -----
c     -----               just computed -----
c
                 if (count.ne.0) then
                     loff=ioff 
                     do 30 icount=1,count
                        eigloc=eigvec+(index(icount)-1)*nbasis
                        call scopy(nbasis,z(eigloc),1,z(salc+loff),1)
                        loff=loff+nbasis
   30                continue
                 endif
                 call scopy(nbasis*numout(irrep),z(t1),1,z(salc+loff),1)
c
c     ----- increment offsets -----
c
                 ioff=ioff+numso(irrep)*nbasis 
                 joff=joff+numout(irrep)*nbasis
                 koff=koff+numout(irrep)
                 numout(irrep)=count+numout(irrep)
                 total=total+numout(irrep)
             endif
c
   10     continue
          call iosys('write integer "number of symmetry orbitals" '//
     1                'to rwf',nirrep,numout,0,' ')
          call iosys('write real "salc transformation matrix" to rwf',
     1               total*nbasis,z(salc),0,' ')
      else
c    
c     ----- no symmetry. just diagonalize overlap -----
c
          call runit(z(t2),nbasis)
c
          if (freeze.ne.0) then
c
c     ----- schmidt and then form new overlap matrix -----
c
             call schmab(z(s),z(eigvec),z(t2),z(ovab),z(t1),
     1                    nbasis,freeze,nbasis,prnt,chklin)
             call aeqbc(z(t1),1,nbasis,
     1                  z(s),1,nbasis,
     2                  z(t2),1,nbasis,
     3                  nbasis,nbasis,nbasis)
             call aeqbc(z(t4),1,nbasis,
     1                  z(t2),nbasis,1,
     2                  z(t1),1,nbasis,
     3                  nbasis,nbasis,nbasis)
          else
             call scopy(nbsq,z(s),1,z(t4),1)
          endif
c
c     ----- diagonalize the overlap matrix -----
c
          call linvec(z(t4),z(tri),z(eigvt),z(eigtmp),z(t1),z(t3),
     1                nbasis,numout(1),nnp,tol,
     2                'no s**-1/2','none',0)
          ntot=freeze+numout(1)
          if (ntot.gt.nbasis) then
              call lnkerr('too many linearly independent vectors')
          endif
c
c     ------ return to original basis -----
c
          call aeqbc(z(t1),1,nbasis,
     1                  z(t2),1,nbasis,
     2                  z(eigvt),1,nbasis,
     3                  nbasis,nbasis,numout(1))
c
c     ----- put in proper slots in eigvec -----  
c
          call scopy(nbasis*numout(1),z(t1),1,z(eigvec+joff),1)
          call scopy(numout(1),z(eigtmp),1,z(eigval+koff),1)
c
      endif
c
 
      call iosys('write real "guess vector" to rwf',ntot*nbasis,
     1            z(eigvec),0,' ')
      call iosys('write real "guess orbital energies" to rwf',ntot,
     $            z(eigval),0,' ')
      call iosys ('write integer "no. linearly independent vectors" '//
     1             'to rwf',1,ntot,0,' ')
c
c    ----- check overlap matrix -----
c
      if (chklin) then
          call aeqbc(z(t1),1,nbasis,
     1               z(s),1,nbasis,
     2               z(eigvec),1,nbasis,
     3               nbasis,nbasis,ntot)
          call aeqbc(z(t2),1,ntot,
     1               z(eigvec),nbasis,1,
     2               z(t1),1,nbasis,
     3               ntot,nbasis,ntot)
          write(iout,*) 'overlap matrix' 
          call matout(z(t2),ntot,ntot,ntot,ntot,iout)
      endif
      call iosys ('write integer maxsiz to rwf',1,maxsiz,0,' ')
c
c
      call iosys('read character "basis function labels" from rwf',
     $            -1,0,0,bflabl)
      call wvec(z(eigvec),z(eigval),nbasis,ntot,bflabl,' ')
c
c     stop timing routines.
c
c     exit.
      call chainx(0)
c
c
      stop
      end

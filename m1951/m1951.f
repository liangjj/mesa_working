*deck @(#)pm1951.f	5.1  11/6/94
      program m1951
c***begin prologue    pm1951.f
c***date written      860722   (yymmdd)
c***revision date     11/6/94
c
c     8  january  1994  rlm at lanl
c        writing the mulliken charges to the rwf for solvent calculations
c     15 february 1993  rlm at lanl
c        adding a warning when one-electron properties are computed
c        for ecp integrals.
c     28 november 1989 rlm at lanl
c        fixed core allocation bug associated with printing natural orbitals
c        in conjunction with a properties run.
c     11 august 1987   pws at lanl
c        fixed core allocation and other junk to handle ions on the end
c        of the list of atoms. nactiv is the number of 'real' atoms,
c        natoms, the total including ions.
c     10 august 1987   pws at lanl
c        correcting cases where number of roots requested is greater
c        than the number of configurations.
c     10 october 1986  rlm at lanl
c        added evaluation of one-electron properties.
c***keywords          mulliken,population,wave-function analysis,properties
c***author            saxe, paul,    (lanl)
c***source            $W%  11/6/94
c***purpose           to calculate and print the mulliken population
c                     and evaluate the one-electron properties.
c***description
c       m1951 will form and print mulliken populations
c       given the appropriate density and overlap matrices for either
c       hartree-fock, mcscf or ci wavefunctions.  it also evaluates the
c       one-electron properties over the various wavefunctions using 
c       property integrals provided by m1902.
c
c       options currently available are:
c          hf                hartree-fock population
c          ci                ci populations included
c          mcscf or casscf   mcscf populations produced
c          nroots=n          number of ci roots to analyze, defaults to 1
c          fullpop           includes orbital overlap populations
c          nopop             disables all population output
c          no_orb_pop        disables gross orbital populations
c          no_bond_order     disables bond-order and spin density output
c          no_total_charges  disables printing of total charges
c          e0                electric monopole(overlap) analysis.
c          e1                electric dipole analysis.
c          e2                electric quadrupole analysis.
c          e3                electric octupole analysis.
c          e4                electric hexadecapole analysis.
c          v0                electrostatic potential analysis.
c          v1                electric field analysis.
c          v2                electric field gradient
c          fermi             fermi contact analysis (dirac delta).
c          del               momentum analysis.
c          mv                mass-velocity analysis.
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
c***end prologue      pm1951.f
      implicit none
c     --- input variables -----
      integer a
      real*8 z
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxnbf,maxatm
      parameter (maxnbf=2000,maxatm=2000)
c
      integer nwks,nroots,nbf,natoms,nbtype,nae,nbe,lenxyz
      integer nnp,nshell,ntodo,nactiv,bfstrt,bfend,ncont,nobf,mintyp
      integer nx,ny,nz
      integer f,s,d,c,eig,totpop,atomp,orbpop
      integer i,j,atomc,temp,top,zan,ddel,prpscr,prpint,prptmp,top1,top2
      integer ngot,nos,occ,t1,t2,root1,root2
      integer shell,nprint,prnos,prnt,root
      integer intkey,min
      integer iadtwp,wpadti
      integer inp,iout
      integer idum
      character*4096 ops
      character*16 names(maxatm)
      character*16 bflabl(maxnbf)
      character*8 symlabl(maxnbf)
      character*8 caltyp
      character*8 citype
      character*4096 refops
      character chrkey*16
      character*4 itoc
      character*3 answer
      character*2 vn
      character*32 xform
      logical logkey,prprun,hf,ci,mcscf,rohf,trans
      logical dolp,dft,usesym
c
      data trans/.false./, mcscf/.false./, ci/.false./
      save trans,mcscf,ci
c
      common /io/ inp,iout
      pointer (p,z(1)), (p,a(1))
c
    6 format (5x,'scf vector:')
    9 format (5x,'mcscf vector:')
    8 format (5x,'averaged natural orbitals:')
   20 format (' m1951:')
   22 format(5x,'warning --- the accuracy of 1e- ecp properties'
     $        //' is largely unknown')
   24 format(5x,'wavefunction analysis for hartree-fock ',
     $           'orbitals')
   25 format(5x,'wavefunction analysis for kohn-sham ',
     $           'orbitals')
   30 format(/5x,'root number:',i3)
   35 format(5x,'wavefunction analysis for ci natural ',
     $              'orbitals')
   36 format(5x,'wavefunction analysis for mcscf ',
     $              'natural orbitals')
 1000 format(5x,'transition analysis for mcscf roots:',2i4,/)
 1010 format(5x,'transition analysis for ci roots:',2i4,/)
 1020 format(5x,'no analysis will be performed on ',
     $          'unconverged mcscf orbitals.')
c
c     --- set up the reference options list.
      call drum
      call iosys('read character options from rwf',-1,0,0,ops)
      refops='hf=tcscf,ci=(roots),mcscf,casscf,'//
     $     'properties=(e0,e1,e2,e3,e4,v0,v1,v2,fermi,del,mv),'//
     $     'populations=(full,orbital,bond,charge),'//
     $     'transition,dft'
      refops=' '
c
c     --- recover the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- check for other options ---
      hf=logkey(ops,'hf',.false.,refops)
      dft=logkey(ops,'dft',.false.,refops)
c     --- for now, treat the dft case the same as hf.
c         the dft vectors are stored as 'scf vectors' as well.
      if(dft) hf=.true.
c     --- perhaps this is an ecp run.
      call iosys('read integer dolp from rwf',1,dolp,0,' ')
c     --- is symmetry in use?
      usesym=.not.logkey(ops,'sym=off',.false.,' ')
c     --- see whether the last many-electron wfn. was generated by 
c         mcscf,ci, or neither.
      call iosys('does mcorci exist on rwf',0,0,0,answer)
      if (answer.eq.'yes') then
         call iosys('read character mcorci from rwf',0,0,0,vn)
         if(vn.eq.'mc') then
c           was this a quadratic scf run?
            call iosys('read character "mcscf: ci used" from rwf',
     $                 0,0,0,citype)
            if (citype.eq.'hf') then
               hf=.true.
            else
               mcscf=.true.
            end if
         else if (vn.eq.'ci') then
            ci=.true.
         else
            call lnkerr('unknown many-electron wavefunction')
         endif
      endif
      if (ci.or.mcscf) then
         call iosys('read integer nwks from rwf',1,nwks,0,' ')
         nroots=min(nwks,intkey(ops,'ci=roots',1,refops))
c        --- check to see if mcscf step failed.
         call iosys('read character "transformation vector" from rwf',
     $               -1,0,0,xform)
c
         if(xform.eq.'"mcscf failed"')  then
            write(iout,1020)
            mcscf=.false.
            nroots=0
         end if
      end if
c
c     --- are we requested to do properties?
      prprun=.false.
      if(      logkey(ops,'properties=e0',.false.,refops)
     $  .or.   logkey(ops,'properties=e1',.false.,refops)
     $  .or.   logkey(ops,'properties=e2',.false.,refops)
     $  .or.   logkey(ops,'properties=e3',.false.,refops)
     $  .or.   logkey(ops,'properties=e4',.false.,refops)
     $  .or.   logkey(ops,'properties=v0',.false.,refops)
     $  .or.   logkey(ops,'properties=v1',.false.,refops)
     $  .or.   logkey(ops,'properties=v2',.false.,refops)
     $  .or.   logkey(ops,'properties=fermi',.false.,refops)
     $  .or.   logkey(ops,'properties=del',.false.,refops)
     $  .or.   logkey(ops,'properties=mv',.false.,refops))prprun=.true.
c
c     --- get the constants for dividing up core ---
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nbf,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $            1,natoms,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $            1,nbtype,0,' ')
      call iosys('read integer "number of alpha electrons" from rwf',
     $            1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $            1,nbe,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      nnp=nbf*(nbf+1)/2
c
c     --- recover the basis function labels.
      call iosys('read character "basis function labels" from rwf',
     $           -1,0,0,bflabl)
c
c     --- set up hf specific information.
      if (hf) then
         rohf=nae.ne.nbe
         call iosys('read integer "number of shells" from rwf',
     $               1,nshell,0,' ')
      else
         rohf=.false.
         nshell=1
      end if
      if (rohf) then
         ntodo=4
      else
         ntodo=1
      end if
c
c     --- retrieve the number of "real" atoms, i.e. those with a
c         basis set.  this is determined and written in m102.
      call iosys('read integer '//
     $     '"number of atoms with basis functions"'//
     $     'from rwf',1,nactiv,0,' ')
c
c     --- allocate core
      bfstrt=1
      bfend=bfstrt+nactiv
      ncont=bfend+nactiv
      nobf=ncont+natoms*nbtype
      mintyp=nobf+nbtype
      nx=mintyp+nbtype
      ny=nx+lenxyz
      nz=ny+lenxyz
      f=iadtwp(nz+lenxyz)
      s=f+nshell
      d=s+nbf**2
      if (rohf) then
         c=d+nbf*nbf*2
      else
         c=d+nbf**2
      end if
      eig=c+nbf**2
      totpop=eig+nbf
      atomp=totpop+nbf**2*ntodo
      orbpop=atomp+nactiv**2*ntodo
      atomc=orbpop+nbf*ntodo
      temp=atomc+nactiv*ntodo
      top1=wpadti(temp+nnp)
c
c     --- properties core.
      if(prprun) then
         zan=totpop
         ddel=zan+natoms
         prpscr=ddel+nactiv
         prpint=prpscr+nbf*nbf
c        need temporary holding area for electrostatic properties.
         prptmp=prpint+nbf*nbf
         top2=wpadti(prptmp+6*nactiv*10)
         top1=max(top1,top2)
         if(dolp) then
            write(iout,22)
         endif
      endif
c
c     --- natural orbitals / ci printing.
      if (ci.or.mcscf) then
         prnos=top1
         prnt=prnos+nroots
         nos=iadtwp(prnt+nroots)
         occ=nos+nbf**2
         t1=occ+nbf
         t2=t1+nbf**2
         top2=wpadti(t2+nbf**2)
      else
         top2=top1
      end if
      top=max(top1,top2)
c
      call getmem(top,p,ngot,'m1951',0)
c
c     --- read in basis set ordering information.
      call iosys('read integer "pointer to cartesians" from rwf',
     $            nbtype,a(mintyp),0,' ')
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
      call iosys('read integer "number of contraction coefficients" '//
     $            'from rwf', natoms*nbtype,a(ncont),0,' ')
      call iosys('read integer "number of pure functions" from rwf',
     $            nbtype,a(nobf),0,' ')
c
c     --- read in the overlap integrals ---
      call iosys('read real "overlap integrals" from rwf',nnp,
     $            z(temp),0,' ')
      call trtosq(z(s),z(temp),nbf,nnp)
c
c     --- read in the z-matrix names ---
      call iosys('read character "z-names w/o dummies" from rwf',
     $            -1,0,0,names)
c
c     --- ok, here we go ---
      write (iout,20)
c
c     --- hartree-fock cases ---
      if (hf) then
         if(dft) then
            write(iout,25)
         else
            write (iout,24)
         endif
c
c        --- read in the density matrix ---
         call iosys('read real f from rwf',nshell,z(f),0,' ')
         call iosys('read real "scf vector" from rwf',
     $              nbf**2,z(c),0,' ')
         call iosys('read real "orbital energies" from rwf',
     $              nbf,z(eig),0,' ')
         if (rohf) then
            caltyp='rohf'
c
c           --- the scf writes closed-shell and open-shell
c               density matrices in that order
            if (nbe.eq.0) then
c              --- in the case of a 1e- calculation only one density 
c                  matrix has been written.
               call rzero(z(d),nbf*nbf)
               call iosys('read real "hf density matrix" from rwf',
     $                    nnp,z(temp),0,' ')
               call trtosq(z(d+nbf**2),z(temp),nbf,nnp)
            else
               call iosys('read real "hf density matrix" from rwf',
     $                     nnp,z(temp),0,' ')
               call trtosq(z(d),z(temp),nbf,nnp)
               call iosys('read real "hf density matrix" from rwf',
     $                    nnp,z(temp),nnp,0)
               call trtosq(z(d+nbf**2),z(temp),nbf,nnp)
            end if
         else
            caltyp='rhf'
            call iosys('read real "hf density matrix" from rwf',nnp,
     $           z(temp),0,' ')
            call trtosq(z(d),z(temp),nbf,nnp)
            if (logkey(ops,'hf=tcscf',.false.,refops)) then
c
c              --- tcscf ---
               do 355 shell=2,3
                  call iosys('read real "hf density matrix" '//
     $                 'from rwf',nnp,z(temp),nnp*(shell-1),' ')
                  call trtosq(z(totpop),z(temp),nbf,nnp)
                  call saxpy(nbf**2,z(f-1+shell),z(totpop),1,z(d),1)
 355           continue
            end if
         end if
c
c        --- print the vectors if wished ---
         if (logkey(ops,'print=scf=vector',.false.,' ')) then
            write (iout,6)
            if (chrkey(ops,'print=scf=vector',' ',' ').eq.'all') then
               nprint=nbf
            else
               nprint=min(nae+4,nbf)
            end if
c
            if(usesym) then
               call iosys('does "orbital symmetry labels"'
     $                    //' exist on rwf',0,0,0,answer)
               if(answer.eq.'yes') then
                  call iosys('read character "orbital symmetry labels"'
     $                       //' from rwf',-1,0,0,symlabl)
               endif
               call matprt(z(c),nbf,nbf,nbf,nprint,1,2,bflabl,symlabl,
     $                     0,z(eig),.true.)
            else 
               call wvec(z(c),z(eig),nbf,nprint,bflabl,' ')
            endif
         end if
c
c        --- form the population matrices ---
         if(logkey(ops,'properties=populations',.true.,refops)) then
            call mulpop(nbf,natoms,z(s),z(d),z(totpop),z(atomp),
     $                  z(atomc),z(orbpop),a(ncont),a(nobf),nbtype,
     $                  a(bfstrt),a(bfend),names,ntodo,a(mintyp),
     $                  a(nx),a(ny),a(nz),caltyp,ops,bflabl,refops,
     $                  nactiv)
            call iosys('write real "scf mulliken charges" to rwf',
     $                  nactiv,z(atomc),0,' ')
         endif
c
c        --- evaluate the properties.
         if(prprun) 
     $      call prpone(nbf,nactiv,z(prpint),z(d),z(prpscr),a(nobf),
     $                  a(mintyp),nbtype,z(zan),z(ddel),z(prptmp),
     $                  a(nx),a(ny),a(nz),caltyp,ops,trans,refops)
      end if
c
c
c     --- ci and mcscf cases ---
      ntodo=1
      if (ci.or.mcscf) then
c
c        --- work out which roots to print no's ---
         call izero(a(prnos),nroots)
         if (logkey(ops,'print=natural-orbitals',.false.,' ')) then
c
            if (chrkey(ops,'print=natural-orbitals',' ',' ').eq.
     $           'all') then
               nprint=nbf
            else
               nprint=min(nae+4,nbf)
            end if
c
            if (logkey(ops,'print=natural-orbitals=roots',.false.,' '))
     $           then
               call izero(a(prnt),nroots)
               call intarr(ops,'print=natural-orbitals=roots',a(prnt),
     $                     nroots,' ')
               do 530 i=1,nroots
                  j=a(prnt+i-1)
                  if (j.gt.0) then
                     a(prnos+j-1)=j
                  end if
 530           continue
            else
               do 540 i=1,nroots
                  a(prnos+i-1)=i
 540           continue
            end if
         end if
c
         if (ci) then
            write (iout,35)
         else if (mcscf) then
            write (iout,36)
            if (logkey(ops,'print=mcscf=vector',.false.,' ')) then
               write (iout,9)
c
               if (chrkey(ops,'print=mcscf=vector',' ',' ').eq.
     $              'all') then
                  nprint=nbf
               else
                  nprint=min(nae+4,nbf)
               end if
c
               call iosys('read real "mcscf vector" from rwf',
     $              nbf**2,z(nos),0,' ')
               call iosys('read real "mcscf orbital energies" '//
     $              'from rwf',nbf,z(occ),0,' ')
               call wvec(z(nos),z(occ),nbf,nprint,bflabl,' ')
            end if
         end if
c
         caltyp='ci'
c
c        --- averaged natural orbitals ---
         if (logkey(ops,'print=natural-orbitals=averaged',.false.,' '))
     $        then
            call iosys('does "no vector 0" exist on rwf',0,0,0,answer)
            if (answer.eq.'yes') then
               write (iout,8)
               call iosys('read real "no vector 0" from rwf',
     $              nbf**2,z(nos),0,' ')
               call iosys('read real "no occ 0" from rwf',
     $              nbf,z(occ),0,' ')
               call wvec(z(nos),z(occ),nbf,nprint,bflabl,' ')
            end if
         end if
c
c        --- loop over roots of the ci ---
         do 100 root=1,nroots
            write (iout,30) root
c
c           --- print the vectors if wished ---
            if (ci.and.a(prnos+root-1).gt.0) then
               write (iout,7)
 7             format (5x,'natural orbitals:')
               call iosys('read real "no vector '//itoc(root)//
     $              '" from rwf',nbf**2,z(nos),0,' ')
               call iosys('read real "no occ '//itoc(root)//
     $              '" from rwf',nbf,z(occ),0,' ')
               call wvec(z(nos),z(occ),nbf,nprint,bflabl,' ')
            end if
c
c           --- read in the ao density matrix.
            if (mcscf) then
               if(root.eq.1) then
                  call iosys('read real "mcscf ao 1pdm" from rwf',
     $                        nnp,z(temp),0,' ')
                  call trtosq(z(d),z(temp),nbf,nnp)
               else
c                 no information for this case
                  go to 100
               end if
            else if (ci) then
c              ci case
               call iosys('does "ao 1pdm '//itoc(root)
     $                    //'" exist on rwf',0,0,0,answer)
               if (answer.eq.'yes') then
c                 note that ci ao 1pdm is stored as square array.
                  call iosys('read real "ao 1pdm '//itoc(root)//
     $                       '" from rwf',-1,z(d),0,' ')
               else
c                 no information for this case
                  go to 100
               endif
            end if
c
c           --- form the population matrices ---
            if(logkey(ops,'properties=populations',.true.,refops)) then
               call mulpop(nbf,natoms,z(s),z(d),z(totpop),z(atomp),
     $                     z(atomc),z(orbpop),a(ncont),a(nobf),nbtype,
     $                     a(bfstrt),a(bfend),names,ntodo,
     $                     a(mintyp),a(nx),a(ny),a(nz),caltyp,ops,
     $                     bflabl,refops,nactiv)
            endif
c
c           --- evaluate the properties.
            if(prprun)
     $         call prpone(nbf,nactiv,z(prpint),z(d),z(prpscr),a(nobf),
     $                     a(mintyp),nbtype,z(zan),z(ddel),z(prptmp),
     $                     a(nx),a(ny),a(nz),caltyp,ops,trans,refops)
  100    continue
      end if
c
c     --- transition properties requested?
      if((chrkey(ops,'properties=transition',' ',refops).ne.' ').and.
     $     prprun) then
         trans=.true.
c        which roots?
         do 200 root1=1,nroots
            if(logkey(ops,'properties=transition='//itoc(root1),
     $           .false.,refops))
     $         then
               do 190 root2=root1+1,nroots
                  if(logkey(ops,'properties=transition='//itoc(root2),
     $                 .false.,refops)
     $                 .or.logkey(ops,'properties=transition=*',
     $                 .false.,refops))
     $               then
c                    look for transition 1-matrix.
                     call iosys('read "ao t1pdm:'//itoc(root1)
     $                    //itoc(root2)//'" from rwf',-1,z(d),0,' ')
                     if(mcscf) then
                        write(iout,1000) root1,root2
                     else if(ci) then
                        write(iout,1010) root1,root2
                     endif
                     call prpone(nbf,nactiv,z(prpint),z(d),z(prpscr),
     $                        a(nobf),a(mintyp),nbtype,z(zan),z(ddel),
     $                        z(prptmp),a(nx),a(ny),a(nz),caltyp,ops,
     $                        trans,refops)
                  endif
  190          continue
            endif
  200    continue
      endif
      call getmem(-ngot,p,idum,'m1951',idum)
c
c     --- terminate this code and start the next ---
      call chainx(0)
c
c
      stop
      end

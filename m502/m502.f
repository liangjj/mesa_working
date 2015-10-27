*deck @(#)pm502.f	5.1  11/6/94
      program m502
c***begin prologue     m502
c***date written       861108   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m502, link 502, ivo orbitals
c***author             schneider, barry (lanl) and byron lengsfield (llnl)
c***source             m502
c***purpose            to form and diagonalize an ivo type hamiltonian
c***description
c     m502 recognizes the options subtrings:
c     timing                 collect and print timing statistics
c     print                  print ivos
c
c***references         original source from paul saxe link 601 to form
c                      generalized fock matrices and scf hamiltonian
c
c***routines called
c***end prologue       m502
      implicit integer(a-z)
c
      parameter (maxrep=14)
      character *4096 ops
      character *8 prtflg, lirrep(maxrep), gengrp
      character*128 namchk,namint
      real*8  f, aij, bij
      logical prnt
      logical logkey, usesym
c
      real*8 z
      integer a
      integer numso(maxrep), lambda(maxrep)
c
      common /io/ inp, iout
      real*8 fn,fd,an,ad,bn,bd
      common /shell/  noshel(10), orbshel(50,10), nvshel(10),
     1 vshel(200,10), f(10), aij(10,10), bij(10,10),
     2 fn(10),fd(10),an(10),ad(10),bn(10),bd(10)
      data memtri /1000000 /
      pointer (p,z(1)), (p,a(1))
c
    2 format(1x,'m502:improved virtual orbitals')
    3 format(5x,'memory use',18x,i9)
    4 format(5x,'all integrals held in core.')
    5 format(5x,'# integral triangles in core',i4)
    6 format (/,5x,'ivo will use symmetry')
    7 format(5x,'reading the chk file')
    8 format(5x,'reading the rwf file')
    9 format(5x,'reordering molecular orbitals; new homo:',i5)
c
      call drum
c     ----- recover the options string -----
c
      call iosys ('read character options from rwf',-1,0,0,ops)
c
c     ----- print turned off externally? -----
c
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      prnt=prtflg.ne.'minimum'
      if(prnt) then
         write(iout,2)
      endif
c
c     ----- check for use of symmetry -----
      usesym=logkey(ops,'ivo=symmetry',.false.,' ')
      
c     ----- open the integral file -----
c
      call iosys('read character "integral filename" from rwf',0,0,0,
     $            namint)
      call iosys ('open ints as old',0,0,0,namint)
c
c     ----- get the dimensions, etc -----
c
      call iosys ('read integer "number of basis functions" from rwf',1,
     $            num,0,' ')
      nnp=(num+1)*num/2
c
c     ----- go get information on no. hamiltonia , f , aij , bij -----
c     ----- parameters from setup                                -----
c
c
      call setup (occmo,virtmo,num,ndmat,prnt,ops)
c
      ncoul=ndmat
      nexch=ndmat
c
c     ----- allocate core -----
c
      trmtnew=1
      eigsav=trmtnew+num*num
      dmatl=eigsav+num
      tplsv=dmatl+nnp*ndmat
c-----                      jmat overlaid with tr2
      tr2=tplsv+nnp
      sq1=tr2+nnp
      sq2=sq1+num*num
      nwords1=sq2+num*num
c-----                      jmat overlaid with tr2
      jmat=tplsv+nnp
      kmat=jmat+nnp*ncoul
      isq=kmat+nnp*nexch
      fockl=isq
c-----                      dsq overlaid with  fockl
      fockc=fockl+num*num
      eig=fockc+num*num
      work1=eig+num
      work2=work1+num
      t1=work2+num
      t2=t1+num*num*nexch
      values=t2+num*num
c-----                      dsq overlaid with  fockl
      dsq=isq+num*num
      nwords2=values+nnp*nnp
      nwords=max(nwords1,nwords2)
c-----
      if(usesym) then
         salc=nwords
         nwords=nwords+num*num
         call iosys('read integer "number of irreducible'//
     1              ' representations" from rwf',1,nirrep,0,' ')
         call iosys('read integer "number of symmetry orbitals" from'//
     1              ' rwf',nirrep,numso,0,' ')
         call iosys('read integer "degeneracies of irreducible'//
     1              ' representations" from rwf',nirrep,lambda,0,' ')
         call iosys('read character "group symbol" from rwf',0,0,0,
     1               gengrp)
         call grpsym(gengrp,lirrep,maxrep)
      endif   
c
c
c
c     ----- get the core needed -----
c
      wish=wpadti(nwords+nnp*nnp)
      call getmem(wish,p,canget,'m502',1)
      if (wish.le.canget) then
          ntriang=nnp
      else
         left=min(memtri,iadtwp(canget)-nwords-nnp)
         ntriang=left/nnp+1
      end if
      lenbuf=ntriang*nnp
c      call getscm (top,z,maxcor,'m502',1)
c
c      left=iadtwp(maxcor)-values
c      ntriang=min(left/nnp,nnp)
c
      if (prnt) then
         write(iout,3) canget
         if(ntriang.eq.nnp) then
            write(iout,4)
         else
            write(iout,5) ntriang
         endif
      endif
c
      if (ntriang.lt.1) call lnkerr ('not enough core in m502')
c
c     ----- read in transformation vectors and re-arrange -----
c     ----- to correspond to shell structure              -----
c
      if(logkey(ops,'ivo=chk',.false.,' ')) then
         write(iout,7)
         call iosys('read character "checkpoint filename" from rwf',
     $               0,0,0,namchk)
         call iosys('open chk as old',0,0,0,namchk)
         call iosys ('read real "scf vector" from chk',num*num,
     #                z(trmtnew),0,' ')
         call iosys('read real "orbital energies" from chk',
     #               num,z(eigsav),0,' ')
      else if(logkey(ops,'ivo=chkmc',.false.,' ')) then
         write(iout,7)
         call iosys('read character "checkpoint filename" from rwf',
     $               0,0,0,namchk)
         call iosys('open chk as old',0,0,0,namchk)
         call iosys ('read real "mcscf vector" from chk',num*num,
     #                z(trmtnew),0,' ')
         call iosys('read real "mcscf orbital energies" from chk',
     #               num,z(eigsav),0,' ')
      else
         write(iout,8)
         call iosys ('read real "scf vector" from rwf',num*num,
     #                z(trmtnew),0,' ')
         call iosys('read real "orbital energies" from rwf',
     #               num,z(eigsav),0,' ')
      end if
c
      if (usesym) then
          call iosys ('read real "salc transformation matrix" from rwf',
     1                -1,z(salc),0,' ')
      endif
c
      if(logkey(ops,'ivo=homo',.false.,' ')) then
         homo=intkey(ops,'ivo=homo',occmo,' ')
         write(iout,9) homo
         call flip(z(trmtnew),num,homo,occmo)
      else
         homo=0
      end if
c
c     ----- now calculate the density matrices for each shell -----
c
      call dmat (z(dmatl),z(trmtnew),num,nnp,ndmat)
c
c     ----- read in the kinetic and potential integrals -----
c     ----- and then add them together. results are stored -----
c     ----- in triangular form in z(tplsv)                 -----
c
      call onel (z(tplsv),z(tr2),nnp)
c
c     ----- go form the coulomb and exchange matrices in the -----
c               ----- a.o. basis -----
c
      call fock (z(values),z(dmatl),nnp,num,z(jmat),z(kmat),ncoul,nexch
     1           ,ntriang,z(isq),ndmat,z(dsq))
c
c     ----- form the shell hamiltonia and diagonalize -----
c     -----     to obtain the new orbitals            -----
c
      call ham (z(fockc),z(fockl),z(trmtnew),z(tplsv),z(jmat),z(kmat),
     1          z(t1),z(t2),z(eig),z(work1),z(work2),z(eigsav),z(dmatl),
     2          ops,num,occmo,nnp,ndmat,homo,prnt)
c
c     ----- and exit gracefully -----
c
      call chainx (0)
c
      stop
      end

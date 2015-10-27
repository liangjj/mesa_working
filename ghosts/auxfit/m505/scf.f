*deck %W%  %G%
      subroutine scf(s,t,v,h,f,u,eigval,t1,t2,t3,triang,values,
     $               c,smhalf,nbf,nnp,lenbuf,itap44,d,t4,t5,error,
     $               focksv,mxiter,finish,stdiis,nocc,iprint,dlast,
     $               flast,ntriang,bmatrx,nshell,
     $               shlmin,shlmax,alpha,beta,pulay,calc,ncoul,nexch,
     $               ndmat,nfock,extrap,diagnl,energs,diissv,prnt,ops,
     $               mxds,ptr,dsptr,level,level1,cdiis,fcoef,jmat,kmat,
     $               bflabl,usesym,salc,numso,lambda,lirrep,nirrep,
     $               orbsym,occsym,temp,kept,nbfx,nnpx)
c
c  02 august  1992    rlm at lanl
c       discovered that eigenvalues are not ordered properly in the
c       fix dated 16 may, 1991. is now fixed.
c  30 may     1991    rlm at lanl
c       print option changed back to print=m503=vector
c       print=scf=vector will print the scf vector only once in m1951 
c       at the end of an optimization, e.g.  print=m503=vector will print
c       it every time you come through here.
c  16 may     1991    bis at lanl
c       fixed bug in unique when using symmetry.
c       diagonalization now done symmetry by symmetry.
c       also added option to write vectors in order expected for kohn
c       calculations.
c  09 january 1988    bhl at brl
c       print option changed form print=m503=vector to
c       print=scf=vector
c  21 april 1987       pws at lanl
c       passing level-shift parameters 'level' and 'level1' through to
c       qmat.
c  21 february 1987    pws at lanl
c       reducing core requirements by not holding all the diis information
c       for long past iterations. to do this, i am introducing a new
c       option (maxdiis-->'mxds') which is the maximum number of
c       iterations to hold for diis purposes (default is 5). this also
c       adds a couple of helper arrays: 'dsptr' which is 'mxds' long
c       and holds the iteration number of the data in that diis storage
c       bin, and 'ptr', 'mxiter' long, which holds which diis bin a
c       particular iterations data is in.
c  13 may 1987     pws at lanl
c      adding tcscf capabilities. the two ci coefficients are 'c1'
c      and 'c2'.
c
c***module to perform scf iterations
c
c paul saxe                 27 july 1984                   lanl
c
      implicit integer (a-z)
c
      character*(*) ops
      character*16 bflabl(nbf)
      character*16 chrkey
      character*(*) lirrep(nirrep)
      character*(*) calc
      real*8 s(nnp),t(nnp),v(nnp),h(nnp),u(nbf,nbf),eigval(nbf)
      real*8 t1(nbf,nbf),t2(nbf,nbf),t3(nbf,nbf),triang(nnp)
      real*8 values(lenbuf),c(nbf,nbf),smhalf(nnp),bmatrx(*)
      real*8 t4(nbf,nbf),t5(nbf,nbf),error(nbf,nbf,mxds)
      real*8 focksv(nnp,mxds),flast(nnp),dlast(nnp)
      real*8 f(nnp),d(nnp,ndmat),salc(nbf,nbf)
      real*8 alpha(nshell,nshell),beta(nshell,nshell),fcoef(nshell)
      real*8 energs(mxiter),diissv(mxiter)
      real*8 energy,ten,enuc,eonel,etwoel
      real*8 diiser,maxerr,level,level1,cdiis(nbf,nbf)
      real*8 gamma,cicoef(2)
      real*8 jmat(nnp,ncoul),kmat(nnp,nexch)
      real*8 damp(50),xmin
      integer ptr(mxiter),dsptr(mxds)
      integer shlmin(nshell),shlmax(nshell)
      integer numso(nirrep),lambda(nirrep),occsym(nirrep,nshell)
      integer orbsym(nbf),temp(nbf,2),kept(nshell)
      logical prnt,logkey,printj,printk,first,usesym,debug
      data damp/50*1./
c
      common /io/     inp,iout
c
      data first /.true./, debug/.false./
      data ten /10.0d+00/
      data maxerr/1.0d-06/
      save damp,first,debug,ten,maxerr
c
c     ----- note t1 and bmatrx are implicitly equivalenced -----
c
c     ----- set the pointers to diis storage -----
c
      call izero(dsptr,mxds)
      do 700 i=1,mxiter
         ptr(i)=-9999999
  700 continue
c
c
      if(logkey(ops,'scf=damp',.false.,' ')) then
         call fparr(ops,'scf=damp',damp,nshell-1,' ')
         write(iout,699) (damp(i),i=1,nshell-1)
 699     format(/,' damping will be invoked during the',
     $            ' gscf iterations',/,' damping factor per shell ',/,
     $          10(2x,f5.3))
      end if
c
c     ----- set up options for printing j and k matrices -----
c
      printj=logkey(ops,'scf=print=j',.false.,' ')
      printk=logkey(ops,'scf=print=k',.false.,' ')
c
c     ----- print t and v if requested -----
c
      if (logkey(ops,'scf=print=kinetic',.false.,' ')) then
         write (iout,103)
  103    format (//,' the kinetic-energy one-electron integrals',/)
         call print(t,nnp,nbf,iout)
      end if
      if (logkey(ops,'scf=print=potential',.false.,' ')) then
         write (iout,104)
  104    format (//,' the potential-energy one-electron integrals',/)
         call print(v,nnp,nbf,iout)
      end if
c
c     ----- form h -----
c
      call vadd(h,t,v,nnp)
c
c     ----- print h if requested ------
c
      if (logkey(ops,'scf=print=one-electron',.false.,' ')) then
         write (iout,105)
  105    format (//,' the total one-electron integrals ',/)
         call print(h,nnp,nbf,iout)
      end if
c
c     ----- retrieve guess vector -----
c
      call iosys('read real "guess vector" from rwf',-1,c,0,' ')
      call iosys('read real "guess orbital energies" from rwf',
     $           -1,eigval,0,' ')
c
c     ----- normalize -----
c
      call schmdt(c,s,t1,t2,t3,nbf,nbf,nnp,maxerr)
c
c
      if (logkey(ops,'scf=print=guess',.false.,' ')) then
         write (iout,1)
    1    format (/,' initial guess vector ',/)
         call vecout(c,eigval,nbf)
      end if
c
c     ----- get the nuclear repulsion -----
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,enuc,0,' ')
      if(prnt) then
         write(iout,106) enuc
  106    format(5x,'nuclear repulsion energy:',f12.6)
      endif
c
c     ----- zero density and fock matrices for last iteration -----
c
      call rzero(dlast,nnp)
      call rzero(flast,nnp)
c     --------------------------------------------------------
c     --------------------------------------------------------
c
c     ----- initialize the loop for iterations -----
c
      ndiis=0
      iter=0
      mndiis=1
      call rzero(energs,mxiter)
      call rzero(diissv,mxiter)
      if (prnt) then
         if (calc.eq.'gvb') then
            write (iout,5002)
 5002       format(5x,'iter              energy    diis convergence',
     $             '        c(0)       c(1)')
         else
            write(iout,33)
         end if
      end if
c
c..bhl 10/12/89 llnl
      noabort=0
c
c
c     ----- the iteration loop -----
    2 continue
      iter=iter+1
c
c     ----- form the density matrices -----
c
      do 124 i=1,nshell-1
         call gdmat(d(1,i),c,nbf,nnp,shlmin(i),shlmax(i))
  124 continue
c
      if (logkey(ops,'scf=print=density',.false.,' ')) then
         write (iout,70)
   70    format (' density matrices')
         do 772 junk=1,ndmat
            call print(d(1,junk),nnp,nbf,iout)
  772    continue
      end if
c
c     ----- change to delta(d) -----
c
c      call vsub(d,d,dlast,nnp)
c      call vadd(dlast,dlast,d,nnp)
c      if (btest(iprint,19)) then
c         scalar=0.0d+00
c         do 901 iq=1,nnp
c            scalar=max(scalar,abs(d(iq,1)))
c  901    continue
c         write (iout,902) scalar
c  902    format (' largest element of d:',e15.3)
c         write (iout,119)
c  119    format (/,' delta(density-matrix)',//)
c         call print(d,nnp,nbf,iout)
c      end if
c
c     ----- form the j and k matrices -----
c
      xs=1
      x2=xs+nnp*nbfx
      t6=x2+nnpx 
      t7=t6+nbf*nbf*nbfx
      t8=t7+nbfx*nbfx
      t9=t8+nbf*nbf*nbfx
      t10=t9+nbf*nbf*nbfx
      t11=t10+nbf*nbf*ncoul
      top=t11+nbf*nbf*nexch
      if(iter.eq.1) then
         call iosys('read real "expansion overlap integrals" from rwf',
     $               nnp*nbfx,values(xs),0,' ')
         call iosys('read real "auxiliary two-electron integrals"'
     $            //' from rwf',nnpx,values(x2),0,' ')
      endif
      call formjk(values(xs),values(x2),nbfx,nnpx,d,nnp,nbf,
     $            jmat,kmat,ncoul,nexch,ndmat,t4,values(t6),values(t7),
     $            values(t8),values(t9),values(t10),values(t11))
c
c     ----- perhaps print the mmatrices -----
      if (printj) then
         write (iout,80)
   80    format (' j matrices')
         do 773 junk=1,ncoul
            call print(jmat(1,junk),nnp,nbf,iout)
  773    continue
      end if
      if (printk) then
         write (iout,85)
   85    format (' k matrices')
         do 774 junk=1,ncoul
            call print(kmat(1,junk),nnp,nbf,iout)
  774    continue
      end if
c
c     ----- calculate the energy -----
c
      call hfenrg(jmat,kmat,d,nnp,nbf,ncoul,nexch,ndmat,h,enuc,
     $            energy,cicoef,gamma,calc,eonel,etwoel,
     $            nshell,fcoef,alpha,beta)
      energs(iter)=energy
c
c     ----- form a unique, pseudocanonical representation and
c           transform j and k to this basis
c
      call pseud(nbf,nnp,nshell,ncoul,nexch,h,t1,jmat,kmat,
     $     c,t2,t3,eigval,fcoef,alpha,beta,shlmin,shlmax,t4,t5,
     $     ops)
c
c
c     ----- find the gvb coefficients -----
c
      if (calc.eq.'gvb') then
         nv=2
         nnv=nv*(nv+1)/2
         root=1
         call gvbcof(jmat,kmat,t1,nnp,ncoul,nexch,t2,nnv,shlmin,
     $               shlmax,nshell,nv,shlmin(2),cicoef,root,t3,t4,
     $               t3(1,2),t3(1,3))
c
c        ----- redo coefficients for gvb cases in particular -----
c
         fcoef(1)=1.0d+00
         fcoef(2)=cicoef(1)**2
         fcoef(3)=cicoef(2)**2
         alpha(1,1)=2.0d+00
         alpha(2,1)=2.0d+00*fcoef(2)
         alpha(1,2)=alpha(2,1)
         alpha(2,2)=fcoef(2)
         alpha(3,1)=2.0d+00*fcoef(3)
         alpha(1,3)=alpha(3,1)
         alpha(3,2)=0.0d+00
         alpha(2,3)=0.0d+00
         alpha(3,3)=fcoef(3)
         beta(1,1)=-1.0d+00
         beta(2,1)=-fcoef(2)
         beta(1,2)=beta(2,1)
         beta(2,2)=0.0d+00
         beta(3,1)=-fcoef(3)
         beta(1,3)=beta(3,1)
         beta(3,2)=cicoef(1)*cicoef(2)
         beta(2,3)=beta(3,2)
         beta(3,3)=0.0d+00
      end if
c
c     ----- form the q matrices (capital and small) -----
c
      call qmat(jmat,kmat,t1,nnp,ncoul,nexch,f,pulay,calc,shlmin,
     $          shlmax,nshell,nbf,energy,extrap,diagnl,
     $          level,level1,gamma,cicoef,fcoef,alpha,beta,damp)
c
c
      if (logkey(ops,'scf=print=q-matrix',.false.,' ')) then
         write (iout,1013)
 1013    format (/,t10,'q-matrix in the m. o. basis')
         call print(f,nnp,nbf,iout)
      end if
c
c     ----- transform the 'fock' matrix to the ao basis -----
c
      call trtosq(t3,s,nbf,nnp)
      call trtosq(t1,f,nbf,nnp)
      call ebct(t2,t1,c,nbf,nbf,nbf)
      call ebc(t1,c,t2,nbf,nbf,nbf)
      call ebc(t2,t1,t3,nbf,nbf,nbf)
      call ebc(t1,t3,t2,nbf,nbf,nbf)
      call sqtotr(f,t1,nbf,nnp)
c
      if (calc.eq.'gvb') then
         fcoef(1)=1.0d+00
         fcoef(2)=cicoef(1)**2
         fcoef(3)=cicoef(2)**2
      end if
c
c     ----- form the total density matrix -----
c
      call rzero(triang,nnp)
      do 340 shell=1,nshell-1
         do 330 i=1,nnp
            triang(i)=triang(i)+fcoef(shell)*d(i,1)
  330    continue
  340 continue
c
c     ----- perform pulay's diis extrapolation -----
c
      call diis(f,t1,t2,t3,t4,t5,triang,error,
     $          diissv,energs,bmatrx,cicoef,
     $          focksv,c,cdiis,diiser,energy,stdiis,iprint,
     $          nnp,nbf,mxds,mxiter,iter,calc,
     $          ndiis,mndiis,ptr,dsptr,
     $          prnt,first,ops)
c
c     ----- transform to either the  m. o. basis or the salc basis -----
c
      call trtosq(t2,f,nbf,nnp)
c
      if(usesym) then
ctemp        call ebc(t1,t2,salc,nbf,nbf,nbf)
ctemp        call ebtc(t2,salc,t1,nbf,nbf,nbf)
      else
         call ebc(t1,t2,c,nbf,nbf,nbf)
         call ebtc(t2,c,t1,nbf,nbf,nbf)
      endif
c
      call sqtotr(f,t2,nbf,nnp)
c..bhl 4/19/90
c
      if (logkey(ops,'scf=shift',.false.,' ')) then
         maxocc=shlmax(nshell-1)
         write(iout,*)' shift option invoked maxooc ',maxocc
         ix=0
         xmin=f(1)
         do 994 i=1,maxocc
            ix=ix+i
            xmin=max(xmin,f(ix))
 994     continue
         xmin=xmin+.05d+00
         do 995 i=maxocc+1,nbf
            ix=ix+i
            if(f(ix).lt.xmin) then
               f(ix)=xmin+.05d+00
            end if
 995     continue
      end if
c
c..bhl
c
c
c     ----- print the fock-matrices if requested -----
c
      if (logkey(ops,'scf=print=mo_fock_matrix',.false.,' ')) then
         write (iout,120)
  120    format (/,' the fock-matrices in the m.o. basis ',/)
         call print(f,nnp,nbf,iout)
      end if
      if (logkey(ops,'scf=print=diagonal_fock',.false.,' ')) then
         maxocc=min(nbf,shlmax(nshell-1)+5)
         ix=0
         do 991 i=1,maxocc
            ix=ix+i
            triang(i)=f(ix)
 991     continue
         write (iout,992)
 992     format (/,' the fock-matrices in the m.o. basis ',/)
         write(iout,993)(triang(i),i=1,maxocc)
 993     format(5(2x,f10.4))
      end if
      do 191 i=1,nnp
         triang(i)=f(i)
  191 continue
cend
c      if (btest(iprint,16)) then
c         write (iout,116)
c  116    format (//,' f(tilde) ',/)
c         call print(triang,nnp,nbf,iout)
c      end if
c
c     ----- diagonalize f(tilde) -----
c
      if(usesym) then
         call symrsp(triang,u,t2,t3,t4,t5,eigval,temp,kept,nbf,nnp,
     $               nshell,shlmin,shlmax,nirrep,numso,lambda,
     $               orbsym,occsym,salc,c)
      else
         call degrsp(nbf,nnp,triang,t4,1,u,t2,t3)
c
c        ----- form new scf vector = s**(-1/2)(transpose) * u -----
c
         call ebc(t1,c,u,nbf,nbf,nbf)
         call vmove(c,t1,nbf**2)
      endif
c
c
      if (logkey(ops,'scf=print=rotation-matrix',.false.,' ')) then
         write (iout,117)
  117    format (//,' eigenvectors and eigenvalues of f(tilde)',/)
         call vecout(u,t4,nbf)
      end if
cend
c
      if (logkey(ops,'scf=print=viters',.false.,' ')) then
         write (iout,74) iter
 74      format (//,' vector at the end of iteration',i5,/)
         if(usesym) then
            write(iout,*) ' orbital symmetries'
            write(iout,*) (orbsym(i),i=1,nbf)
         end if
         call vecout(c,eigval,nbf)
      end if
c
      if (diiser.lt.ten**(-finish)) go to 90
c
      if (iter.lt.mxiter) go to 2
c     -------------------------------------------------------------------
c     -------------------------------------------------------------------
c     ----- end iteration loop -----
c
      do 32 i=1,iter
         write(iout,31) i,energs(i),diissv(i)
   32 continue
   31 format(5x,i4,2(5x,f15.9))
   33 format(5x,'iter              energy    diis convergence')
      if(.not.logkey(ops,'scf=noabort',.false.,' ')) then
         call lnkerr('convergence not met.')
      end if
      noabort=1
c
   90 continue
c
c     ----- print the number of shells and the alpha and beta matrices --
c
      if (logkey(ops,'scf=print=coefficients',.false.,' ')) then
         write (iout,1701) (iq,shlmin(iq),shlmax(iq),
     $        shlmax(iq)-shlmin(iq)+1,iq=1,nshell)
 1701    format (/,5x,' orbital shell     first   last   number',/,
     $        (t11,i3,t15,i3,t25,i3,t35,i3))
         write (iout,1702)
 1702    format (/,5x,'alpha matrix')
         call matout(alpha,nshell,nshell,nshell,nshell,iout)
         write (iout,703)
 703     format (/,5x,'beta matrix')
         call matout(beta,nshell,nshell,nshell,nshell,iout)
         write (iout,704)
 704     format (/)
      end if
c
c     ----- form the average fock matrix and rotate independent pairs
c           of orbitals under this fock matrix
c
c     ----- and print the final energy and vector if desired -----
c
      if (prnt.and..not.logkey(ops,'print=scf=energy',.false.,' ')) then
         if (calc.eq.'gvb') then
            write (iout,5001) iter,energy,diiser,cicoef(1),cicoef(2)
 5001       format(5x,i4,2(5x,f15.9),2x,2f11.4)
         else
            write(iout,31) iter,energs(iter),diissv(iter)
         end if
      endif
c
      if (calc.eq.'closed') then
         nmat=1
      else if (calc.eq.'open') then
         nmat=2
      else if (calc.eq.'gvb') then
         nmat=3
      end if
c
c..ci opt     if (logkey(ops,'scf=core-fock',.false.,' ')) then
      if (.not.logkey(ops,'scf=gscf',.false.,' ')
     $    .and. (calc.ne.'general')) then
         if(prnt) then
            write (iout,801)
 801        format(5x,'rotating scf orbitals under core fock operator')
         endif
c        rotate the orbitals under the core-fock operator.
c        this means the virtual eigenvalues will be appropriate
c        for the closed shell operator.
         call unique(c,eigval,f,jmat,kmat,values,d,u,t1,t2,h,
     $               nbf,nnp,nshell,ntriang,iout,nmat,
     $               shlmin,shlmax,ops,s,smhalf,fcoef,
     $               triang,t4,t5,temp,kept,usesym,nirrep,numso,
     $               lambda,occsym,orbsym,salc,values(xs),values(x2),
     $               nbfx,nnpx,values(t6),values(t7),values(t8),
     $               values(t9),values(t10),values(t11))
      end if
c
c..bhl 8/10/89
c
      if(logkey(ops,'linear',.false.,' ')) then
         call enforc(nbf,nnp,s,smhalf,eigval,c,t1,t2)
      end if
c
c----------------------------------------------------------------------c
c                 put out vectors in different order                   c
c                 from putvec. all occupied first with no              c
c                 regard for symmetry. then all the rest               c
c                 symmetry by symmetry                                 c
c----------------------------------------------------------------------c
      if (logkey(ops,'scf=symmetry=kohn-order',.false.,' ')) then
         nsmall=0
         if(logkey(ops,'scf=nsmall',.false.,' ')) then
            nsmall=intkey(ops,'scf=nsmall',1,' ')
         endif      
         call khnvec(c,eigval,orbsym,t1,t2,temp,shlmin,shlmax,nbf,
     $               nshell,nirrep,nsmall,noabort)
      endif
 
c
c..bhl 8/10/89
c
      if(logkey(ops,'linear',.false.,' ')) then
         call enforc(nbf,nnp,s,smhalf,eigval,c,t1,t2)
      end if
c
      if (logkey(ops,'print=m503=vector',.false.,' ')) then
         call iosys('read character "basis function labels" from rwf',
     $               len(bflabl(1))*nbf,0,0,bflabl)
         write (iout,6)
    6    format (5x,'final vector:')
         if (chrkey(ops,'print=scf=vector',' ',' ').eq.'all') then
            nprint=nbf
         else
            nprint=min(nocc+4,nbf)
         end if
         call wvec(c,eigval,nbf,nprint,bflabl,' ')
      end if
c
      if(logkey(ops,'linear',.false.,' ')) then
         call sortbl(c,eigval,t1,t2,t3,t3(1,2),nbf,ops)
         if (logkey(ops,'print=m503=vector',.false.,' ')) then
            write (iout,61)
  61        format (/,5x,'final symmetrized vector :')
            call wvec(c,eigval,nbf,nprint,bflabl,' ')
         end if
      end if
c
c     print the vectors in a manner in which m401 can read them with
c     the rdinp option.
      if(logkey(ops,'scf=punch',.false.,' '))then
         write(iout,62)
  62     format (/,5x,'punching occupied vectors :')
         write(iout,*)'$vectors'
         do 63 i=1,nocc
            write(iout,*)' vector ',i
            write(iout,64)(c(j,i),j=1,nbf)
  63     continue
         write(iout,*)'$end'
  64     format(5(2x,f8.3))
      end if
c
c     ----- store the energy, vector, eigenvalues, lagrangian and density
c           matrices
c
      call iosys('write real f to rwf',nshell,fcoef,0,' ')
      call iosys('write real alpha to rwf',nshell**2,alpha,0,' ')
      call iosys('write real beta to rwf',nshell**2,beta,0,' ')
      call iosys('write real energy to rwf',1,energy,0,' ')
      call iosys('write real "hf energy" to rwf',1,energy,0,' ')
      call iosys('write real "hf 1e energy" to rwf',1,eonel,0,' ')
      call iosys('write real "hf 2e energy" to rwf',1,etwoel,0,' ')
      call putvec(c,eigval,nbf,noabort)
c
c
      return
      end

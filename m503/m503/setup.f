*deck @(#)setup.f	5.1  11/6/94
      subroutine setup(shlnbf,shlmin,shlmax,nshell,nbf,alpha,beta,f,
     $                 calc,ops,xn,xd,nocc,occsym,numso,lambda,
     $                 lirrep,nirrep,usesym,orbsym,salc)
c
c***begin prologue     setup
c***date written       850601  (yymmdd)
c***revision date      920802  (yymmdd)
c
c  02 august 1992    rlm at lanl
c     getting default occsym from the rwf.
c  01 august 1992    rlm at lanl
c     adding symmetry restricted scf capabilities
c  05 oct 1988       bhl at llnl
c     adding general scf code
c  12 may 1987       pws at lanl
c     adding tcscf portion to the code.
c  11 december 1986  pws at lanl
c     writing out the 'shlmin' and 'shlmax' arrays.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)setup.f	5.1   11/6/94
c***purpose
c
c            to initialize the symmetry, min and shlmax arrays and the alpha and
c            beta matrices for an scf calculation.
c
c***description
c
c***references
c***routines called
c***end prologue      setup
c
      implicit integer (a-z)
c
      common /io/ inp,iout
      logical logkey
      logical usesym
c
      character*(*) calc
      character*(*) ops
      character*(*) lirrep(nirrep)
      real*8 xn(*),xd(*)
      real*8 c1,c2
      real*8 alpha(nshell,nshell),beta(nshell,nshell),f(nshell)
      real*8 salc(nbf,nbf)
      integer shlnbf(nshell),shlmin(nshell),shlmax(nshell)
      integer numso(nirrep),lambda(nirrep),occsym(nirrep,nshell)
      integer orbsym(nbf)
c
c     ----- retrieve the number of alpha and beta electrons, and define
c           the number of doubly and singly occupied orbitals
      call iosys('read integer "number of alpha electrons" from rwf',
     $            1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $            1,nbe,0,' ')
c
c     ----- setup the shell data -----
      if (calc.eq.'gvb') then
c
c        ----- tcscf -----
         shlnbf(1)=nbe-1
         shlnbf(2)=1
         shlnbf(3)=1
      else if (calc.eq.'general') then
c
c        ----- general scf -----
         if(nbe.ne.0) then
            shlnbf(1)=nbe
            shlnbf(2)=nae-nbe
         else
            shlnbf(1)=nae
         end if
         if(logkey(ops,'scf=shlnbf',.false.,' ')) then
            call intarr(ops,'scf=shlnbf',shlnbf,nshell-1,' ')
         end if
      else
c
c        ----- closed-shell and high-spin open shell -----
         shlnbf(1)=nbe
         shlnbf(2)=nae-nbe
      end if
c
      nocc=0
      do 1 i=1,nshell-1
         nocc=nocc+shlnbf(i)
    1 continue
      shlnbf(nshell)=nbf-nocc
c
      shlmin(1)=1
      shlmax(nshell)=nbf
      do 2 i=2,nshell
         shlmin(i)=shlmin(i-1)+shlnbf(i-1)
         shlmax(i-1)=shlmin(i)-1
    2 continue
c
c     ----- store shlnbf on the rwf for later use in drt -----
      call iosys('write integer "number in shell" to rwf',nshell,
     $           shlnbf,0,' ')
      call iosys('write integer "first of shell" to rwf',nshell,
     $            shlmin,0,' ')
      call iosys('write integer "last of shell" to rwf',nshell,
     $            shlmax,0,' ')
c
c     ----- setup occupation numbers and coupling coefficients -----
      call rzero(f,nshell)
      call rzero(alpha,nshell**2)
      call rzero(beta,nshell**2)
c
      if (calc.eq.'closed') then
         f(1)=1.0d+00
         alpha(1,1)=2.0d+00
         beta(1,1)=-1.0d+00
      else if (calc.eq.'open') then
         f(1)=1.0d+00
         f(2)=0.5d+00
         alpha(1,1)=2.0d+00
         alpha(2,1)=1.0d+00
         alpha(2,2)=0.5d+00
         beta(1,1)=-1.0d+00
         beta(2,1)=-0.5d+00
         beta(2,2)=-0.5d+00
      else if (calc.eq.'general') then
         if(nshell.ne.2) then
            f(1)=1.0d+00
            f(2)=0.5d+00
            alpha(1,1)=2.0d+00
            alpha(2,1)=1.0d+00
            alpha(2,2)=0.5d+00
            beta(1,1)=-1.0d+00
            beta(2,1)=-0.5d+00
            beta(2,2)=-0.5d+00
         else
            f(1)=.50d0
            alpha(1,1)=.5d0
            beta(1,1)=-.5d0
         end if
c
         ntshl=nshell*(nshell-1)/2
         if(logkey(ops,'scf=f',.false.,' ')) then
            call fparr(ops,'scf=f',f,nshell-1,' ')
         end if
         if(logkey(ops,'scf=fn',.false.,' ')) then
            call fparr(ops,'scf=fn',xn,nshell-1,' ')
            call fparr(ops,'scf=fd',xd,nshell-1,' ')
            do 51 i=1,nshell-1
               f(i)=xn(i)/xd(i)
  51        continue
         end if
         if(logkey(ops,'scf=alpha',.false.,' ')) then
            call fparr(ops,'scf=alpha',xn,ntshl,' ')
            ij=0
            do 52 i=1,nshell-1
               do 53 j=1,i
                  ij=ij+1
                  alpha(i,j)=xn(ij)
  53           continue
  52        continue
         end if
         if(logkey(ops,'scf=an',.false.,' ')) then
            call fparr(ops,'scf=an',xn,ntshl,' ')
            call fparr(ops,'scf=ad',xd,ntshl,' ')
            ij=0
            do 62 i=1,nshell-1
               do 63 j=1,i
                  ij=ij+1
                  alpha(i,j)=xn(ij)/xd(ij)
  63           continue
  62        continue
         end if
         if(logkey(ops,'scf=beta',.false.,' ')) then
            call fparr(ops,'scf=beta',xn,ntshl,' ')
            ij=0
            do 72 i=1,nshell-1
               do 73 j=1,i
                  ij=ij+1
                  beta(i,j)=xn(ij)
  73           continue
  72        continue
         end if
         if(logkey(ops,'scf=bn',.false.,' ')) then
            call fparr(ops,'scf=bn',xn,ntshl,' ')
            call fparr(ops,'scf=bd',xd,ntshl,' ')
            ij=0
            do 82 i=1,nshell-1
               do 83 j=1,i
                  ij=ij+1
                  beta(i,j)=xn(ij)/xd(ij)
  83           continue
  82        continue
         end if
      else if (calc.eq.'gvb') then
c
c        ----- tcscf ... assume some coefficients to stop pseud blowing up
         c1=0.9d+00
         c2=-sqrt(1.0d+00-c1**2)
         f(1)=1.0d+00
         f(2)=c1**2
         f(3)=c2**2
         alpha(1,1)=2.0d+00
         alpha(2,1)=2.0d+00*f(2)
         alpha(2,2)=f(2)
         alpha(3,1)=2.0d+00*f(3)
         alpha(3,2)=0.0d+00
         alpha(3,3)=f(3)
         beta(1,1)=-1.0d+00
         beta(2,1)=-f(2)
         beta(2,2)=0.0d+00
         beta(3,1)=-f(3)
         beta(3,2)=-c1*c2
         beta(3,3)=0.0d+00
      end if
c
      do 20 i=2,nshell
         do 10 j=1,i-1
            alpha(j,i)=alpha(i,j)
            beta(j,i)=beta(i,j)
 10      continue
 20   continue
c
      call iosys('write real f to rwf',nshell,f,0,' ')
      call iosys('write real alpha to rwf',nshell**2,alpha,0,' ')
      call iosys('write real beta to rwf',nshell**2,beta,0,' ')
c
c     ----- print general coupling coefficients -----
      if(calc.eq.'general') then
         write(iout,901)
         write(iout,902)(shlnbf(i),i=1,nshell)
         write(iout,903)
         call matout(f,1,nshell,1,nshell,iout)
         write(iout,904)
         call matout(alpha,nshell,nshell,nshell,nshell,iout)
         write(iout,905)
         call matout(beta,nshell,nshell,nshell,nshell,iout)
  901    format(/,'       --- general scf calculation ---',
     $          /,'       number of orbitals in each shell')
  902    format(5x,10(2x,i5))
  903    format(/,' general scf coupling constants ',/,
     $            '             f           ')
  904    format(/,'           alpha ')
  905    format(/,'           beta ')
c
      end if
c
c     ----- retrieve some symmetry information -----
      if(usesym) then
         call iosys('read integer "number of symmetry orbitals"'
     $              //' from rwf',nirrep,numso,0,' ')
         call iosys('read integer "degeneracies of irreducible'//
     $              ' representations" from rwf',nirrep,lambda,0,' ')
         call iosys('read character "labels of irreducible '//
     $              'representations" from rwf',nirrep*len(lirrep(1)),
     $               0,0,lirrep)
         call iosys('read real "orthogonal salc transformation matrix"'
     $              //' from rwf',-1,salc,0,' ')
         call iosys('read integer "guess vector symmetries" from rwf',
     $               nbf,orbsym,0,' ')
c
c        ----- if requested, read the number of occupied orbitals to be
c              kept in each symmetry, otherwise get them from the
c              guess vector
         call izero(occsym,nirrep*nshell)
         if(logkey(ops,'scf=occsym',.false.,' ')) then
            call intarr(ops,'scf=occsym',occsym,nirrep*nshell,' ')
         else
            do 92 shell=1,nshell
               do 91 i=shlmin(shell),shlmax(shell)
                  occsym(orbsym(i),shell)=occsym(orbsym(i),shell)+1
   91          continue
   92       continue
         endif
         do 93 i=1,nshell
            write(iout,1000) i
            write(iout,1005) (lirrep(j),j=1,nirrep)
            write(iout,1010) (occsym(j,i),j=1,nirrep)
 1000       format(5x,'shell',i2,' occupations by symmetry')
 1005       format(8x,18a4)
 1010       format(8x,18i4)
   93    continue
      endif
c
c
      return
      end

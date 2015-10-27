*deck @(#)setup.f	5.1  11/28/95
      subroutine setup(shlnbf,shlmin,shlmax,nshell,nbf,nae,nbe,
     $                 ops,calc,occsym,numso,lambda,
     $                 lirrep,nirrep,usesym,orbsym,salc,
     $                 alpha,beta,f)
c
c***begin prologue     setup.f
c***date written       930511  (yymmdd) 
c***revision date      11/6/94      
c   12 june, 1993      rlm at lanl
c      adding open shell and general scf coupling coefficients
c
c***keywords           symmetry, shells, dft 
c***author             martin, richard (lanl) 
c***source             @(#)setup.f	5.1   11/28/95
c***purpose            to initialize the symmetry and shell
c                      occupation arrays
c***description
c     
c    
c
c***references
c
c***routines called
c                      intarr,iosys,izero,logkey
c***end prologue       setup.f
      implicit none
c     --- input variables -----
      integer nshell,nbf,nae,nbe,nirrep
      logical usesym
      character*(*) ops
      character*(*) calc
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer shlnbf(nshell),shlmin(nshell),shlmax(nshell)
      integer numso(nirrep),lambda(nirrep),occsym(nirrep,nshell)
      integer orbsym(nbf)
      real*8 salc(nbf,nbf)
      real*8 alpha(nshell,nshell),beta(nshell,nshell),f(nshell)
c     --- output variables ---
      character*(*) lirrep(nirrep)
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,ij,nocc,ntshl,shell
      integer mxshlt
      parameter (mxshlt=3*4/2)
      integer inp,iout
      real*8 one,two,half
      real*8 xn(mxshlt),xd(mxshlt)
      logical logkey
c
c
      include 'msgtypesf.h'
      integer mynodeid,nodeid, mitob, mdtob
      logical ispar
      common /tcgmesa/ ispar
c 
      common /io/ inp,iout
      parameter (one=1.0d0,two=2.0d0,half=0.5d0)
c
c
 1010 format(5x,'--- general scf calculation ---',
     $       /8x,' number of orbitals in each shell')
 1020 format(5x,10(2x,i5))
 1030 format(8x,' general scf coupling constants ',
     $      /9x,'f    ')
 1040 format(9x,'alpha ')
 1050 format(9x,'beta ')
 2000 format(5x,'shell',i2,' occupations by symmetry')
 2005 format(8x,3x,18a5)
 2010 format(8x,18i5)
c
      mynodeid=nodeid()
c     ----- setup the shell data -----
      if (calc.eq.'general') then
c if we get here we're hosed, cos we were supposed to die before this
c$$$c
c$$$c        --- general scf ---
c$$$         if(nbe.ne.0) then
c$$$            shlnbf(1)=nbe
c$$$            shlnbf(2)=nae-nbe
c$$$         else
c$$$            shlnbf(1)=nae
c$$$         end if
c$$$         if(logkey(ops,'scf=shlnbf',.false.,' ')) then
c$$$            call intarr(ops,'scf=shlnbf',shlnbf,nshell-1,' ')
c$$$         end if
      else
c
c        --- closed-shell and high-spin open shell ---
         shlnbf(1)=nbe
         shlnbf(2)=nae-nbe
      end if
c
c     --- determine the number of virtuals ---
      nocc=0
      do 10 i=1,nshell-1
         nocc=nocc+shlnbf(i)
   10 continue
      shlnbf(nshell)=nbf-nocc
c
c     --- setup the shell min/max arrays ---
      shlmin(1)=1
      shlmax(nshell)=nbf
      do 20 i=2,nshell
         shlmin(i)=shlmin(i-1)+shlnbf(i-1)
         shlmax(i-1)=shlmin(i)-1
   20 continue
c
      if (mynodeid .eq. 0) then
c     --- store shlnbf on the rwf for later use in drt ---
         call iosys('write integer "number in shell" to rwf',nshell,
     $        shlnbf,0,' ')
         call iosys('write integer "first of shell" to rwf',nshell,
     $        shlmin,0,' ')
         call iosys('write integer "last of shell" to rwf',nshell,
     $        shlmax,0,' ')
      endif
c
c     --- setup occupation numbers and coupling coefficients ---
      call rzero(f,nshell)
      call rzero(alpha,nshell**2)
      call rzero(beta,nshell**2)
c
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
c$$$         if(nshell.ne.2) then
c$$$            f(1)=1.0d+00
c$$$            f(2)=0.5d+00
c$$$            alpha(1,1)=2.0d+00
c$$$            alpha(2,1)=1.0d+00
c$$$            alpha(2,2)=0.5d+00
c$$$            beta(1,1)=-1.0d+00
c$$$            beta(2,1)=-0.5d+00
c$$$            beta(2,2)=-0.5d+00
c$$$         else
c$$$            f(1)=.50d0
c$$$            alpha(1,1)=.5d0
c$$$            beta(1,1)=-.5d0
c$$$         end if
      end if
c
c        --- perhaps read the occupations and coupling coefficients ---
      ntshl=nshell*(nshell-1)/2
      if (mynodeid .eq. 0) then
         if(logkey(ops,'scf=f',.false.,' ')) then
            call fparr(ops,'scf=f',f,nshell-1,' ')
         end if
         if(logkey(ops,'scf=fn',.false.,' ')) then
            call fparr(ops,'scf=fn',xn,nshell-1,' ')
            call fparr(ops,'scf=fd',xd,nshell-1,' ')
            do 50 i=1,nshell-1
               f(i)=xn(i)/xd(i)
 50         continue
         end if
         if(logkey(ops,'scf=alpha',.false.,' ')) then
            call fparr(ops,'scf=alpha',xn,ntshl,' ')
            ij=0
            do 60 i=1,nshell-1
               do 55 j=1,i
                  ij=ij+1
                  alpha(i,j)=xn(ij)
 55            continue
 60         continue
         end if
         if(logkey(ops,'scf=an',.false.,' ')) then
            call fparr(ops,'scf=an',xn,ntshl,' ')
            call fparr(ops,'scf=ad',xd,ntshl,' ')
            ij=0
            do 70 i=1,nshell-1
               do 65 j=1,i
                  ij=ij+1
                  alpha(i,j)=xn(ij)/xd(ij)
 65            continue
 70         continue
         end if
         if(logkey(ops,'scf=beta',.false.,' ')) then
            call fparr(ops,'scf=beta',xn,ntshl,' ')
            ij=0
            do 80 i=1,nshell-1
               do 75 j=1,i
                  ij=ij+1
                  beta(i,j)=xn(ij)
 75            continue
 80         continue
         end if
         if(logkey(ops,'scf=bn',.false.,' ')) then
            call fparr(ops,'scf=bn',xn,ntshl,' ')
            call fparr(ops,'scf=bd',xd,ntshl,' ')
            ij=0
            do 90 i=1,nshell-1
               do 85 j=1,i
                  ij=ij+1
                  beta(i,j)=xn(ij)/xd(ij)
 85            continue
 90         continue
         end if
c     
c     symmetrize the coupling coefficients
         do 120 i=2,nshell
            do 110 j=1,i-1
               alpha(j,i)=alpha(i,j)
               beta(j,i)=beta(i,j)
 110        continue
 120     continue
c     
         call iosys('write real f to rwf',nshell,f,0,' ')
         call iosys('write real alpha to rwf',nshell**2,alpha,0,' ')
         call iosys('write real beta to rwf',nshell**2,beta,0,' ')
      endif
      if (ispar) then
         call brdcst(200+MSGDBL,f,mdtob(nshell),0)
         call brdcst(201+MSGDBL,alpha,mdtob(nshell*nshell),0)
         call brdcst(202+MSGDBL,beta,mdtob(nshell*nshell),0)
      endif
c
c     --- print general coupling coefficients ---
      if(calc.eq.'general') then
c$$$         write(iout,1010)
c$$$         write(iout,1020)(shlnbf(i),i=1,nshell)
c$$$         write(iout,1030)
c$$$         call matout(f,1,nshell,1,nshell,iout)
c$$$         write(iout,1040)
c$$$         call matout(alpha,nshell,nshell,nshell,nshell,iout)
c$$$         write(iout,1050)
c$$$         call matout(beta,nshell,nshell,nshell,nshell,iout)
c
      end if
c
c
c     --- retrieve some symmetry information ---
      if(usesym) then
         if (mynodeid .eq. 0) then
            call iosys('read integer "number of symmetry orbitals"'
     $           //' from rwf',nirrep,numso,0,' ')
            call iosys('read integer "degeneracies of irreducible'//
     $           ' representations" from rwf',nirrep,lambda,0,' ')
            call iosys('read character "labels of irreducible '//
     $           'representations" from rwf',nirrep*len(lirrep(1)),
     $           0,0,lirrep)
            call iosys(
     $           'read real "orthogonal salc transformation matrix"'
     $           //' from rwf',-1,salc,0,' ')
            call iosys(
     $           'read integer "guess vector symmetries" from rwf',
     $           nbf,orbsym,0,' ')
c     
c       --- if requested, read the number of occupied orbitals to be
c           kept in each symmetry, otherwise get them from the
c           guess vector
            call izero(occsym,nirrep*nshell)
            if(logkey(ops,'scf=occsym',.false.,' ')) then
               call intarr(ops,'scf=occsym',occsym,nirrep*nshell,' ')
            else
               do 150 shell=1,nshell
                  do 160 i=shlmin(shell),shlmax(shell)
                     occsym(orbsym(i),shell)=occsym(orbsym(i),shell)+1
 160              continue
 150           continue
            endif
            do 180 i=1,nshell
               write(iout,2000) i
               write(iout,2005) (lirrep(j),j=1,nirrep)
               write(iout,2010) (occsym(j,i),j=1,nirrep)
 180        continue
         endif
      endif
c
c
      return
      end

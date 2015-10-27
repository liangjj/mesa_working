*deck %W%  %G%
      subroutine pm1010(z,a)
c
c***begin prologue     m1010
c***date written       860819   (yymmdd)
c***revision date      890119   (yymmdd)
c
c  19 january  1989    bhl at llnl
c     test inserted in subroutine lagrng for numi=0
c
c  29 march    1988    bhl at llnl
c     checking if files exist ( mo derivative lagrangian ) before
c     attempting to create then so mesa can perform ci optimizations
c
c  15 november 1987    pws at lanl
c      adding formation of derivative core lagrangian.
c
c   8 december 1986    pws at lanl
c      modifying m601, the scf lagrangian code, to use derivative ao integrals
c      and form the 'derivative' lagrangian matrices needed for the cphf
c      equations.
c
c***keywords           m1010, link 1010, scf lagrangian, hf lagrangian
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c***purpose            to form the lagrangian for general-fock scf
c***description
c     m1010 recognizes the options subtrings:
c     timing                 collect and print timing statistics
c
c***references         y. osamura, y. yamaguchi, p. saxe, m. a. vincent,
c                      j. f. gaw and h. f. schaefer iii, "unified theoretical
c                      treatment of analytic first and second energy
c                      derivatives in open-shell hartree-fock theory",
c                      chemical physics 72 (1982) 131-139.
c
c***routines called
c***end prologue       m1010
c
      implicit integer (a-z)
c
      character*4096 ops
      character*8 prtflg
      character*128 namdnt
      character*2 ians
      real*8 z(*)
      integer a(*)
      logical prnt
      logical logkey
c
      common /io/     inp,iout
c
c
      data maxcor /1/
c
    2 format(1x,'m1010:')
    3 format(5x,'memory use',18x,i9)
    4 format(5x,'all integrals held in core.')
    5 format(5x,'# integral triangles in core',i4)
c
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- print turned off externally? -----
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      prnt=prtflg.ne.'minimum'
c
c     ----- start timing routines if enabled -----
c
c
c     ----- open the integral file -----
c
      call iosys('read character "raw derivative integral filename"'
     $         //' from rwf',0,0,0,namdnt)
      call iosys('open rdints as old',0,0,0,namdnt)
c
c     ----- get the dimensions, etc -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $            1,num,0,' ')
      nnp=(num+1)*num/2
      call iosys('read integer "spin multiplicity" from rwf',
     $           1,multip,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
c
c     ----- set up some parameters depending on multip -----
c
      call iosys('read integer "number of shells" from rwf',
     $           1,nshell,0,' ')
c
      if (multip.eq.1) then
         scfnum=0
         nfock=1
         ncoul=1
         nexch=1
         ndmat=1
      else
         scfnum=1
         nfock=2
         ncoul=2
         nexch=2
         ndmat=2
      end if
c
c      call iosys('write integer ndmat to rwf',1,ndmat,0,' ')
c      call iosys('write integer no_types_orbitals to rwf',1,nshell,0,' ')
c
c     ----- allocate core -----
c
      numshl=1
      minshl=numshl+nshell
      maxshl=minshl+nshell
      grad=iadtwp(maxshl+nshell)
      zan=grad+3*natoms
      c=zan+natoms
      f=c+3*natoms
      alpha=f+nshell
      beta=alpha+nshell**2
      h=beta+nshell**2
      c=h+nnp
      d=c+num**2
      j=d+nnp*ndmat
      k=j+nnp*ncoul
      t1=k+nnp*nexch
      t2=t1+num**2
      lag=t2+num**2
      clag=lag+num**2
      values=clag+num**2
c
c     ----- find out how much core is available -----
c
      call getscm(0,z,canget,'m1010: how much core',0)
c
      top=min(canget,wpadti(values+nnp**2)+100)
c
c     ----- get the core needed -----
c
      call getscm(top,z,maxcor,'m1010',1)
c
      left=iadtwp(maxcor)-values
      ntriang=min(left/nnp,nnp)
      lenbuf=ntriang*nnp
c
      if (prnt) then
         write(iout,2)
         write(iout,3) maxcor
         if(ntriang.eq.nnp) then
            write(iout,4)
         else
            write(iout,5) ntriang
         endif
      endif
      if (ntriang.lt.1) call lnkerr('not enough core in m1010')
c
c     ----- read in atomic charges and coordinates for
c           nuclear repulsion term
c
      call iosys('read real "nuclear charges" from rwf',-1,z(zan),0,
     $           ' ')
      call iosys('read real coordinates from rwf',-1,z(c),0,' ')
c
c     ----- nuclear repulsion contribution to the gradient -----
c
      call nucrep(z(zan),z(c),natoms,z(grad))
c
c     ----- initialize the arrays -----
c
      call iosys('read real f from rwf',nshell,z(f),0,' ')
      call iosys('read real alpha from rwf',nshell**2,z(alpha),0,' ')
      call iosys('read real beta from rwf',nshell**2,z(beta),0,' ')
      call iosys('read integer "number in shell" from rwf',nshell,
     $           a(numshl),0,' ')
      call iosys('read integer "first of shell" from rwf',nshell,
     $            a(minshl),0,' ')
      call iosys('read integer "last of shell" from rwf',nshell,
     $            a(maxshl),0,' ')
c
c     ----- create space to write the lagrangians -----
c
c..bhl
      call iosys('does "mo derivative hf lagrangian" exist on rwf',
     #           0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create real "mo derivative hf lagrangian" on rwf',
     $              num**2*3*natoms,0,0,' ')
      else
         call iosys('rewind "mo derivative hf lagrangian" on rwf',
     #               0,0,0,' ')
      end if
c
      call iosys('does "ao derivative hf lagrangian" exist on rwf',
     #            0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create real "ao derivative hf lagrangian" on rwf',
     #               num**2*3*natoms,0,0,' ')
      else
         call iosys('rewind "ao derivative hf lagrangian" on rwf',
     #               0,0,0,' ')
      end if
c
      call iosys('does "mo derivative core lagrangian" exist on rwf',
     #            0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create real "mo derivative core lagrangian" '
     $              //'on rwf',num**2*3*natoms,0,0,' ')
      else
         call iosys('rewind "mo derivative core lagrangian" on rwf',
     #               0,0,0,' ')
      end if
c..bhl
c
c     ----- form the lagrangian -----
c
      do 100 der=1,3*natoms
        call lagrng(a(numshl),a(minshl),a(maxshl),z(f),z(alpha),
     #        z(beta),z(h),z(c),z(d),z(j),z(k),z(values),num,nnp,
     #        nshell,ndmat,ncoul,nexch,scfnum,ntriang,z(t1),z(t2),
     #        z(lag),ops,der,z(grad),3*natoms,z(clag))
  100 continue
c
c     ----- and the overlap contribution to the scf gradients -----
c
      call sder(z(h),z(d),num,nnp,z(t1),z(t2),z(grad),3*natoms)
c
      if (logkey(ops,'print=gradients=scf-from-lagrangians',
     $     .false.,' ')) then
         write (iout,110)
 110     format (//,t10,'scf derivatives from lagrangians')
         call matout(z(grad),3,natoms,3,natoms,iout)
      end if
c
c
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      stop
      end

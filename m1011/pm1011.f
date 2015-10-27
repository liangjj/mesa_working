*deck @(#)pm1011.f	5.1  11/6/94
      subroutine pm1011(z,a)
c***begin prologue     m1011
c***date written       860819   (yymmdd)
c***revision date      861208   (yymmdd)
c
c   8 december 1986   pws at lanl
c     changing 'namint' and iosys open to character.
c
c***keywords           m1011, link 1011, scf lagrangian, hf lagrangian
c***author             saxe, paul (lanl)
c***source             @(#)pm1011.f	5.1   11/6/94
c***purpose            to form the generalised lagrangians for general-fock scf
c***description
c     m1011 recognizes the options subtrings:
c     timing                 collect and print timing statistics
c
c***references         y. osamura, y. yamaguchi, p. saxe, m. a. vincent,
c                      j. f. gaw and h. f. schaefer iii, "unified theoretical
c                      treatment of analytic first and second energy
c                      derivatives in open-shell hartree-fock theory",
c                      chemical physics 72 (1982) 131-139.
c
c***routines called
c***end prologue       m1011
c
      implicit integer (a-z)
c
      character*4096 ops
      character*8 prtflg
      character*128 namint
      character*8 calc
      real*8 z(*)
      integer a(*)
      logical prnt
c
      common /io/     inp,iout
c
c
    2 format(1x,'m1011:')
    3 format(5x,'memory use',18x,i9)
    4 format(5x,'all integrals held in core.')
    5 format(5x,'# integral triangles in core',i4)
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
c     ----- open the integral file -----
c
      call iosys('read character "integral filename" from rwf',
     $     0,0,0,namint)
      call iosys('open ints as old',0,0,0,namint)
c
c     ----- get the dimensions, etc -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
      call iosys('read integer "spin multiplicity" from rwf',
     $     1,multip,0,' ')
c
c     ----- set up some parameters depending on multip -----
c
      call iosys('read integer "number of shells" from rwf',
     $     1,nshell,0,' ')
c
      if (multip.eq.1) then
         calc='closed'
         ncoul=1
         nexch=1
         ndmat=1
      else
         calc='open'
         ncoul=2
         nexch=2
         ndmat=2
      end if
c
c     ----- allocate core -----
c
      numshl=1
      minshl=numshl+nshell
      maxshl=minshl+nshell
      f=iadtwp(maxshl+nshell)
      alpha=f+nshell
      beta=alpha+nshell**2
      h=beta+nshell**2
      c=h+nnp
      d=c+nbf**2
      j=d+nnp*ndmat
      k=j+nnp*ncoul
      t1=k+nnp*nexch
      t2=t1+nbf**2
      lag=t2+nbf**2
      values=lag+nbf**2
      top=wpadti((values+nnp*nnp)+100)
c
c     ----- find out how much core is available -----
c
      call getscm(0,z,maxcor,'m1011: how much core',0)
c
c NONONONONO.  This ASSumes that we can get all the triangles into core.
c That will be an oddity.  Postpone the puke until AFTER we determine how
c many triangles we can fit.
c$$$      if(top.gt.maxcor) then
c$$$         write(iout,*)"m1011: available=",maxcor
c$$$         write(iout,*)"Need=",top
c$$$         call lnkerr('need more core in m1011')
c$$$      endif
c
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
      if (ntriang.lt.1) call lnkerr('not enough core in m1011')
c
c     ----- read in t and v one-electron integrals -----
c
      call iosys('read real "kinetic integrals" from rwf',
     $           -1,z(values),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $           -1,z(h),0,' ')
      call vadd(z(h),z(h),z(values),nnp)
c
c     ----- initialize the arrays -----
c
      call iosys('read real f from rwf',nshell,z(f),0,' ')
      call iosys('read real alpha from rwf',nshell**2,z(alpha),0,' ')
      call iosys('read real beta from rwf',nshell**2,z(beta),0,' ')
      call iosys('read integer "number in shell" from rwf',nshell,
     $            a(numshl),0,' ')
      call iosys('read integer "first of shell" from rwf',nshell,
     $            a(minshl),0,' ')
      call iosys('read integer "last of shell" from rwf',nshell,
     $            a(maxshl),0,' ')
c
c     ----- form the lagrangian -----
c
      scfnum=1
      call lagrng(a(numshl),a(minshl),a(maxshl),z(f),z(alpha),z(beta),
     #            z(h),z(c),z(d),z(j),z(k),z(values),nbf,nnp,nshell,
     #            ndmat,ncoul,nexch,ntriang,z(t1),z(t2),
     #            z(lag),ops)

c
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      return
      end

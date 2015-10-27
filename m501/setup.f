*deck @(#)setup.f	5.1  11/6/94
      subroutine setup(shlnbf,shlmin,shlmax,nshell,nbf,alpha,beta,f,
     $                 scfnum,ops,xn,xd,nocc)
c
c***begin prologue     setup
c***date written       850601  (yymmdd)
c***revision date      881005  (yymmdd)
c
c  05 oct 1988       bhl at llnl
c   adding general scf code
c
c  12 may 1987       pws at lanl
c   adding tcscf portion to the code.
c
c  11 december 1986  pws at lanl
c     writing out the 'shlmin' and 'shlmax' arrays.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)setup.f	5.1   11/6/94
c***purpose
c
c            to initialize the min and shlmax arrays and the alpha and
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
c
      character*(*) ops
      real*8 xn(*),xd(*)
      real*8 c1,c2
      real*8 alpha(nshell,nshell),beta(nshell,nshell),f(nshell)
      integer shlnbf(nshell),shlmin(nshell),shlmax(nshell)
c
c     ----- start timing -----
c
c
c     ----- retrieve the number of alpha and beta electrons, and define
c           the number of doubly and singly occupied orbitals
c
      call iosys('read integer "number of alpha electrons" from rwf',
     $            1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $            1,nbe,0,' ')
c
c
      if (scfnum.eq.2.or.scfnum.eq.3) then
c
c        ----- tcscf -----
c
         shlnbf(1)=nbe-1
         shlnbf(2)=1
         shlnbf(3)=1
      else if (scfnum.eq.6.) then
c
c        ----- general scf -----
c
c..bhl 7/11/89
         if(nbe.ne.0) then
            shlnbf(1)=nbe
            shlnbf(2)=nae-nbe
         else
            shlnbf(1)=nae
         end if
c..bhl 7/11/89
         if(logkey(ops,'scf=shlnbf',.false.,' ')) then
            call intarr(ops,'scf=shlnbf',shlnbf,nshell-1,' ')
         end if
c
      else
c
c        ----- closed-shell and high-spin open shell -----
c
         shlnbf(1)=nbe
         shlnbf(2)=nae-nbe
      end if
c
c
      nocc=0
      do 1 i=1,nshell-1
         nocc=nocc+shlnbf(i)
    1 continue
      shlnbf(nshell)=nbf-nocc
c
c
      shlmin(1)=1
      shlmax(nshell)=nbf
      do 2 i=2,nshell
         shlmin(i)=shlmin(i-1)+shlnbf(i-1)
         shlmax(i-1)=shlmin(i)-1
    2 continue
c
c     ----- store shlnbf on the rwf for later use in drt -----
c
      call iosys('write integer "number in shell" to rwf',nshell,
     $           shlnbf,0,' ')
      call iosys('write integer "first of shell" to rwf',nshell,
     $            shlmin,0,' ')
      call iosys('write integer "last of shell" to rwf',nshell,
     $            shlmax,0,' ')
c
      call rzero(f,nshell)
      call rzero(alpha,nshell**2)
      call rzero(beta,nshell**2)
c
      if (scfnum.eq.0.or.scfnum.eq.5) then
         f(1)=1.0d+00
         alpha(1,1)=2.0d+00
         beta(1,1)=-1.0d+00
      else if (scfnum.eq.1.or.scfnum.eq.4) then
         f(1)=1.0d+00
         f(2)=0.5d+00
         alpha(1,1)=2.0d+00
         alpha(2,1)=1.0d+00
         alpha(2,2)=0.5d+00
         beta(1,1)=-1.0d+00
         beta(2,1)=-0.5d+00
         beta(2,2)=-0.5d+00
      else if (scfnum.eq.6) then
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
c
      else if (scfnum.eq.2.or.scfnum.eq.3) then
c
c        ----- tcscf ... assume some coefficients to stop pseud blowing up
c
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
      if(scfnum.eq.6) then
         write(iout,901)
         write(iout,902)(shlnbf(i),i=1,nshell)
         write(iout,91)
         call matout(f,1,nshell,1,nshell,iout)
         write(iout,92)
         call matout(alpha,nshell,nshell,nshell,nshell,iout)
         write(iout,93)
         call matout(beta,nshell,nshell,nshell,nshell,iout)
  91     format(/,'  general scf coupling constants ',/,
     #            '              f           ')
  92     format(/,'            alpha ')
  93     format(/,'            beta ')
 901     format(//,'        --- general scf calculation ---',/,
     $             '        number of orbitals in each shell')
 902     format(5x,10(2x,i5))
c
      end if
c
c
c     ----- stop timing -----
c
c
c
      return
      end

*deck @(#)hfden.f	1.1  11/30/90
      subroutine hfden(c,d,num,nnp,ndoc,nalp,opnscf,ops)
c
c***begin prologue     hfden
c***date written       871108   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           hf density
c***author             saxe, paul (lanl)
c***source             @(#)hfden.f	1.1   11/30/90
c
c***purpose            to construct high-spin hf density matrices.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       hfden
c
      implicit integer (a-z)
c
      character*(*) ops
      logical logkey
      real*8 c(num,num)
      real*8 d(nnp,2)
      integer shlnbf(3)
      integer shlmin(3)
      integer shlmax(3)
      real*8 f(3)
      real*8 alpha(9)
      real*8 beta(9)
c
      common /io/ inp,iout
c
      if (opnscf.eq.0) then
         call iosys('read real "mcscf vector" from rwf',num**2,c,0,' ')
      else
         call iosys('read real "scf vector" from rwf',num**2,c,0,' ')
      end if
c
      call rzero(d,nnp*2)
c
c     ----- closed shell density matrix -----
c
      ij=0
      do 3 i=1,num
         do 2 j=1,i
            ij=ij+1
            do 1 k=1,ndoc
               d(ij,1)=d(ij,1)+c(i,k)*c(j,k)
 1          continue
 2       continue
 3    continue
c
c     ----- write out either the mcscf core density, or hf junk -----
c
      if (opnscf.eq.0) then
c.bhl
c        density matrix with a consistent scale is written by
c        integral transformation program
c.bhl
c         call iosys('write real "mcscf ao core density" to rwf',
c     $        nnp,d,0,' ')
c
c         if (logkey(ops,'print=mcscf=core-density',.false.,' ')) then
c            write (iout,100)
c 100        format ('1',//,t20,'ao mcscf core density matrix')
c            call print(d,nnp,num,iout)
c..bhl         end if
      else
c
c        ----- open shell density matrix -----
c
         ij=0
         do 6 i=1,num
            do 5 j=1,i
               ij=ij+1
               do 4 k=ndoc+1,ndoc+nalp
                  d(ij,2)=d(ij,2)+c(i,k)*c(j,k)
 4             continue
 5          continue
 6       continue
c
         if (logkey(ops,'print=scf=density',.false.,' ')) then
            write (iout,110)
 110        format ('1',//,t20,'ao scf closed-shell density matrix')
            call print(d,nnp,num,iout)
         end if
c
         if (nalp.le.0) then
            call iosys('write real "hf density matrix" to rwf',
     $           nnp,d,0,' ')
            nhfden=1
         else
            call iosys('write real "hf density matrix" to rwf',
     $           nnp*2,d,0,' ')
            nhfden=2
c
            if (logkey(ops,'print=scf=density',.false.,' ')) then
               write (iout,120)
 120           format ('1',//,t20,'ao scf open-shell density matrix')
               call print(d(1,2),nnp,num,iout)
            end if
c
         end if
c
         call iosys('write integer "number of hf density matrices" '//
     $        'to rwf',1,nhfden,0,' ')
c
c        ----- fix up numbers of shells, etc. -----
c
         call iosys('read integer "number of alpha electrons" from rwf',
     $        1,nae,0,' ')
         call iosys('read integer "number of beta electrons" from rwf',
     $        1,nbe,0,' ')
c
         if (nalp.le.0) then
            nshell=2
         else
            nshell=3
         end if
c
         shlnbf(1)=nbe
         shlnbf(2)=nae-nbe
c
         junk=0
         do 7 i=1,nshell-1
            junk=junk+shlnbf(i)
 7       continue
         shlnbf(nshell)=num-junk
c
         shlmin(1)=1
         shlmax(nshell)=num
         do 8 i=2,nshell
            shlmin(i)=shlmin(i-1)+shlnbf(i-1)
            shlmax(i-1)=shlmin(i)-1
 8       continue
c
c        ----- store shlnbf on the rwf for later use in drt -----
c
         call iosys('write integer "number of shells" to rwf',
     $        1,nshell,0,' ')
         call iosys('write integer "number in shell" to rwf',nshell,
     $        shlnbf,0,' ')
         call iosys('write integer "first of shell" to rwf',nshell,
     $        shlmin,0,' ')
         call iosys('write integer "last of shell" to rwf',nshell,
     $        shlmax,0,' ')
c
         call rzero(f,nshell)
         call rzero(alpha,nshell**2)
         call rzero(beta,nshell**2)
c
         if (nshell.eq.2) then
c
c           ----- closed shell -----
c
            f(1)=1.0d+00
            alpha(1)=2.0d+00
            beta(1)=-1.0d+00
         else
c
c           ----- open-shell -----
c
            f(1)=1.0d+00
            f(2)=0.5d+00
            alpha(1)=2.0d+00
            alpha(2)=1.0d+00
            alpha(4)=1.0d+00
            alpha(5)=0.5d+00
            beta(1)=-1.0d+00
            beta(2)=-0.5d+00
            beta(4)=-0.5d+00
            beta(5)=-0.5d+00
         end if
c
         call iosys('write real f to rwf',nshell,f,0,' ')
         call iosys('write real alpha to rwf',nshell**2,alpha,0,' ')
         call iosys('write real beta to rwf',nshell**2,beta,0,' ')
c
      end if
c
c
      return
      end

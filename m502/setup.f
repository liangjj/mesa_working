*deck @(#)setup.f	5.1  11/6/94
      subroutine setup (occmo,virtmo,num,ndmat,prnt,ops)
      implicit integer(a-z)
      character*(*) ops
      logical logkey,prnt
      real*8 f, aij, bij
      real*8 fn,fd,an,ad,bn,bd
      common /shell/ noshel(10), orbshel(50,10), nvshel(10),
     1 vshel(200,10), f(10), aij(10,10), bij(10,10),
     2 fn(10),fd(10),an(10),ad(10),bn(10),bd(10)
c
      common /io/ inp, iout
c
c*** purpose to initialize data for ivo calculation
c
c     the program assumes that the occupied orbitals are solutions
c     of one or more one electron eigenvalue equations. these are
c     then used to construct an effective hamiltonian in
c     the unoccupied space.
c
c     read in hamiltonian information for scf shells.
c     the last shell read in is the ivo shell.
c
      call iosys('read integer "number of alpha electrons" from rwf',
     #           1,na,0,' ')
c
      occmo=na
      virtmo=num-na
      nclos=na-1
c
      if(logkey(ops,'ivo=singlet',.false.,' ')) then
         spin=1
         if(prnt) write(iout,80)
      else if(logkey(ops,'ivo=neutral',.false.,' ')) then
         spin=0
         if(prnt) write(iout,90)
         nclos=na
      else
         spin=3
         if(prnt) write(iout,100)
      end if
c
c
      call rzero (f,10)
      call rzero (aij,100)
      call rzero (bij,100)
      if(logkey(ops,'ivo=nshell',.false.,' ')) then
c        ..general input
c         nshell=intkey(ops,'ivo=nshell',3,' ')
c       
         ndmat=nshell-1
         if(logkey(ops,'ivo=fn',.false.,' ')) then
            call fparr(ops,'ivo=fn',fn,ndmat,' ')
            call fparr(ops,'ivo=an',an,ndmat,' ')
            call fparr(ops,'ivo=bn',bn,ndmat,' ')
            call fparr(ops,'ivo=fd',fd,ndmat,' ')
            call fparr(ops,'ivo=ad',ad,ndmat,' ')
            call fparr(ops,'ivo=bd',bd,ndmat,' ')
            do 1 i=1,ndmat
               f(i)=fn(i)/fd(i)
               aij(i,1)=an(i)/ad(i)
               bij(i,1)=bn(i)/bd(i)
  1         continue
         else
            call fparr(ops,'ivo=f',f,ndmat,' ')
            call fparr(ops,'ivo=alpha',aij,ndmat,' ')
            call fparr(ops,'ivo=beta',bij,ndmat,' ')
         end if
         call intarr(ops,'ivo=noshell',noshel,nshell,' ')
         occmo=0
         do 2 i=1,ndmat
            occmo=occmo+noshel(i)
  2      continue
         virtmo=noshel(nshell)
c        ..end general input
      else
c        ..default input
         if(nclos.eq.0) then
c           ...special section for 2-e systems
            noshel(1)=1
            noshel(2)=num-1
            ndmat=1
            nshell=2
            f(1)=1.0+00
            aij(1,1)=1.d00
c
            if(spin.eq.3) then
               bij(1,1)=-1.d00
            else
               bij(1,1)=+1.d00
            end if
c
c           ...end special section for 2-e systems
         else if(spin.eq.0) then
c           ...standard input
            nshell=2
            ndmat=1
            noshel(1)=na
            noshel(2)=virtmo
            f(1)=1.0d+00
            aij(1,1)=2.d+00
            bij(1,1)=-1.d+00
         else
            noshel(1)=na-1
            noshel(2)=1
            noshel(3)=virtmo
            ndmat=2
            nshell=3
            f(1)=1.0d+00
            aij(1,1)=2.d+00
            bij(1,1)=-1.d+00
            aij(2,1)=1.d00
c
            if(spin.eq.3) then
               bij(2,1)=-1.d00
            else
               bij(2,1)=+1.d00
            end if
c
c           ...standard coupling
         end if
c     ..end default input
      end if
c
      if(prnt) then
         write (iout,55) f(1)
         write (iout,60)
c
         do 30 i=1,nshell
            write (iout,70) i,noshel(i),aij(i,1),bij(i,1)
 30      continue
      endif
c
c
      return
   55 format (5x,'f vector-coupling-constant  ',f10.5)
   60 format (5x,'shell    no. orbitals ',
     #           '    alpha        beta ')
   70 format (5x,2(3x,i4,3x),2(3x,f10.5))
   80 format (5x,'singlet ivo calculation')
   90 format (5x,'neutral ivo calculation')
  100 format (5x,'triplet ivo calculation')
      end

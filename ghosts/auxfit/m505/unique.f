*deck %W%  %G%
      subroutine unique(c,eigval,f,jmat,kmat,values,d,u,t1,t2,h,
     $                  nbf,nnp,nshell,ntriang,iout,nmat,
     $                  shlmin,shlmax,ops,s,smhalf,fcoef,
     $                  triang,t4,t5,temp,kept,usesym,nirrep,numso,
     $                  lambda,occsym,orbsym,salc,xs,x2,nbfx,nnpx,
     $                  t6,t7,t8,t9,t10,t11)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),eigval(nbf),f(nnp),jmat(nnp,nmat),kmat(nnp,nmat)
      real*8 values(*),d(nnp,nmat),u(nbf,nbf),t1(nbf,nbf)
      real*8 t2(nbf,nbf,nmat),h(nnp),s(nnp),smhalf(nnp),fcoef(nshell)
      real*8 triang(nnp),t4(nbf,nbf),t5(nbf,nbf),salc(nbf,nbf)
      real*8 xs(nnp,nbfx),x2(nbfx,nbfx)
      real*8 t6(nbf*nbf,nbfx),t7(nbfx,nbfx),t8(nbf*nbf,nbfx)
      real*8 t9(nbf*nbf,nbfx),t10(nbf*nbf,nmat),t11(nbf*nbf,nmat)
      integer shlmin(nshell),shlmax(nshell)
      integer temp(nbf),orbsym(nbf),kept(nshell)
      integer numso(nirrep),lambda(nirrep),occsym(nirrep,nshell)
      character*(*) ops
      logical logkey,usesym
c
c     02 august, 1992    rlm at lanl
c        replaced code to diagonalize matrix symmetry with a call to
c        symrsp.
c
c     ----- form the density matrix or matrices -----
c
      do 124 i=1,nshell-1
         call gdmat(d(1,i),c,nbf,nnp,shlmin(i),shlmax(i))
  124 continue
c
c     ----- and the coulomb and exchange matrices -----
c
      call formjk(xs,x2,nbfx,nnpx,d,nnp,nbf,jmat,kmat,nmat,nmat,
     $            nmat,t2,t6,t7,t8,t9,t10,t11)
c
c     ----- transform the j and k matrices to the current mo basis ----
c
      do 1 i=1,nmat
         call trtosq(t1,jmat(1,i),nbf,nnp)
         call ebc(t2,t1,c,nbf,nbf,nbf)
         call ebtc(t1,c,t2,nbf,nbf,nbf)
         call sqtotr(jmat(1,i),t1,nbf,nnp)
c
         call trtosq(t1,kmat(1,i),nbf,nnp)
         call ebc(t2,t1,c,nbf,nbf,nbf)
         call ebtc(t1,c,t2,nbf,nbf,nbf)
         call sqtotr(kmat(1,i),t1,nbf,nnp)
    1 continue
c
c     ----- and construct the averaged fock matrix
c                f=h+f(i)[2j(i)-k(i)]
c
      call trtosq(t1,h,nbf,nnp)
      call ebc(t2,t1,c,nbf,nbf,nbf)
      call ebtc(t1,c,t2,nbf,nbf,nbf)
      call sqtotr(f,t1,nbf,nnp)
c
      do 3 shell=1,nshell-1
         do 2 i=1,nnp
            f(i)=f(i)+fcoef(shell)*(2*jmat(i,shell)-kmat(i,shell))
    2    continue
    3 continue
c
c     ----- zero out the sub-blocks mixing orbital shells -----
c
      do 7 ishell=2,nshell
         do 6 jshell=1,ishell-1
            do 5 i=shlmin(ishell),shlmax(ishell)
               ia=i*(i-1)/2
               do 4 j=shlmin(jshell),shlmax(jshell)
                  f(ia+j)=0.0d+00
    4          continue
    5       continue
    6    continue
    7 continue
c
c
      call rzero(u,nbf*nbf)       
      if (usesym) then
c
c        ----- transform f back into a.o. basis to satisfy the
c        requirements of symrsp for the transformation to
c        the salc basis
         call trtosq(t4,s,nbf,nnp)
         call trtosq(t1,f,nbf,nnp)
         call ebct(t2,t1,c,nbf,nbf,nbf)
         call ebc(t1,c,t2,nbf,nbf,nbf)
         call ebc(t2,t1,t4,nbf,nbf,nbf)
         call ebc(t1,t4,t2,nbf,nbf,nbf)
         call sqtotr(f,t1,nbf,nnp)
         call symrsp(f,u,t1,t2,t4,t5,eigval,temp,kept,nbf,nnp,
     $               nshell,shlmin,shlmax,nirrep,numso,lambda,
     $               orbsym,occsym,salc,c)
      else
             call degrsp(nbf,nnp,f,eigval,1,u,t1,t2)
c
c        ----- transform back to the ao basis -----
         call ebc(t1,c,u,nbf,nbf,nbf)
         call vmove(c,t1,nbf**2)
      endif
c
      if (logkey(ops,'print=scf=core-lagrangian',.false.,' ')) then
         write (iout,1000)
 1000    format (5x,'mo core lagrangian')
         call matout(eigval,nbf,1,nbf,1,iout)
      end if
c
      call iosys('write real "scf mo core lagrangian" to rwf',
     $           nbf,eigval,0,' ')
c
c
c
      return
      end

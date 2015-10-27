*deck %W%  %G%
      subroutine unisym(c,eigval,f,jmat,kmat,values,d,u,t1,t2,h,
     $                  nbf,nnp,nshell,scfnum,ntriang,iout,nmat,
     $                  shlmin,shlmax,ops,s,smhalf,fcoef,
     $                  csym,symorb,symtyp,nsym,nobs)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),eigval(nbf),f(nnp),jmat(nnp,nmat),kmat(nnp,nmat)
      real*8 values(nnp,ntriang),d(nnp,nmat),u(nbf,nbf),t1(nbf,nbf)
      real*8 t2(nbf,nbf,nmat),h(nnp),s(nnp),smhalf(nnp),fcoef(nshell)
      real*8 csym(nbf,nbf),big
      integer symorb(nbf,8),symtyp(nbf)
      integer shlmin(nshell),shlmax(nshell),nobs(nsym)
      integer sympt(8)
      character*(*) ops
      logical logkey
c
c
c     ----- form the density matrix or matrices -----
c
c
      do 124 i=1,nshell-1
         call gdmat(d(1,i),c,nbf,nnp,shlmin(i),shlmax(i))
  124 continue
c
c     ----- and the coulomb and exchange matrices -----
c
      call jandks(values,d,nnp,nbf,jmat,kmat,nmat,nmat,ntriang,
     $            t1,nmat,t2)
c      call fock(values,d,f,nnp,nbf,jmat,kmat,nmat,nmat,ntriang,-nmat,
c     #          1,t1,h,nmat,t2,printj,printk,c1,c2,haa,hbb,hab)
c
c     ----- transform the J and K matrices to the current MO basis ----
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
c                F=h+f(i)[2J(i)-K(i)]
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
c ---- determine the symmetry of the MOs
c
      call trtosq(t1,s,nbf,nnp)
      call ebc(t2,t1,c,nbf,nbf,nbf)
      call ebtc(t1,csym,t2,nbf,nbf,nbf)
c
c
      jof=0
      do 9 i=1,nsym
      sympt(i)=0
       nob=nobs(i)
       do 8 j=jof+1,jof+nob
        symtyp(j)=i
  8    continue
      jof=jof+nob
  9   continue
c
c      write(iout,*)' SYMTYP '
c      write(iout,902)(symtyp(i),i=1,nbf)
c
      call izero(symorb,nbf*8)
c
c      write(iout,*)' Orbital Symmetry NOBS jmax jmin jtest '
c
      do 10 i=1,nbf
c
       jtest=isamax(nbf,t1(1,i),1)
c
       isym=symtyp(jtest)
       js=sympt(isym)+1
       symorb(js,isym)=i
       sympt(isym)=js
c
c       write(iout,902)i,isym,js,jtest
c
 10   continue
c
c      write(iout,*)' SYMPT  '
c      write(iout,902)(sympt(i),i=1,nsym)
c
      do 905 i=1,nsym
       if(sympt(i).ne.nobs(i)) then
        write(iout,*)' SYM sympt nob ',i,sympt(i),nobs(i)
        call lnkerr(' M501: bug unisym ')
       end if
 905  continue
c
c      write(iout,*)' SYMORB '
c      do 901 i=1,nsym
c       nob=nobs(i)
c       write(iout,*) ' SYM NOB ',i,nob
c       write(iout,902)(symorb(j,i),j=1,nob)
c 901  continue
c
 902  format(2x,10i5)
c
c
      call trtosq(u,f,nbf,nnp)
c
      ij=0
      do 11 isym=1,nsym
       nob=nobs(isym)
       do 12 i=1,nob
        is=symorb(i,isym)
        do 13 j=1,i
         js=symorb(j,isym)
         ij=ij+1
         f(ij)=u(is,js)
  13    continue
  12   continue
  11  continue
c
c     ----- and diagonalize this matrix, transforming the scf vectors -
c
      ix=1
      ixx=1
      jx=1
      do 14 isym=1,nsym
       nob=nobs(isym)
       nnps=nob*(nob+1)/2
       call degrsp(nob,nnps,f(ix),eigval(ixx),1,u(jx,1),t1,t2)
       jx=jx+nob*nob
       ix=ix+nnps
       ixx=ixx+nob
  14  continue
c
      ix=1
      jx=1
      do 15 isym=1,nsym
       nob=nobs(isym)
       do 16 i=1,nob
        is=symorb(i,isym)
        call scopy(nbf,c(1,is),1,t2(1,i,1),1)
  16   continue
       call ebc(t1(1,ix),t2,u(jx,1),nbf,nob,nob)
       jx=jx+nob*nob
       ix=ix+nob
  15  continue
c
c
      big=1.d+28
      do 17 i=1,nbf
       jtest=ismin(nbf,eigval,1)
       t2(i,1,1)=eigval(jtest)
       eigval(jtest)=big
       call scopy(nbf,t1(1,jtest),1,c(1,i),1)
  17  continue
c
      call scopy(nbf,t2,1,eigval,1)
c
      if (logkey(ops,'PRINT=SCF=CORE-LAGRANGIAN',.false.,' ')) then
         write (iout,88)
88       format (5x,'MO core lagrangian')
         call matout(eigval,nbf,1,nbf,1,iout)
      end if
c
      call iosys('write real "scf mo core lagrangian" to rwf',
     $     nbf,eigval,0,' ')
c
c
c     ----- rotate degenerate eigenvectors -----
c
c..bhl      call enforc(nbf,nnp,s,smhalf,eigval,c,t1,t2)
c
c
      return
      end

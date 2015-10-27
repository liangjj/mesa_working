*deck @(#)unique.f	5.1  11/6/94
      subroutine unique(c,eigval,f,jmat,kmat,values,d,u,t1,t2,h,
     $                  nbf,nnp,nshell,scfnum,ntriang,iout,nmat,
     $                  shlmin,shlmax,ops,s,smhalf,fcoef)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),eigval(nbf),f(nnp),jmat(nnp,nmat),kmat(nnp,nmat)
      real*8 values(nnp,ntriang),d(nnp,nmat),u(nbf,nbf),t1(nbf,nbf)
      real*8 t2(nbf,nbf,nmat),h(nnp),s(nnp),smhalf(nnp),fcoef(nshell)
      integer shlmin(nshell),shlmax(nshell)
      character*(*) ops
      logical logkey
c
c     ----- form the density matrix or matrices -----
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
c      do 2 i=1,nnp
c         f(i)=f(i)+2*jmat(i,1)-kmat(i,1)
c    2 continue
c      if (scfnum.eq.1) then
c         do 3 i=1,nnp
c            f(i)=f(i)+jmat(i,2)-0.5d+00*kmat(i,2)
c    3    continue
c      end if
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
c     ----- and diagonalize this matrix, transforming the scf vectors -
c
      call degrsp(nbf,nnp,f,eigval,1,u,t1,t2)
c
      if (logkey(ops,'print=scf=core-lagrangian',.false.,' ')) then
         write (iout,8)
 8       format (5x,'mo core lagrangian')
         call matout(eigval,nbf,1,nbf,1,iout)
      end if
c
      call iosys('write real "scf mo core lagrangian" to rwf',
     $     nbf,eigval,0,' ')
c
      call ebc(t1,c,u,nbf,nbf,nbf)
      call vmove(c,t1,nbf**2)
c
c     ----- rotate degenerate eigenvectors -----
c
      call enforc(nbf,nnp,s,smhalf,eigval,c,t1,t2)
c
c
      return
      end

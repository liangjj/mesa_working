*deck @(#)symrsp.f	3.1  11/20/92
      subroutine symrsp(triang,u,t2,t3,t4,t5,eigval,temp,kept,nbf,nnp,
     $                  nshell,shlmin,shlmax,nirrep,numso,lambda,
     $                  orbsym,occsym,bftoso,c)

      implicit integer(a-z)
c
c     nirrep   number of irreducible representations
c     numso    number of salcs in each irreducible representation
c     lambda   dimension of each irrep
c     orbsym   the irreducible representation of the mo's.
c     occsym   the number of occupied orbitals in each irrep and shell.
c
      integer numso(nirrep),lambda(nirrep),orbsym(nbf),temp(nbf)
      integer occsym(nirrep,nshell)
      integer shlmin(nshell),shlmax(nshell),kept(nshell)
      real*8 triang(nnp),u(nbf,nbf),t2(nbf,nbf),t3(nbf,nbf),eigval(nbf)
      real*8 t4(nbf,nbf),t5(nbf,nbf),bftoso(nbf,nbf),c(nbf,nbf)
c
      real*8 big
      logical debug
      common/io/inp,iout
c
      parameter (debug=.false.)
 1000 format(1x,'eigenvectors from symrsp')
c
c
      call trtosq(t2,triang,nbf,nnp)
      salc=1
      call izero(kept,nshell)
      do 100 irrep=1,nirrep
         nsf=numso(irrep)*lambda(irrep)
         if(nsf.ne.0) then
c
c           transform fock matrix to the salc basis
            call ebc(t3,t2,bftoso(1,salc),nbf,nbf,nsf)
            call ebtc(t5,bftoso(1,salc),t3,nsf,nbf,nsf)
c
c           diagonalize
            ltri=nsf*(nsf+1)/2
            call sqtotr(triang,t5,nsf,ltri)
            call degrsp(nsf,ltri,triang,t4,1,u,t3,t5)
c
c           transform back to the original basis.
            call ebc(t3,bftoso(1,salc),u,nbf,nsf,nsf)
c           keep the correct number of vectors of this symmetry in each shell.
            mo=0
            do 80 shell=1,nshell
               target=shlmin(shell)-1+kept(shell)
               do 70 i=1,occsym(irrep,shell)
                  mo=mo+1
                  call vmove(c(1,target+i),t3(1,mo),nbf)
                  eigval(target+i)=t4(mo,1)
                  orbsym(target+i)=irrep
                  kept(shell)=kept(shell)+1
   70             continue
   80          continue
            salc=salc+nsf
         end if
  100 continue
c
c
c
c     ---- reorder the eigenvectors and eigenvalues according
c          to ascending eigenvalue in each shell
c          return vectors in c, eigenvalues in eigval,
c          and symmetries in orbsym
      call vmove(t3,c,nbf*nbf)
      call vmove(triang,eigval,nbf)
      call vmove(temp,orbsym,nbf)
      big=1.0d+44
c
      do 200 shell=1,nshell
         do 190 i=shlmin(shell),shlmax(shell)
            offset=shlmin(shell)-1
            mineig=offset+ismin(shlmax(shell)-offset,triang(offset+1),1)
            call vmove(c(1,i),t3(1,mineig),nbf)
            eigval(i)=triang(mineig)
            orbsym(i)=temp(mineig)
            triang(mineig)=big
 190     continue
 200  continue
c
c
      return
      end

*deck @(#)pseud.f	5.1  11/6/94
      subroutine pseud(nbf,nnp,nshell,ncoul,nexch,hao,h,jmat,kmat,
     $     c,lag,u,eigval,f,alpha,beta,shlmin,shlmax,t1,t2,ops)
c
c***begin prologue     pseud
c***date written       870720   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           scf orbitals, pseudocanonical orbitals
c***author             saxe, paul (lanl)
c***source             @(#)pseud.f	5.1   11/6/94
c
c***purpose            to form a psuedocanonical representation of the
c                      scf orbitals.
c
c***description
c     this routine forms the diagonal blocks of the lagrangian matrix
c     and diagonalizes them to form the transformation matrix to
c     give a pseudocanonical representation for the scf orbitals.
c
c***references
c
c***routines called    (none)
c
c***end prologue       pseud
c
      implicit integer (a-z)
c
      character*(*) ops
      integer shlmin(nshell)
      integer shlmax(nshell)
      logical logkey
      real*8 jmat(nnp,ncoul)
      real*8 kmat(nnp,ncoul)
      real*8 t1(nbf,nbf)
      real*8 t2(nbf,nbf)
      real*8 hao(nnp)
      real*8 h(nnp)
      real*8 c(nbf,nbf)
      real*8 lag(nnp)
      real*8 u(nbf,nbf)
      real*8 eigval(nbf)
      real*8 f(nshell)
      real*8 alpha(nshell,nshell)
      real*8 beta(nshell,nshell)
      real*8 fc
      real*8 ac
      real*8 bc
      real*8 small
c
      parameter (small=0.1d+00)
c
      common /io/ inp,iout
c
c     ----- statement functions -----
c
      offset(i,j)=i*(i-1)/2+j
c
c
c..bhl
      do 1 i=1,nnp
         lag(i)=hao(i)+2.0d+00*jmat(i,1)-kmat(i,1)
    1 continue
      call iosys('write real "scf ao fock operator" to rwf',
     $            nnp,lag,0,' ')
c..bhl
c     ----- transform the one-electron integrals to the mo basis -----
c
      call trtosq(t2,hao,nbf,nnp)
      call ebc(t1,t2,c,nbf,nbf,nbf)
      call ebtc(t2,c,t1,nbf,nbf,nbf)
      call sqtotr(h,t2,nbf,nnp)
c
c     ----- transform the j and k matrices to the mo basis -----
c
      do 10 i=1,ncoul
         call trtosq(t2,jmat(1,i),nbf,nnp)
         call ebc(t1,t2,c,nbf,nbf,nbf)
         call ebtc(t2,c,t1,nbf,nbf,nbf)
         call sqtotr(jmat(1,i),t2,nbf,nnp)
 10   continue
      do 20 i=1,nexch
         call trtosq(t2,kmat(1,i),nbf,nnp)
         call ebc(t1,t2,c,nbf,nbf,nbf)
         call ebtc(t2,c,t1,nbf,nbf,nbf)
         call sqtotr(kmat(1,i),t2,nbf,nnp)
 20   continue
c
c     ----- form the diagonal blocks of the lagrangian matrix -----
c               (use closed-shell operator for virtuals)
c
      call rzero(lag,nnp)
c
      do 100 shell=1,nshell
c
c        ----- the one-electron term -----
c
         if (shell.ne.nshell) then
cps            fc=f(shell)
            fc=1.0d+00
         else
            fc=1.0d+00
         end if
         do 40 i=shlmin(shell),shlmax(shell)
            do 30 j=shlmin(shell),i
               lag(offset(i,j))=fc*h(offset(i,j))
 30         continue
 40      continue
c
c        ----- the two-electron part -----
c
         do 80 lshell=1,nshell-1
c
            if (shell.ne.nshell) then
               ac=alpha(shell,lshell)/f(shell)
               bc=beta(shell,lshell)/f(shell)
            else
               ac=2.0d+00*f(lshell)
               bc=-1.0d+00*f(lshell)
            end if
c
            do 70 i=shlmin(shell),shlmax(shell)
               ia=offset(i,0)
               do 60 j=shlmin(shell),i
                  ij=ia+j
                  lag(ij)=lag(ij)+ac*jmat(ij,lshell)+bc*kmat(ij,lshell)
 60            continue
 70         continue
 80      continue
 100  continue
c
      if (logkey(ops,'scf=print=pseudolagrangian',.false.,' ')) then
         write (iout,105)
 105     format (/,t5,'the pseudolagrangian matrix')
         call print(lag,nnp,nbf,iout)
      end if
c
c     ----- now diagonalize this blocked matrix -----
c
      call degrsp(nbf,nnp,lag,eigval,1,u,t1,t2)
      call vmove(t1,eigval,nbf)
      call vmove(t2,u,nbf**2)
c
c     ----- make sure the transformation matrix does not mix up blocks
c
      do 200 shell=1,nshell
         n=shlmin(shell)-1
         do 190 i=1,nbf
            do 180 j=shlmin(shell),shlmax(shell)
               if (abs(t2(j,i)).gt.small) then
                  n=n+1
                  if (n.gt.shlmax(shell)) then
                     call lnkerr('found too many vectors')
                  end if
                  eigval(n)=t1(i,1)
                  call vmove(u(1,n),t2(1,i),nbf)
                  go to 190
               end if
 180        continue
 190     continue
         if (n.ne.shlmax(shell)) then
            call lnkerr('not enough vectors found')
         end if
 200  continue
c
c     ----- and transform the scf vector to this basis -----
c
      call ebc(t1,c,u,nbf,nbf,nbf)
      call vmove(c,t1,nbf**2)
c
c     ----- transform the one-electron integrals to the new mo basis --
c
      call trtosq(t2,h,nbf,nnp)
      call ebc(t1,t2,u,nbf,nbf,nbf)
      call ebtc(t2,u,t1,nbf,nbf,nbf)
      call sqtotr(h,t2,nbf,nnp)
c
c     ----- transform the j and k matrices to the new mo basis -----
c
      do 110 i=1,ncoul
         call trtosq(t2,jmat(1,i),nbf,nnp)
         call ebc(t1,t2,u,nbf,nbf,nbf)
         call ebtc(t2,u,t1,nbf,nbf,nbf)
         call sqtotr(jmat(1,i),t2,nbf,nnp)
 110  continue
      do 120 i=1,nexch
         call trtosq(t2,kmat(1,i),nbf,nnp)
         call ebc(t1,t2,u,nbf,nbf,nbf)
         call ebtc(t2,u,t1,nbf,nbf,nbf)
         call sqtotr(kmat(1,i),t2,nbf,nnp)
 120  continue
c
c
      return
      end

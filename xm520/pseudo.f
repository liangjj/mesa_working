*deck  @(#)pseudo.f	5.1 11/6/94
      subroutine pseudo(nbf,nnp,nshell,ncoul,nexch,hao,h,jmat,kmat,
     $                  c,lag,u,eigval,f,alpha,beta,shlmin,shlmax,
     $                  t1,t2,ops)
c***begin prologue     pseudo.f
c***date written       930610  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin,richard and saxe,paul (lanl) 
c***source             @(#)pseudo.f	5.1   11/6/94
c***purpose            to transform the one-electron,j,and k matrices
c                      to a pseudo-canonical orbital representation. 
c***description
c     this routine transforms the one-electron, coulomb and
c     exchange matrices to a pseudo-canonical representation.
c
c     the pseudo-canonical representation is generated by forming
c     and diagonalizing the diagonal shell blocks of the lagrangian
c     matrix.  this "pre-diagonalization" mixes occupied orbitals
c     within shell blocks only, and can sometimes accelerate
c     convergence. 
c     
c***references
c                      m.page and j.w.mciver,jr., j.chem.phys. 79,4985(1983).
c
c***routines called
c
c***end prologue       pseudo.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,nshell,ncoul,nexch
c     --- input arrays (unmodified) ---
      character*(*) ops
      integer shlmin(nshell),shlmax(nshell)
      real*8 f(nshell),alpha(nshell,nshell),beta(nshell,nshell)
      real*8 hao(nnp)
c     --- input arrays (scratch) ---
      real*8 lag(nnp),u(nbf,nbf),eigval(nbf)
      real*8 t1(nbf,nbf),t2(nbf,nbf)
c     --- output arrays ---
      real*8 h(nnp),jmat(nnp,nshell),kmat(nnp,nshell)
      real*8 c(nbf,nbf)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      logical logkey
      integer inp,iout
      integer i,j,ia,ij,shell,lshell,n
      integer offset
      real*8 fc
      real*8 ac
      real*8 bc
      real*8 small
c
      parameter (small=0.1d+00)
c
      common /io/ inp,iout
c
c     --- statement functions ---
      offset(i,j)=i*(i-1)/2+j
c
c     --- transform the one-electron,j and k matrices to the mo basis ---
      call tomo(nbf,nnp,nshell,ncoul,nexch,hao,h,jmat,kmat,c,t1,t2)
c
c     --- form the diagonal blocks of the lagrangian matrix ---
c            (use closed-shell operator for virtuals)
      call rzero(lag,nnp)
      do 100 shell=1,nshell
c        --- the one-electron term ---
         if (shell.ne.nshell) then
cps         fc=f(shell)
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
c        --- the two-electron part ---
         do 80 lshell=1,nshell-1
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
c     --- now diagonalize this blocked matrix ---
      call degrsp(nbf,nnp,lag,eigval,1,u,t1,t2)
      call vmove(t1,eigval,nbf)
      call vmove(t2,u,nbf**2)
c
c     --- make sure the transformation matrix does not reorder
c         the constitution of blocks
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
c     --- and transform the scf vector to this basis ---
      call ebc(t1,c,u,nbf,nbf,nbf)
      call vmove(c,t1,nbf**2)
c
c     --- transform the one-electron,j and k matrices to 
c         the new mo basis ---
      call vmove(lag,h,nnp)
      call tomo(nbf,nnp,nshell,ncoul,nexch,lag,h,jmat,kmat,u,t1,t2)
c
c
      return
      end

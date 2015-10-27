*deck @(#)ham.f	5.1  11/6/94
      subroutine ham(fockc,fock,c,h,jmat,kmat,t1,t2,eig,work1,work2,
     1               eigsav,d,ops,num,occmo,nnp,ndmat,homo,prnt)
c
c***begin prologue      ham
c***date written        861110  (yymmdd)
c***revision date       yymmdd  (yymmdd)
c***keywords
c***author              barry schneider (lanl) and byron lengsfield (llnl)
c***purpose             form one electron operators of fock type for ivo
c***                    calculations then diagonalize to get orbitals
c***description
c***reference
c***routines called      none
c***end prologue
c
      implicit integer(a-z)
      parameter (maxnbf=2000)
      logical logkey,usesym,prnt
      character*16 bflabl(maxnbf)
      character *(*) ops
      character itoc*4
      real*8 c(num,num), fock(num,num), h(nnp), jmat(nnp,ndmat)
      real*8 kmat(nnp,ndmat), t1(num*num), t2(num,num), fockc(*)
      real*8 eig(num),work1(num),work2(num),eigsav(num),d(nnp,*)
      real*8 f,aij,bij,fn,fd,an,ad,bn,bd
      real*8 eclos,eopen,ecore,enuc,sdot
      common /shell/  noshel(10), orbshel(50,10), nvshel(10),
     1 vshel(200,10), f(10), aij(10,10), bij(10,10),
     2 fn(10),fd(10),an(10),ad(10),bn(10),bd(10)
      common /io/ inp, iout
c
c     ----- loop over the shell hamiltonia -----
c
      call iosys('read character "basis function labels" from rwf',
     $            len(bflabl(1))*num,0,0,bflabl)
c
c
      pointer=occmo+1
      dimi=noshel(ndmat+1)
      call rzero (fock,num*num)
      call rzero (fockc,nnp)
c
c
c     ----- put one electron term in fock ------
c
      ik=0
      do 10 j=1,num
         do 10 k=1,j
            ik=ik+1
            fock(j,k)=fock(j,k)+f(1)*h(ik)
   10 continue
c
c     ----- add in the coulomb and exchange contribution -----
c
      do 30 j=1,ndmat
         kl=0
         do 20 k=1,num
            do 20 l=1,k
               kl=kl+1
               fock(k,l)=fock(k,l)
     $                  +aij(j,1)*jmat(kl,j)+bij(j,1)*kmat(kl,j)
   20    continue
   30 continue
c
c     build the closed shell fock operator
c
      if(ndmat.ne.1) then
         do 35 k=1,nnp
            fockc(k)=h(k)+2.d0*jmat(k,1)-kmat(k,1)
  35     continue
      end if
c
      do 40 j=1,num
         do 40 k=1,j
            fock(k,j)=fock(j,k)
   40 continue
c
c      call iosys ('write real "ivo ao ham'//itoc(i)//'" to rwf',num*num,
c     1         fock,0,' ')
c
      if (index(ops,'print_ao_fock_matrix').ne.0) then
         write (iout,70) 1
         call matout (fock,num,num,num,num,iout)
      endif
c
c     ----- the ith fock operator has been formed in the ao -----
c     -----                  representation -----
c     ----- now transform it to the mo basis -----
c
      call ebc (t2,fock,c(1,pointer),num,num,dimi)
      call ebtc (t1,c(1,pointer),t2,dimi,num,dimi)
c     call iosys ('write real "ivo mo ham'//itoc(i)//'" to rwf',
c     1            dimi*dimi,t1,0,' ')
c
      if (index(ops,'print_mo_fock_matrix').ne.0) then
         write (iout,80) i
         call matout (t1,dimi,dimi,dimi,dimi,iout)
      endif
c
c     -----              diagonalize              -----
c
c     temporary change noted below.
      usesym=.false.
      if (usesym) then
c        ----- if symmetry is turned on -----
c         not yet implemented
c        call dagsym
         call lnkerr(' symmetry not yet implemented in m502')
      else
c        ----- if symmetry is turned off -----
         nnpdim=dimi*(dimi+1)/2
         call sqtotr(t2,t1,dimi,nnpdim)
         call degrsp(dimi,nnpdim,t2,eig,1,t1,work1,work2)
      endif
c
c     ----- transform vectors back to ao representation -----
c     -----     and put back into c                     -----
      call ebc (t2,c(1,pointer),t1,num,dimi,dimi)
c
      call vmove(eigsav(pointer),eig,dimi)
c
      call vmove (c(1,pointer),t2,num*dimi)
c
c     compute the total energy
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $           1,enuc,0,' ')
c
      if(ndmat.eq.1) then
         call trtosq(t2,h,num,nnp)
         call trtosq(t1,d(1,1),num,nnp)
         eopen=sdot(num*num,t1,1,t2,1)
         eclos=0.d0
      else
         call trtosq(t2,h,num,nnp)
         call trtosq(t1,d(1,1),num,nnp)
         eclos=sdot(num*num,t1,1,t2,1)
         call trtosq(t2,fockc,num,nnp)
         eclos=eclos+sdot(num*num,t1,1,t2,1)
         call trtosq(t1,d(1,2),num,nnp)
         eopen=sdot(num*num,t1,1,t2,1)
      end if
c
      mm=min(5,dimi)
      ecore=eopen+eclos+enuc
      if(prnt) then
         write(iout,91) enuc,eclos,eopen,ecore
         write(iout,92) (i,eig(i),eig(i)+ecore,i=1,mm)
      endif
c
c
      if(homo.ne.0) then
         if(prnt) then
            write(iout,93)
         endif
         call flip(c,num,occmo,homo)
      end if
c
c
      if(logkey(ops,'print=ivo',.false.,' ')) then
         call wvec(c,eigsav,num,num,bflabl,' ')
      endif
c
      call iosys('write real "ivo vector" to rwf',num*num,c,0,' ')
      call iosys('write real "ivo energies" to rwf',num,eigsav,0,' ')
c
      if(prnt) write(iout,100)
      call iosys('write real "scf vector" to rwf',num*num,c,0,' ')
      call iosys('write real "orbital energies" to rwf',
     #            num,eigsav,0,' ')
c
      return
c
   70 format (/,10x,i3,'th',1x,'ao fock matrix')
   80 format (/,10x,i3,'th',1x,'mo fock matrix')
   90 format (/,15x,'ivo eigenvalues for',1x,i2,1x,'shell')
   91 format(5x,'*** ivo energy evaluation ***',/,
     #       5x,'nuclear repulsion energy     ',f18.9,/,
     #       5x,'closed  shell     energy     ',f18.9,/,
     #       5x,'open    shell     energy     ',f18.9,/,
     #       5x,'total   core      energy     ',f18.9)
   92 format(/,'        state        ivo eigenvalue ',
     #          '       total energy   ',
     # 5(/,5x,i5,5x,2(2x,f18.9)),/)
   93 format(5x,'restoring original mo order')
  100 format(5x,'overwriting scf vector with ivos')
      end

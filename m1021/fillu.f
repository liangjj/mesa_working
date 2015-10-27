*deck @(#)fillu.f	5.1  11/6/94
      subroutine fillu(values,nnp,num,ntriang,isq,ds,u,c,dclag,
     #                 pt,nindep,f,nshell,orbshl,
     #                 minshl,maxshl,numshl,clag,t1,t2,nder,udep,ops,
     #                 nco,fj,fk,intj,uk,nconum)
c
c***begin prologue     fillu
c***date written       861211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           u vectors, coupled perturbed hartree-fock, cphf
c***author             saxe, paul (lanl)
c***source             @(#)fillu.f	5.1   11/6/94
c***purpose            form the right hand sides (b0's) of the cphf equations
c***description
c
c  this routine uses the following equation to fill in the dependent
c  solutions to the cphf equations. nb that this implies that the hf
c  orbitals have been rotated to diagonalize a core fock matrix defined
c  as follows:
c                          occ
c       clag(ij) = h(ij) + sum f(k) * {2j(ij) - k(ij)}
c                           k
c
c
c
c                      1           ind pairs
c      u(ij,a) = --------------- { sum sum   u(kl,a) * ( 4[ij;kl] -
c                clag(j)-clag(i)    k > l
c
c                                     [il;jk] - [ik;jl] )
c
c           occ occ
c         - sum sum ds(kl,a) * f(k) * ( 4[ij;kl] - [il;jk] - [ik;jl] )
c            k > l
c
c           occ
c         - sum ds(kk,a) * f(k) * ( 2[ij;kk] - [ik;jk] ) + dclag(ij,a) -
c            k
c
c                                          clag(j) * ds(ij,a) }
c
c***references         "unified theoretical treatment of analytic first
c                       and second derivatives in open-shell hartree-
c                       fock theory", y. osamura, y. yamaguchi, p. saxe,
c                       m. a. vincent, j. f. gaw and h. f. schaefer iii,
c                       chemical physics 72 (1983) 131-139.
c
c***routines called    iosys
c
c***end prologue       fillu
c
      implicit integer (a-z)
c
      character*(*) ops
      logical logkey
      real*8 values(nnp,ntriang)
      real*8 isq(num,num)
      real*8 f(nshell)
      real*8 clag(num)
      real*8 ds(nnp,nder),udep(nindep,nder),c(num,num)
      real*8 t1(num,num),t2(num,num),dclag(num,num,nder)
      real*8 u(num,num,nder)
      real*8 fac,fshell
      real*8 intj(nconum,ntriang)
      real*8 fj(nnp,nder),fk(num,nder,num)
      real*8 uk(num,nder,nco)
      integer pt(nnp),orbshl(num)
      integer minshl(nshell),maxshl(nshell),numshl(nshell)
c
      common /io/ inp,iout
c
      numnum=num*num
c
c     ----- read in the dependent elements of u and put into full u ---
c
      call iosys('read real "independent cphf solutions" from rwf',
     $     nindep*nder,udep,0,' ')
c
      call rzero(u,num*num*nder)
      do 5 ishell=2,nshell
         do 4 jshell=1,ishell-1
            do 3 i=minshl(ishell),maxshl(ishell)
               ia=i*(i-1)/2
               do 2 j=minshl(jshell),maxshl(jshell)
                  ij=ia+j
                  pij=pt(ij)
                  do 1 der=1,nder
                     u(i,j,der)=udep(pij,der)
 1                continue
 2             continue
 3          continue
 4       continue
 5    continue
c
c     ----- create the array of orbital shells (shells) -----
c
      do 7 ishell=1,nshell
         do 6 i=minshl(ishell),maxshl(ishell)
            orbshl(i)=ishell
 6       continue
 7    continue
c
c     ----- read in and transform the derivative overlap integrals -----
c
      call iosys('read real "scf vector" from rwf',num**2,c,0,' ')
      call iosys('read real "ao derivative overlap integrals" '//
     $     'from dints',nnp*nder,ds,0,' ')
c
      do 9 der=1,nder
         call trtosq(t1,ds(1,der),num,nnp)
         call ebc(t2,t1,c,num,num,num)
         call ebtc(t1,c,t2,num,num,num)
         call sqtotr(ds(1,der),t1,num,nnp)
 9    continue
c
c     ----- u(i,j)+ds(ij)+u(j,i) = 0 -----
c

      ii=0
      do 927 i=1,num
      ii=ii+i
        do 926 der=1,nder
         u(i,i,der)=-.5d0*ds(ii,der)
 926    continue
 927  continue
c
      do 130 i=2,num
         ia=i*(i-1)/2
         do 129 j=1,i-1
            ij=ia+j
            do 128 der=1,nder
               u(j,i,der)=-u(i,j,der)-ds(ij,der)
 128        continue
 129     continue
 130  continue
c
c  scale the u-matrix   u(j,i) --> u(j,i)*f(i)
c
      do 820 i=1,nco
        ishl=orbshl(i)
        fshell=f(ishl)
       do 821 der=1,nder
        do 822 j=1,num
         u(j,i,der)=u(j,i,der)*fshell
 822    continue
        do 823 j=1,num
         uk(j,der,i)=u(j,i,der)
 823    continue
 821   continue
 820  continue
c
c     ----- read in the derivative lagrangians and form
c
      call iosys('read real "mo derivative core lagrangian" from rwf',
     #            num**2*nder,dclag,0,' ')
c
c
c     ----- get the core lagrangian ------
c
      call iosys('read real "scf mo core lagrangian" from rwf',
     $     num,clag,0,' ')
c
c     ----- rewind the mo integrals -----
c
      call iosys('rewind "mo two-electron integrals" on tints',
     $            0,0,0,' ')
c
c     ----- read through integrals -----
c
      n=0
c
      call rzero(fk,num*num*nder)
c
      maxkl=0
      kl=0
c
      do 110 k=1,nco
         do 100 l=1,k
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               nread=maxkl-minkl+1
               lnread=nread*nnp
               call iosys('read real "mo two-electron integrals" '//
     $              'from tints without rewinding',lnread,values,0,' ')
            end if
c
            call trtosq(t1,values(1,kl-minkl+1),num,nnp)
            call scopy(nconum,t1,1,intj(1,kl-minkl+1),1)
c
            if(kl.eq.maxkl) then
c
c           ----- work on coulomb part of u -----
c
c               fj(kl,der)=intj(ml,kl)*u(ml,der)
c
       call mxma(intj,nconum,1,u,1,numnum,fj(minkl,1),1,nnp,
     $           nread,nconum,nder)
            end if
c
c           ----- work on exchange part of u -----
c
c               fk(i,der,k)=fk(i,der,k)+int(kl,i,m)*uk(m,der,l)
c
         if(k.ne.l) then
             call apbc(fk(1,1,k),t1,uk(1,1,l),num,num,nder)
             call apbc(fk(1,1,l),t1,uk(1,1,k),num,num,nder)
         else
             call apbc(fk(1,1,k),t1,uk(1,1,l),num,num,nder)
         end if
c
 100     continue
 110   continue
c
c
      do 210 k=nco+1,num
         do 200 l=1,nco
c
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               nread=maxkl-minkl+1
               lnread=nread*nnp
               call iosys('read real "mo two-electron integrals" '//
     $              'from tints without rewinding',lnread,values,0,' ')
            end if
c
            call trtosq(t1,values(1,kl-minkl+1),num,nnp)
            call scopy(nconum,t1,1,intj(1,kl-minkl+1),1)
c
            if(kl.eq.maxkl) then
c
c           ----- work on coulomb part of u -----
c
c               fj(kl,der)=intj(ml,kl)*u(ml,der)
c
      call mxma(intj,nconum,1,u,1,numnum,fj(minkl,1),1,nnp,
     $                 nread,nconum,nder)
            end if
c
c           ----- work on exchange part of u -----
c
c               fk(i,der,k)=fk(i,der,k)+int(kl,i,m)*uk(m,der,l)
c
             call apbc(fk(1,1,k),t1,uk(1,1,l),num,num,nder)
c
 200     continue
c
         do 205 l=nco+1,k
c
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               nread=maxkl-minkl+1
               lnread=nread*nnp
               call iosys('read real "mo two-electron integrals" '//
     $              'from tints without rewinding',lnread,values,0,' ')
            end if
c
            call trtosq(t1,values(1,kl-minkl+1),num,nnp)
            call scopy(nconum,t1,1,intj(1,kl-minkl+1),1)
c
            if(kl.eq.maxkl) then
c
c           ----- work on coulomb part of u -----
c
c               fj(kl,der)=intj(ml,kl)*u(ml,der)
c
             call mxma(intj,nconum,1,u,1,numnum,fj(minkl,1),1,nnp,
     $                 nread,nconum,nder)
            end if
c
 205     continue
 210   continue
c
c  un-do scaling of u-matrix   u(j,i)*f(i) --->  u(j,i)
c
      do 701 i=1,nco
       ishl=orbshl(i)
       fac=1.d0/f(ishl)
       do 702 der=1,nder
        do 703 j=1,num
         u(j,i,der)=u(j,i,der)*fac
 703    continue
 702   continue
 701  continue
c
c     ----- add in the lagrangian term:
c
      do 126 shell=1,nshell
         do 125 i=minshl(shell),maxshl(shell)
            ia=i*(i-1)/2
            do 124 j=minshl(shell),i
               ij=ia+j
               fac=clag(j)-clag(i)
               if (abs(fac).lt.1.0d-06) then
                  if (i.ne.j) write (iout,121) i,j
 121           format (5x,'warning: orbitals ',i3,' and ',i3,
     $                 'are degenerate')
                  do 122 der=1,nder
                     u(i,j,der)=-0.5d+00*ds(ij,der)
 122              continue
               else
                  do 123 der=1,nder
                   u(i,j,der)=(dclag(i,j,der)-ds(ij,der)*clag(j)+
     $              4.d0*fj(ij,der)-fk(j,der,i)-fk(i,der,j))/fac
 123              continue
                  do 223 der=1,nder
                     u(j,i,der)=-u(i,j,der)-ds(ij,der)
 223              continue
               end if
 124        continue
 125     continue
 126  continue
c
c      write(iout,*)' dclag ndf=1 '
c      call matout(dclag,num,num,num,num,iout)
c
      do 991 i=1,num
       do 992 j=1,num
       dclag(i,j,1)=fk(i,1,j)
 992   continue
 991  continue
c
c      write(iout,*)' fk ndf=1 '
c      call matout(dclag,num,num,num,num,iout)
c      write(iout,*)' fj ndf=1 '
c      call print(fj,nnp,num,iout)
c      write(iout,*)' ds ndf=1 '
c      call print(ds,nnp,num,iout)
c23456
      if (logkey(ops,'print=cphf=solutions',.false.,' ')) then
         do 132 der=1,nder
            write (iout,131) der
 131        format(5x,'cphf solution vector, derivative ',i3)
            call matout(u(1,1,der),num,num,num,num,iout)
 132     continue
      end if
c
c     ----- and put the solution vectors out -----
c
      call iosys('write real "cphf solutions" to rwf',
     $     num**2*nder,u,0,' ')
c
c
      return
      end

*deck @(#)dumbu.f	5.1  11/6/94
      subroutine dumbu(values,nnp,num,ntriang,isq,ds,u,c,dclag,
     #                 pt,nindep,f,nshell,orbshl,
     #                 minshl,maxshl,numshl,clag,t1,t2,nder,udep,ops,
     #                 nco,fj,fk,intj,uk,ints)
c
c***begin prologue     dumbu
c***date written       861211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           u vectors, coupled perturbed hartree-fock, cphf
c***author             saxe, paul (lanl)
c***source             @(#)dumbu.f	5.1   11/6/94
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
c         - sum ds(kk,a) * f(k) * ( 2[ij;kk] - [ik;jk] + dclag(ij,a) -
c            k
c
c                                          clag(j) * ds(ij,a) ) }
c
c***references         "unified theoretical treatment of analytic first
c                       and second derivatives in open-shell hartree-
c                       fock theory", y. osamura, y. yamaguchi, p. saxe,
c                       m. a. vincent, j. f. gaw and h. f. schaefer iii,
c                       chemical physics 72 (1983) 131-139.
c
c***routines called    iosys
c
c***end prologue       dumbu
c
      implicit integer (a-z)
c
      character*(*) ops
      real*8 values(nnp,ntriang)
      real*8 isq(num,num)
      real*8 f(nshell)
      real*8 clag(num)
      real*8 ds(nnp,nder),udep(nindep,nder),c(num,num)
      real*8 t1(num,num),t2(num,num),dclag(num,num,nder)
      real*8 u(num,num,nder)
      real*8 fac
      real*8 intj(num,nco,ntriang)
      real*8 fj(nnp,nder),fk(num,num,nder)
      real*8 uk(num,nder,nco),ints(num,num,num,num)
      integer pt(nnp),orbshl(num)
      integer minshl(nshell),maxshl(nshell),numshl(nshell)
c
      common /io/ inp,iout
c
      numnum=num*num
      nconum=nco*num
      write(iout,*)'  nco  ',nco
      write(iout,*)'  f  '
      write(iout,*)(f(i),i=1,nshell)
      write(iout,*)'   numshl   '
      write(iout,*)(numshl(i),i=1,nshell)
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
      write(iout,*)'  u-input  '
      call matout(u,num,num,num,num,iout)
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
     $     'from rdints',nnp*nder,ds,0,' ')
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
      write(iout,*)'  u-input sq with ds  '
      call matout(u,num,num,num,num,iout)
c
c     ----- read in the derivative lagrangians and form
c
      call iosys('read real "mo derivative core lagrangian" from rwf',
     #            num**2*nder,dclag,0,' ')
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
      call iosys('read real "mo two-electron integrals" from tints',
     $            nnp**2,values,0,' ')
c
      kl=0
      do 3001 k=1,num
        do 3002 l=1,k
          kl=kl+1
          call trtosq(ints(1,1,k,l),values(1,kl),num,nnp)
          call trtosq(ints(1,1,l,k),values(1,kl),num,nnp)
3002    continue
3001  continue
c
c
c           ----- work on coulomb part of u -----
c
c               fj(ij,der)=intj(ml,ij)*u(ml,der)
c
           call rzero(fj,nnp*nder)
c
           kl=0
           do 1001 k=1,num
            do 1002 l=1,k
            kl=kl+1
             do 1003 i=1,nco
              do 1004 j=1,num
               do 1005 der=1,nder
                fj(kl,der)=fj(kl,der)+ints(k,l,j,i)*u(j,i,der)
 1005          continue
 1004         continue
 1003        continue
 1002       continue
 1001      continue
c
cc-cc-cc
           call rzero(fk,num*num*nder)
cc-cc-cc
c
           do 1501 k=1,nco
            do 1502 l=1,k
             do 1503 i=1,num
              do 1505 j=1,num
               do 1506 der=1,nder
                fk(k,i,der)=fk(k,i,der)+ints(k,l,i,j)*u(j,l,der)
 1506          continue
 1505         continue
             if(k.ne.l) then
              do 1504 j=1,num
               do 1508 der=1,nder
                fk(l,i,der)=fk(l,i,der)+ints(k,l,i,j)*u(j,k,der)
 1508          continue
 1504         continue
             end if
 1503        continue
 1502       continue
 1501      continue
c
           do 2001 k=nco+1,num
            do 2002 l=1,nco
             do 2003 i=1,num
              do 2004 j=1,num
               do 2005 der=1,nder
                fk(k,i,der)=fk(k,i,der)+ints(k,l,i,j)*u(j,l,der)
 2005          continue
 2004         continue
 2003        continue
 2002       continue
 2001      continue
c
c
c
c      write(iout,*)' dclag ndf=1 '
c      call matout(dclag,num,num,num,num,iout)
c      write(iout,*)' fk ndf=1 '
c      call matout(fk,num,num,num,num,iout)
c      write(iout,*)' fj ndf=1 '
c      call print(fj,nnp,num,iout)
c      write(iout,*)' ds ndf=1 '
c      call print(ds,nnp,num,iout)
c
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
     $               4.d0*fj(ij,der)-fk(i,j,der)-fk(j,i,der))/fac
 123              continue
                  do 223 der=1,nder
                     u(j,i,der)=-u(i,j,der)-ds(ij,der)
 223              continue
               end if
 124        continue
 125     continue
 126  continue
c
c
c         do 1132 der=1,3
c            write (iout,1131) der
c1131        format(5x,'cphf solution vector, derivative ',i3)
c            call matout(u(1,1,der),num,num,num,num,iout)
c1132     continue
c
c     ----- and put the solution vectors out -----
c
      call iosys('write real "cphf solutions" to rwf',
     $     num**2*nder,u,0,' ')
c
c
      return
      end

*deck @(#)mcscfu.f	5.1  11/6/94
      subroutine mcscfu(values,nnp,num,ntriang,isq,ds,u,c,dclag,
     #                 clag,t1,t2,nder,ops,
     #                 nco,fj,fk,intj,uk,nconum,nao)
c
c***begin prologue     mcscfu
c***date written       861211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           u vectors, coupled perturbed mcscf, cpmcscf
c***author             lengsfield, byron (brl)
c***source             @(#)mcscfu.f	5.1   11/6/94
c***purpose            form the core-core and virtual-virtual portions of
c                      the u-matrix so one can truncate core or virtual
c                      orbitals in a multi-reference ci gradient calculation
c***description
c
c  this routine uses the following equation to fill in the dependent
c  solutions to the cpmcscf equations. nb that this implies that the mcscf
c  orbitals have been rotated to diagonalize a core fock matrix defined
c  as follows:
c                          occ
c       clag(ij) = h(ij) + sum  {2j(ij) - k(ij)}
c                           k
c
c
c
c                      1           all  occ
c      u(ij,a) = --------------- { sum  sum   u(kl,a) * ( 4[ij;kl] -
c                clag(j)-clag(i)    k    l
c
c                                     [il;jk] - [ik;jl] )
c
c               + dclag(ij,a) -  clag(j) * ds(ij,a) }
c
c***references         "common sense in quantum chemistry",
c                       b. h. lengsfield iii, encyclopedia galactica,
c                       19072 (2083) 1-2.
c
c***routines called    iosys
c
c***end prologue       fillu
c
      implicit integer (a-z)
c
      character*(*) ops
      logical debug
      real*8 values(nnp,ntriang)
      real*8 isq(num,num)
      real*8 clag(num)
      real*8 ds(nnp,nder),c(num,num)
      real*8 t1(num,num),t2(num,num),dclag(num,num,nder)
      real*8 u(num,num,nder)
      real*8 fac
      real*8 intj(nconum,ntriang)
      real*8 fj(nnp,nder),fk(num,nder,num)
      real*8 uk(num,nder,nco)
c
      parameter (debug=.false.)
c
      common /io/ inp,iout
c
      numnum=num*num
c
c     ----- read in the dependent elements of u and put into full u ---
c
      call iosys('read real "cphf solutions" from rwf',
     $     num*num*nder,u,0,' ')
c..nco
      if(nco.ne.0) then
c
      call iosys('read real "mcscf vector" from rwf',num**2,c,0,' ')
c
c     ----- read in the ao derivative core fock matrix and form
c     -----             mo derivative core fock matrix
c
      call iosys('read real "ao derivative core lagrangian" '//
     $     'from dints',nnp*nder,ds,0,' ')
c
      do 99 der=1,nder
         call trtosq(t1,ds(1,der),num,nnp)
         call ebc(t2,t1,c,num,num,num)
         call ebtc(dclag(1,1,der),c,t2,num,num,num)
99    continue
c
c     ----- read in and transform the derivative overlap integrals -----
c
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
c  make uk
c
      do 820 i=1,nco
       do 821 der=1,nder
        do 823 j=1,num
         uk(j,der,i)=u(j,i,der)
 823    continue
 821   continue
 820  continue
c
c
c     ----- get the core lagrangian ------
c
      call iosys('read real "mcscf orbital energies" from rwf',
     $     num,clag,0,' ')
c
      if(debug) then
         write(iout,*)' mcscf orbital energies '
         write(iout,88001)(clag(i),i=1,num)
88001    format(5(2x,f12.5))
      endif
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
c     ----- add in the lagrangian term:
c
         ij=0
         do 125 i=1,nco
            do 124 j=1,i
               ij=ij+1
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
         noc=nco+nao
         do 425 i=noc+1,num
            ia=i*(i-1)/2
            do 424 j=noc+1,i
               ij=ia+j
               fac=clag(j)-clag(i)
               if (abs(fac).lt.1.0d-06) then
                  if (i.ne.j) write (iout,121) i,j
                  do 422 der=1,nder
                     u(i,j,der)=-0.5d+00*ds(ij,der)
 422              continue
               else
                  do 423 der=1,nder
                   u(i,j,der)=(dclag(i,j,der)-ds(ij,der)*clag(j)+
     $              4.d0*fj(ij,der)-fk(j,der,i)-fk(i,der,j))/fac
 423              continue
                  do 523 der=1,nder
                     u(j,i,der)=-u(i,j,der)-ds(ij,der)
 523              continue
               end if
 424        continue
 425     continue
 426  continue
c
c
      if (debug) then
         do 99001 der=1,3
            write(iout,*)' dclag ndf= ',der
            call matout(dclag(1,1,der),num,num,num,num,iout)
            do 991 i=1,num
               do 992 j=1,num
                  dclag(i,j,der)=fk(i,der,j)
 992           continue
 991        continue

            write(iout,*)' fk ndf= ',der
            call matout(dclag(1,1,der),num,num,num,num,iout)
            write(iout,*)' fj ndf= ',der
            call print(fj(1,der),nnp,num,iout)
            write(iout,*)' ds ndf= ',der
            call print(ds(1,der),nnp,num,iout)
99001    continue
      end if

      if (debug)then
         do 132 der=1,nder
               write (iout,131) der
 131           format(5x,'cphf solution vector, derivative ',i3)
               call matout(u(1,1,der),num,num,num,num,iout)
 132     continue
      end if
c..nco
      end if
c
c
c     ----- and put the solution vectors out -----
c
      call iosys('write real "cphf solutions" to rwf',
     $     num**2*nder,u,0,' ')
c
c
      return
      end

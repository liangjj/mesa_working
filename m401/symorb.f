*deck @(#)symorb.f	5.1  11/6/94
      subroutine symorb(c,s,v,t,u,vs,num,nnp,nsmall,ops)
c***begin prologue     symorb
c***date written       890307  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           symmetrize initial guess for e-scattering
c***author             lengsfield, byron (llnl)
c***source             @(#)symorb.f	5.1   11/6/94
c***purpose            order mos by symmetry
c***description
c
c***references
c***routines called    iosys(io)
c***end prologue       symorb
c
      implicit integer(a-z)
      real*8 vs(num,nsmall),sdot,xnorm,v(num),small
      real*8 c(num,num),u(num,num),t(num,num),s(nnp),csym(64)
      character*(*) ops
      dimension nbfs(8),nbfspt(64)
      data small/1.d-4/
      save small
      common /io/ inp,iout
c
c
c
      call iosys('read real "overlap integrals" from rwf',nnp,s,0,' ')
c
      call trtosq(u,s,num,nnp)
c
      call matout(u,num,num,num,num,iout)
c
      nsym=intkey(ops,'symmetry=nsym',0,' ')
      call intarr(ops,'nbfs',nbfs,8,' ')
      call intarr(ops,'nbfspt',nbfspt,64,' ')
      call fparr(ops,'csym',csym,64,' ')
c
c      write(iout,*)' input symmetry orbitals  nsym = ',nsym
c
      if(nsym.eq.0) call lnkerr(' m401:symorb: nsym=0 ')
c      call matout(c,num,num,num,num,iout)
      iks=0
      ix=num-nsmall+1
      is=0
      do 1 i=1,nsym
c       write(iout,*)' symmetry ',i
       numsym=nbfs(i)
       call rzero(s,num)
       do 2 j=1,numsym
       s(nbfspt(is+j))=csym(is+j)
  2    continue
c       write(iout,*)' s-vector '
c       call matout(s,num,1,num,1,iout)
       is=is+numsym
       call ebc(t,u,s,num,num,1)
       xnorm=sdot(num,t,1,s,1)
       xnorm=1.d0/sqrt(xnorm)
       call sscal(num,xnorm,t,1)
c
c       write(iout,*)' t-vector '
c       call matout(t,num,1,num,1,iout)
c
       call ebtc(v,c(1,ix),t,nsmall,num,1)
c
c       write(iout,*)' v-vector '
c       call matout(v,nsmall,1,nsmall,1,iout)
c
       do 3 j=1,nsmall
        if(abs(v(j)).gt.small) then
         iks=iks+1
         call scopy(num,c(1,ix-1+j),1,vs(1,iks),1)
        end if
   3  continue
c
   1  continue
c
      if(iks.ne.nsmall) then
      write(iout,*)' symmetry orbitals found and total ',iks,nsmall
      call lnkerr(' m401: symmetry error ')
      end if
c
      ntot=nsmall*num
      call scopy(ntot,vs,1,c(1,ix),1)
c
      return
      end

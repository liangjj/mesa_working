*deck @(#)cilag.f	5.1  11/6/94
      subroutine cilag(num,nnp,lag,ints,tpdm,ntriang,t1,t2,ops)
c
c***begin prologue     cilag
c***date written       871201   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)cilag.f	5.1   11/6/94
c
c***purpose            construct the ci lagrangian
c
c***description
c
c               all                     all
c    x(i,j) = 2 sum h(i,k) * q(k,j) + 4 sum i(i,m,k,l) * g(m,j,k,l)
c                k                       m
c
c***references         "generalization of analytic configuration
c   interaction (ci) gradient techniques for potential energy
c   hypersurfaces, including a solution to the coupled perturbed
c   hartree-fock equations for multiconfiguration scf molecular
c   wave functions", yoshihiro osamura, yukio yamaguchi and henry f.
c   schaefer iii, j. chem. phys. 77(1), p.383-390, 1 july 1982.
c
c***routines called    (none)
c
c***end prologue       cilag
c
      implicit integer (a-z)
c
      character*(*) ops
      logical logkey
      real*8 lag(num,num)
      real*8 ints(nnp,ntriang)
      real*8 tpdm(nnp,ntriang)
      real*8 t1(num,num)
      real*8 t2(num,num)
c
      common /io/ inp,iout
c
c     ----- one-electron part of lagrangian -----
c
c
      write(iout,1)
  1   format(1x,'m1031:ci lagrangian construction ')
c
      call rzero(lag,num**2)
c
      call iosys('read real "mo one-electron integrals" from tints',
     $     nnp,ints,0,' ')
      call iosys('read real "mo 1pdm" from moden',nnp,tpdm,0,' ')
      call trtosq(t1,ints,num,nnp)
      call trtosq(t2,tpdm,num,nnp)
c
c       write(iout,*)' '
c       write(iout,*)' 1e integrals '
c       call matout(t1,num,num,num,num,iout)
c       write(iout,*)
c       write(iout,*)' 1e density   '
c       call matout(t2,num,num,num,num,iout)
c
c
      call ebc(lag,t1,t2,num,num,num)
c
c     ----- multiply by factor of 2 -----
c
      call sscal(num**2,2.0d+00,lag,1)
c
      if (logkey(ops,'print=ci=lagrangian',.false.,' ')) then
         write (iout,290)
 290     format (5x,'the one-electron contribution to the ',
     $        'ci lagrangian:')
         call matout(lag,num,num,num,num,iout)
      end if
c
c     ----- rewind the mo integrals and tpdm -----
c
      call iosys('rewind "mo two-electron integrals" on tints',
     $           0,0,0,' ')
      call iosys('rewind "mo 2pdm" on moden',0,0,0,' ')
c
c     ----- read through integrals and tpdm -----
c
      maxkl=0
      kl=0
      do 110 k=1,num
         do 100 l=1,k
            kl=kl+1
c
c           ----- check that this triangle of integrals is in core ----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "mo two-electron integrals" '//
     $              'from tints without rewinding',lnread,ints,0,' ')
               call iosys('read real "mo 2pdm" '//
     $              'from moden without rewinding',lnread,tpdm,0,' ')
            end if
c
c           ----- square up and multiply together -----
c
            call trtosq(t1,ints(1,kl-minkl+1),num,nnp)
            call trtosq(t2,tpdm(1,kl-minkl+1),num,nnp)
c
c           ----- factor of 4 (8 for k!=l for kl and lk) -----
c
            if (k.eq.l) then
               call sscal(num**2,4.0d+00,t1,1)
            else
               call sscal(num**2,8.0d+00,t1,1)
            end if
c
c           ----- the works -----
c
            call apbc(lag,t1,t2,num,num,num)
c
c
  100    continue
  110 continue
c
      call iosys('write real "ci lagrangian" to rwf',num**2,lag,0,' ')
c
      if (logkey(ops,'print=ci=lagrangian',.false.,' ')) then
         write (iout,300)
 300     format (5x,'ci lagrangian:')
         call matout(lag,num,num,num,num,iout)
      end if
c
c
      return
      end

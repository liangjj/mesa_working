*deck @(#)trlag.f	5.1  11/6/94
      subroutine trlag(c,ao,mo,t1,noc,num,opnscf,ops,nact,ncor)
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      921016   (yymmdd)
c
c   16 october 1992    rlm at lanl
c      testing nact so that closed shell quadratic scf will behave.
c   18 december 1987   bhl at brl
c      opnscf option is ignored so the mcscf routes can
c      be tested against the openshell scf routes
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)trlag.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit integer (a-z)
c
      character*(*) ops
      logical logkey
      real*8 c(num,num)
      real*8 ao(num,num)
      real*8 mo(num,noc)
      real*8 t1(num,noc)
c
      common /io/ inp,iout
c
      nnp=num*(num+1)/2
c
      if (opnscf.ne.0) then
         call iosys('read real "scf vector" from rwf',num**2,c,0,' ')
      else
         call iosys('read real "mcscf vector" from rwf',num**2,c,0,' ')
      endif
c
c
      call iosys('read real mcscf_mo_lagrangian from rwf',
     $            num*noc,mo,0,' ')
c
      call sscal(num*noc,0.5d+00,mo,1)
c
      if (logkey(ops,'print=scf=lagrangian=mo',.false.,' ')) then
         write (iout,120)
 120     format('1',//,t20,'mo hartree-fock lagrangian',/)
         call matout(mo,num,noc,num,noc,iout)
      end if
c
      call ebc(t1,c,mo,num,num,noc)
      call ebct(ao,t1,c,num,noc,num)
      call sqtotr(t1,ao,num,nnp)
c
      if (logkey(ops,'print=scf=lagrangian=ao',.false.,' ')) then
         write (iout,130)
 130     format ('1',//,t20,'ao hartree-fock lagrangian',/)
         call print(t1,nnp,num,iout)
      end if
c
      call iosys('write real "scf ao lagrangian" to rwf',
     $           nnp,t1,0,' ')
c
      call iosys('read real mcscf_mo_lagrangian from rwf',
     $            num*noc,mo,0,' ')
c
      call sscal(num*noc,0.5d+00,mo,1)
c
      if (logkey(ops,'print=mcscf=lagrangian=mo',.false.,' ')) then
          write (iout,140)
 140     format('1',//,t20,'mo mcscf lagrangian',/)
         call matout(mo,num,noc,num,noc,iout)
      end if
c
      call ebc(t1,c,mo,num,num,noc)
      call ebct(ao,t1,c,num,noc,num)
      call sqtotr(t1,ao,num,nnp)
c
      if (logkey(ops,'print=mcscf=lagrangian=ao',.false.,' ')) then
         write (iout,150)
 150     format ('1',//,t20,'ao mcscf lagrangian',/)
         call print(t1,nnp,num,iout)
      end if
c
      call iosys('write real "mcscf ao lagrangian" to rwf',
     $            nnp,t1,0,' ')
c
c     ----- active density matrix -----
c
      nnpact=nact*(nact+1)/2
      if(nact.eq.0) then
         call rzero(t1,nnpact)
      else
         call iosys('read real "mo 1pdm" from mcscr',nnpact,t1,0,' ')
c
         if (logkey(ops,'print=mcscf=active-density=mo',.false.,' '))
     $        then
            write (iout,160)
 160        format ('1',//,t20,'ao mcscf active density matrix',/)
            call print(t1,nnpact,nact,iout)
         end if
c
         call trtosq(mo,t1,nact,nnpact)
         call ebc(t1,c(1,ncor+1),mo,num,nact,nact)
         call ebct(ao,t1,c(1,ncor+1),num,nact,num)
         call sqtotr(t1,ao,num,nnp)
c
         if (logkey(ops,'print=mcscf=active-density=ao',.false.,' '))
     $        then
            write (iout,170)
 170        format ('1',//,t20,'ao mcscf active density matrix',/)
            call print(t1,nnp,num,iout)
         end if
c
         call iosys('write real "mcscf ao active density" to rwf',
     $              nnp,t1,0,' ')
c
      endif
c
c
      if(ncor.ne.0) then
         call iosys('read real "mcscf ao core density" from rwf',
     $               nnp,ao,0,' ')
         call vadd(t1,t1,ao,nnp)
      end if
      call iosys('write real "mcscf ao 1pdm" to rwf',nnp,t1,0,' ')
c
c
      return
      end

*deck @(#)prpone.f	5.2 11/28/95
      subroutine prpone(nbf,natoms,prpint,d,scr,nobf,mintyp,nbtype,
     $                  zan,ddel,temp,nx,ny,nz,caltyp,ops,trans,refops)
c***begin prologue     prpone.f
c***date written       860722  
c***revision date      11/28/95      
c   august 17, 1995    rlm at lanl
c      replacing definition of one and four with data statements.
c      see routine phyfil in m1.
c
c***keywords           
c***author             martin, richard(lanl)
c***source             @(#)prpone.f	5.2   11/28/95
c***purpose            
c***description
c   traces the property integrals with the density.
c***references
c
c***routines called
c
c***end prologue       prpone.f
      implicit none
c     --- input variables -----
      integer nbf,natoms,nbtype
      logical trans
      character*(*) caltyp,ops
      character*(*) refops
c     --- input arrays (unmodified) ---
      integer mintyp(nbtype),nobf(nbtype),nx(*),ny(*),nz(*)
      real*8 d(nbf,nbf,*),zan(natoms)
c     --- input arrays (scratch) ---
      real*8 prpint(nbf,nbf),scr(nbf,nbf)
      real*8 temp(6,natoms,*),ddel(natoms)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer lmult,vderiv,nnp,prpmom,minprp,prp,powx,powy,powz
      integer at,deriv,nprps
      logical domult,dov,fermi,mv,logkey
      real*8 origin(3),nucprp,eprp,totprp,conprp,etoesu,tenten,toang
      real*8 sdot,fines,finesq,massv,darwin,cowan,pi,pi3
      real*8 zero,one,two,three,four,eight
      real*8 efg(6),dipolar(6),eigvec(3,3),eigval(3),t1(3),t2(3)
      real*8 trace,dpsum
      character*4 funcnm,prpnam
      character*1 itoc
      character*6 mult
      character*12 mltnam(0:4)
      character*14 vnam(0:2)
c
      parameter (zero=0.0d+00,two=2.0d+00,eight=8.0d+00)
      parameter (three=3.0d+00,tenten=1.0d+10)
      data one/1.0d0/,four/4.0d0/
      data domult/.false./, dov/.false./
      data fermi/.false./, mv/.false./
      data mltnam/'monopole    ','dipole      ','quadrupole  ',
     $            'octupole    ','hexadecapole'/
      data vnam/'potential     ','field         ','field gradient'/
      save domult,dov,fermi,mv,mltnam,vnam
c
      common/io/inp,iout
c
 1000 format(5x,'origin centered properties:')
 1010 format(8x,'component',7x,'electronic',6x,'nuclear',4x,'total(au)',
     $           2x,'total(esu-ang**)')
 1015 format(8x,a12)
 1020 format(11x,a6,4x,3f13.5,4x,f14.5)
 1022 format(8x,a6,7x,3f13.5)
 1025 format(8x,'field gradient eigenvectors (traceless tensor):')
 1026 format(8x,a11,2x,3f13.5)
 1027 format(8x,'spin field gradient eigenvectors (traceless tensor):')
 1030 format(5x,'atom centered properties:')
 1040 format(8x,'component',7x,'electronic',6x,'nuclear',4x,'total(au)',
     $           2x,'total(esu**2/ang**)')
 1042 format(/,8x,'atom',i2)
 1045 format(8x,a14)
c
c     --- evaluate pi.
      pi=four*atan(one)
      pi3=pi**3
c
c     --- check for the multipole integrals.
      lmult=-1
      if(logkey(ops,'properties=e0',.false.,refops)) lmult=0
      if(logkey(ops,'properties=e1',.false.,refops)) lmult=1
      if(logkey(ops,'properties=e2',.false.,refops)) lmult=2
      if(logkey(ops,'properties=e3',.false.,refops)) lmult=3
      if(logkey(ops,'properties=e4',.false.,refops)) lmult=4
      if(lmult.ge.0) domult=.true.
c
c     --- check for electrostatic properties.
      vderiv=-1
      if(logkey(ops,'properties=v0',.false.,refops)) vderiv=0
      if(logkey(ops,'properties=v1',.false.,refops)) vderiv=1
      if(logkey(ops,'properties=v2',.false.,refops)) vderiv=2
      if(vderiv.ge.0) dov=.true.
c
c     --- check for the relativistic integrals.
      fermi=logkey(ops,'properties=fermi',.false.,refops)
      mv=logkey(ops,'properties=mv',.false.,refops)
c
c     --- prepare the total density matrix.
c         for hartree-fock runs the closed shell density should be 
c         doubled.
      if(caltyp.eq.'rhf'.or.caltyp.eq.'rohf') then
         call vadd(d(1,1,1),d(1,1,1),d(1,1,1),nbf*nbf)
         if(caltyp.eq.'rohf') then
            call vadd(d(1,1,1),d(1,1,1),d(1,1,2),nbf*nbf)
         endif
      endif
c
c     --- evaluate the multipole integrals.
      nnp=nbf*(nbf+1)/2
      call iosys('read real esu/e- from rwf',1,etoesu,0,' ')
      etoesu=etoesu*tenten
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
      if(domult) then
         write(iout,1000)
         write(iout,1010)
      endif
      do 20 prpmom=0,lmult
         write(iout,1015) mltnam(prpmom)
         minprp=mintyp(prpmom+1)
         do 10 prp=minprp,minprp+nobf(prpmom+1)-1
            powx=nx(prp)
            powy=ny(prp)
            powz=nz(prp)
            mult='e'//itoc(prpmom)//funcnm(powx,powy,powz)
            call iosys('read real '//mult//' from rwf after rewinding',
     $                  3,origin,0,' ')
            call iosys('read real '//mult//' from rwf without rewinding'
     $                  ,1,nucprp,0,' ')
            call iosys('read real '//mult//' from rwf without rewinding'
     $                  ,nnp,scr,0,' ')
            call trtosq(prpint,scr,nbf,nnp)
c
c           --- if this is a transition operator, 
c               zero the nuclear contribution.
            if(trans) nucprp=zero
c
c           --- the integrals have been computed for the operator x, e.g.,
c               they should be converted to -ex, etc.
            eprp=-sdot(nbf*nbf,d(1,1,1),1,prpint,1)
            totprp=eprp+nucprp
c
c           --- convert from atomic units to esu-ang,esu-ang**2,etc.
            conprp=totprp*etoesu*(toang**prpmom)
            write(iout,1020) mult,eprp,nucprp,totprp,conprp
   10    continue
   20 continue
c
c     --- evaluate the relativistic corrections.
      call iosys('read real fine-structure from rwf',1,fines,0,' ')
      finesq=fines*fines
      if(mv) then
         call iosys('read real "del4 integrals" from rwf',
     $               3,temp(1,1,1),0,' ')
         call iosys('read real "del4 integrals" from rwf '//
     $              'without rewinding',1,nucprp,0,' ')
c
c        if this is a transition operator, zero nuclear contribution.
         if(trans) nucprp=zero
         call iosys('read real "del4 integrals" from rwf '//
     $              'without rewinding',nnp,scr,0,' ')
         call trtosq(prpint,scr,nbf,nnp)
         massv=sdot(nbf*nbf,d(1,1,1),1,prpint,1)
         massv=-massv*finesq/eight
         mult='mass-v'
         write(iout,1022) mult,massv,nucprp,massv+nucprp
      endif
      if(fermi) then
         call rzero(ddel,natoms)
c        read the atomic numbers.
         call iosys('read real "nuclear charges" from rwf',-1,zan,0,' ')
         call iosys('rewind dirac_delta on rwf',0,0,0,' ')
         darwin=zero
         do 23 at=1,natoms
            call iosys('read real dirac_delta from rwf '
     $               //'without rewinding',3,temp(1,1,1),0,' ')
            call iosys('read real dirac_delta from rwf '
     $               //'without rewinding',1,nucprp,0,' ')
c
c           if this is a transition operator, zero nuclear contribution.
            if(trans) nucprp=zero
            call iosys('read real dirac_delta from rwf '
     $               //'without rewinding',nnp,scr,0,' ')
            call trtosq(prpint,scr,nbf,nnp)
c           the darwin term is computed with the total density.
            darwin=darwin+zan(at)*sdot(nbf*nbf,d(1,1,1),1,prpint,1)
c           the fermi contact term with the open-shell density.
            if(caltyp.eq.'rohf') ddel(at)=sdot(nbf*nbf,d(1,1,2),1,
     $           prpint,1)
   23    continue
         darwin=darwin*two*pi3*finesq
         mult='darwin'
         write(iout,1022) mult,darwin,nucprp,darwin+nucprp
      endif
c
c     --- print the cowan-griffin relativistic correction to the energy.
      if(mv.and.fermi) then
         mult='cowan'
         cowan=massv+darwin
         nucprp=zero
         write(iout,1022) mult,cowan,nucprp,cowan+nucprp
      endif
c
c     --- evaluate the electrostatic properties.
      do 40 deriv=0,vderiv
         minprp=mintyp(deriv+1)
         do 30 prp=minprp,minprp+nobf(deriv+1)-1
            powx=nx(prp)
            powy=ny(prp)
            powz=nz(prp)
            mult='v'//itoc(deriv)//funcnm(powx,powy,powz)
            call iosys('rewind '//mult//' on rwf',0,0,0,' ')
            do 25 at=1,natoms
               call iosys('read real '//mult//' from rwf '
     $                  //'without rewinding',3,temp(1,at,prp),0,' ')
               call iosys('read real '//mult//' from rwf '
     $                  //'without rewinding',1,temp(4,at,prp),0,' ')
               call iosys('read real '//mult//' from rwf '
     $                  //'without rewinding',nnp,scr,0,' ')
               call trtosq(prpint,scr,nbf,nnp)
c
c              the integrals have been computed for the operator -z/r, e.g.
               temp(5,at,prp)=sdot(nbf*nbf,d(1,1,1),1,prpint,1)
c              --- if this is an open-shell run, evaluate the
c                  spin density contribution.
               if(caltyp.eq.'rohf') then
                  temp(6,at,prp)=sdot(nbf*nbf,d(1,1,2),1,prpint,1)
               endif
   25       continue
   30    continue
   40 continue
c
c     --- print the atom centered properties.
      if(dov.or.fermi) then
         write(iout,1030)
         write(iout,1040)
      endif
      do 60 at=1,natoms
         if(dov) write(iout,1042) at
         call rzero(efg,6)
         call rzero(dipolar,6)
         trace=zero
         dpsum=zero
         do 50 deriv=0,vderiv
            write(iout,1045) vnam(deriv)
            minprp=mintyp(deriv+1)
            nprps=minprp+nobf(deriv+1)-1
            do 45 prp=minprp,nprps
               powx=nx(prp)
               powy=ny(prp)
               powz=nz(prp)
               prpnam=funcnm(powx,powy,powz)
               mult='v'//itoc(deriv)//funcnm(powx,powy,powz)
               nucprp=temp(4,at,prp)
c
c              if this is a transition operator, zero nuclear contribution.
               if(trans) nucprp=zero
               eprp=temp(5,at,prp)
               totprp=nucprp+eprp
c
c              --- save the efg tensor as a lower triangle.
c                  also accumulate the trace.
               if(deriv.eq.2) then
                  if(prpnam.eq.'xx') then
                     efg(1)=totprp
                     trace=trace+totprp
                  else if(prpnam.eq.'xy') then
                     efg(2)=totprp
                  else if(prpnam.eq.'yy') then
                     efg(3)=totprp
                     trace=trace+totprp
                  else if(prpnam.eq.'xz') then
                     efg(4)=totprp
                  else if(prpnam.eq.'yz') then
                     efg(5)=totprp
                  else if(prpnam.eq.'zz') then
                     efg(6)=totprp
                     trace=trace+totprp
                  endif
                  if(caltyp.eq.'rohf') then
c                    save the spin dipolar hyperfine tensor
                     if(prpnam.eq.'xx') then
                        dipolar(1)=temp(6,at,prp)
                        dpsum=dpsum+temp(6,at,prp)
                     else if(prpnam.eq.'xy') then
                        dipolar(2)=temp(6,at,prp)
                     else if(prpnam.eq.'yy') then
                        dipolar(3)=temp(6,at,prp)
                        dpsum=dpsum+temp(6,at,prp)
                     else if(prpnam.eq.'xz') then
                        dipolar(4)=temp(6,at,prp)
                     else if(prpnam.eq.'yz') then
                        dipolar(5)=temp(6,at,prp)
                     else if(prpnam.eq.'zz') then
                        dipolar(6)=temp(6,at,prp)
                        dpsum=dpsum+temp(6,at,prp)
                     endif
                  endif
               endif
               conprp=totprp*etoesu*etoesu/(toang*(toang**deriv))
               write(iout,1020) mult,eprp,nucprp,totprp,conprp
   45       continue
c
c           --- if we have done the electric field gradient,
c               do it again in the diagonal representation.
c               note that we also diagonalize the traceless tensor,
c               i.e. 3xixj-r**2.  the term in r**2 has already been
c               taken out of the nuclear contribution.)
            if(deriv.eq.2) then
               efg(1)=efg(1)-(trace/three)
               efg(3)=efg(3)-(trace/three)
               efg(6)=efg(6)-(trace/three)
               call degrsp(3,6,efg,eigval,1,eigvec,t1,t2)
               write(iout,1025)
               write(iout,1026) 'eigenvalues',eigval
               write(iout,1026) '          x',eigvec(1,1),
     $                          eigvec(1,2),eigvec(1,3)
               write(iout,1026) '          y',eigvec(2,1),
     $                          eigvec(2,2),eigvec(2,3)
               write(iout,1026) '          z',eigvec(3,1),
     $                          eigvec(3,2),eigvec(3,3)
c
c              --- if this is an open-shell, repeat for the spin-dipolar 
c              hyperfine tensor.
               if(caltyp.eq.'rohf') then
                  dipolar(1)=dipolar(1)-(dpsum/three)
                  dipolar(3)=dipolar(3)-(dpsum/three)
                  dipolar(6)=dipolar(6)-(dpsum/three)
                  call degrsp(3,6,dipolar,eigval,1,eigvec,t1,t2)
                  write(iout,1027)
                  write(iout,1026) 'eigenvalues',eigval
                  write(iout,1026) '          x',eigvec(1,1),
     $                             eigvec(1,2),eigvec(1,3)
                  write(iout,1026) '          y',eigvec(2,1),
     $                             eigvec(2,2),eigvec(2,3)
                  write(iout,1026) '          z',eigvec(3,1),
     $                             eigvec(3,2),eigvec(3,3)
               endif
            endif
   50    continue
c
c        finish atomic properties with fermi term if requested.
         if(fermi) then
            mult='fermi'
            nucprp=zero
            totprp=ddel(at)
            write(iout,1022) mult,ddel(at),nucprp,totprp
         endif
   60 continue
c
c
      return
      end

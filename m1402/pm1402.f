*deck @(#)pm1402.f	5.1  11/6/94
      subroutine pm1402(z,a,ncore)
c
c old_routine taorb
c
      implicit real*8(a-h,o-z)
cc
      real*8 z(*)
      integer a(ncore)
      integer nbf(20),nob(20)
      integer isoff(8),itoff(8),iloff(8)
      integer ilaoff(8),itaoff(8),nco(8),nao(8),nvo(8),
     $        nocc(8)
      integer nacpair(2,20)
      integer wpadti,orbtbf
      real*8 grad(300)
      real*8 ugrad(300)
      real*8 eci(10)
      real*8 origin(3),nucprp(3)
c
      character*4096 ops
      character*6 mult
      character*4 label(3)
      character*4 funcnm,itoc
      logical prtmo,prtds,prtdu,prtcphf
      logical prtden,prtall
      logical logkey
c
      common /io/ inp,iout
c
      data label/'  x ','  y ','  z '/
      save label
c..
cc
      small=1.d-6
      smltst=1.d-4
      mxpndf=1
c
c..bhl
      write(iout,1)
 1    format(1x,'m1402:process half derivative overlap integrals ')
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
      prtmo=logkey(ops,'halfs=print=mo',.false.,' ')
      prtds=logkey(ops,'halfs=print=overlap',.false.,' ')
      prtden=logkey(ops,'halfs=print=density',.false.,' ')
      prtdu=logkey(ops,'halfs=print=du',.false.,' ')
      prtcphf=logkey(ops,'halfs=print=cphf',.false.,' ')
      if(logkey(ops,'halfs=print=all',.false.,' ')) then
         prtall=.true.
         prtmo=.true.
         prtds=.true.
         prtden=.true.
         prtdu=.true.
         prtcphf=.true.
      else
         prtall=.false.
      end if
c
      call iosys('read integer "number of atoms" from rwf',
     #           1,natoms,0,' ')
c
cc
      ndf=natoms*3
cc
      nsym=1
      call iosys('read integer mc_nbasis from rwf',1,nbf(1),0,' ')
      call iosys('read integer mc_ncore from rwf',1,nco(1),0,' ')
      call iosys('read integer mc_nactive from rwf',1,nao(1),0,' ')
      call iosys('read integer mc_norbs from rwf',1,nob(1),0,' ')
      nvo(1)=nob(1)-nco(1)-nao(1)
      lnlag=nob(1)*(nco(1)+nao(1))
c
      nbft = 0
      do 5 i = 1, nsym
         nbft = nbft + nbf(i)
         nob(i) = nco(i) + nao(i) + nvo(i)
    5 continue
      mxcoef=nbft*nbft
      nso = nbft
      ntot=nso*nbft
c
    7 ix=1
      jx=1
      il=1
      ita=1
      ila=1
      noct=0
      maxmrs=0
      naos2=0
      naos4=0
      naot=0
      do 10 i=1,nsym
         isoff(i)=ix
         itoff(i)=jx
         iloff(i)=il
         ilaoff(i)=ila
         itaoff(i)=ita
         nbfi=nbf(i)
         ncoi=nco(i)
         naoi=nao(i)
         naot=naot+naoi
         naoi2=(naoi*(naoi+1))/2
         naoi4=(naoi2*(naoi2+1))/2
         naos2=naos2+naoi2
         naos4=naos4+naoi4
         noci=ncoi+naoi
         noct=noct+noci
         nocc(i)=noci
         il=il+noci*nbfi
         jx=jx+nbfi*nbft
         ix=ix+nbfi*nbfi
         ila=ila+noci*nbfi*ndf
         ita=ita+nbfi*nbfi*ndf
         mrs=nbfi*nbfi
         maxmrs=max(maxmrs,mrs)
  10  continue
      mtots=ix-1
      lagtot=il-1
cc
      mtot=nbft*nbft
      nbf2=(nbft*(nbft+1))/2
      naot2=(naot*(naot+1))/2
      naot4=(naot2*(naot2+1))/2
cc
      if(prtall) then
         write(iout,*) '  natoms nbft mtot ',natoms,nbft,mtot
      end if
cc
      iorb=1
      iden=iorb+nbft
      ismo=iden+mtot
      iscr=ismo+mtot
      isoao=iscr+mtot
      imo=isoao+mtot
      orbtbf=wpadti(imo+mtot)
      ineed=orbtbf+nbft
c
c
      if(ineed.gt.ncore) go to 1000
c
c
      call iosys('read real "guga square 1pdm" from rwf',
     $            mtot,z(imo),0,' ')
      call iosys('read integer orbtbf from rwf',
     $            nbft,a(orbtbf),0,' ')
c
c reorder density from guga to mo order
c
      call fixden(z(iden),z(imo),a(orbtbf),nbft)
c
      if(prtden) then
         write(iout,*)' '
         write(iout,*)' square density matrix '
         call vecout(z(iden),nbft,nbft)
      endif
c
c--------------------------c
c     read the mcscf vectors
c--------------------------c
c
      call iosys('read real "mcscf vector" from rwf',mtot,
     1           z(imo),0,' ')
c
      if(prtmo) then
         write(iout,9020) 
 9020    format(5x,'orbitals ')
         call vecout(z(imo),nbft,nbft)
      end if
c
      if(logkey(ops,'properties=e1',.false.,' ')) then
c
         idip=iadtwp(ineed)
         ineed=wpadti(idip+mtot)
         if(ineed.gt.ncore)go to 1000
         nnp=nbft*(nbft+1)/2
c
         mult='e1'//funcnm(1,0,0)
         call iosys('read real '//mult//' from rwf after rewinding',
     $               3,origin,0,' ')
         call iosys('read real '//mult//' from rwf without rewinding'
     $              ,1,nucprp(1),0,' ')
         call iosys('read real '//mult//' from rwf without rewinding'
     $              ,nnp,z(iscr),0,' ')
         call trtosq(z(idip),z(iscr),nbft,nnp)
         call ebc(z(iscr),z(idip),z(imo),nbft,nbft,nbft)
         call ebtc(z(idip),z(imo),z(iscr),nbft,nbft,nbft)
c
         xmom=sdot(mtot,z(idip),1,z(iden),1)
c
         mult='e1'//funcnm(0,1,0)
         call iosys('read real '//mult//' from rwf after rewinding',
     $              3,origin,0,' ')
         call iosys('read real '//mult//' from rwf without rewinding'
     $              ,1,nucprp(2),0,' ')
         call iosys('read real '//mult//' from rwf without rewinding'
     $              ,nnp,z(iscr),0,' ')
         call trtosq(z(idip),z(iscr),nbft,nnp)
         call ebc(z(iscr),z(idip),z(imo),nbft,nbft,nbft)
         call ebtc(z(idip),z(imo),z(iscr),nbft,nbft,nbft)
c
         ymom=sdot(mtot,z(idip),1,z(iden),1)
c
         mult='e1'//funcnm(0,0,1)
         call iosys('read real '//mult//' from rwf after rewinding',
     $              3,origin,0,' ')
         call iosys('read real '//mult//' from rwf without rewinding'
     $              ,1,nucprp(3),0,' ')
         call iosys('read real '//mult//' from rwf without rewinding'
     $              ,nnp,z(iscr),0,' ')
         call trtosq(z(idip),z(iscr),nbft,nnp)
         call ebc(z(iscr),z(idip),z(imo),nbft,nbft,nbft)
         call ebtc(z(idip),z(imo),z(iscr),nbft,nbft,nbft)
c
         zmom=sdot(mtot,z(idip),1,z(iden),1)
c
         write(iout,33211) xmom,ymom,zmom
33211    format(5x,'dipole transition moments:',/,
     $            '     x         y         z    ',/,3(1x,f8.5))
c
      end if
c
cc----------------------------------------------cc
c     read derivative overlap integrals
cc----------------------------------------------cc
c
c
      nnp=nbft*(nbft+1)/2
      nbf2=nbft*nbft
c
c
c
      numdf=0
      do 1200 ncent=1,natoms
         do 1100 ixyz=1,3
c
            numdf=numdf+1
c
c=====================================================c
c     fetch half-derivative a.o. overlap integrals    c
c=====================================================c
c
            call iosys('read real "ao half_deriv overlap integrals"'//
     #                 ' from rwf',nbf2,z(isoao),(numdf-1)*nbf2,' ')
c
            if(prtds) then
               write(iout,11233) numdf
11233          format(5x,'perturbation ',i4)
               write(iout,*)' halfs in the ao basis '
               call vecout(z(isoao),nbft,nbft)
            end if
c
            call ebc(z(iscr),z(isoao),z(imo),nbft,nbft,nbft)
            call ebtc(z(isoao),z(imo),z(iscr),nbft,nbft,nbft)
c
c
            if(prtds) then
               write(iout,*)' '
               write(iout,*)' halfs in the mo basis '
               call vecout(z(isoao),nbft,nbft)
            endif
c
c  get cphf contribution
c
            call iosys('read real "cphf solutions" from rwf',
     $                 mtot,z(iscr),(numdf-1)*mtot,' ')
c
            if(prtcphf) then
               write(iout,*)' '
               write(iout,*)' cphf solutions '
               call vecout(z(iscr),nbft,nbft)
            end if
c
c  substract (add) cphf and half-derivatve overlap contributions
c  the result, dij,  should be an antisymmetric matrix
c
            call vseamb(z(isoao),z(iscr),mtot,nbft)
c
c           call vaeapb(z(isoao),z(iscr),mtot,nbft)
c
            if(prtdu) then
               write(iout,*)' '
               write(iout,*)' dij   in the mo basis '
               call vecout(z(isoao),nbft,nbft)
            end if
c
c    trace the square density with the orbital derivatives, dij
c
            result=sdot(mtot,z(isoao),1,z(iden),1)
            ugrad(numdf)=result
c
 1100    continue
 1200 continue
c
      write(iout,1202)
 1202 format(5x,'orbital contribution to the nacme')
      write(iout,1201)(ugrad(i),i=1,ndf)
c
      call iosys('write real orbital_nacme to rwf',
     $            ndf,ugrad,0,' ')      
c
 1201  format(3(1x,f12.8))
c
c
      call iosys('read real "cartesian first derivatives" from rwf',
     $            ndf,grad,0,' ')      
c
      write(iout,6906)
      write(iout,1201)(grad(ii),ii=1,ndf)
 6906 format(5x,'transition ci-gradients from rwf ')
c
        call iosys('read integer "nacme pairs" from rwf',
     $              1,numnac,0,' ')
        call iosys('read integer "nacme counter" from rwf',
     $              1,knacme,0,' ')
        call iosys('read integer "nacme states" from rwf',
     $              2*numnac,nacpair,0,' ')
c
      iroot1=nacpair(1,knacme)
      iroot2=nacpair(2,knacme)
      call iosys('read real "ci energy '//itoc(iroot1)//'" from rwf',
     #           1,eci(1),0,' ')
      call iosys('read real "ci energy '//itoc(iroot2)//'" from rwf',
     #           1,eci(2),0,' ')
c
      write(iout,5907) iroot1,iroot2
 5907 format(5x,'ci energies for  roots ',2i5)
      write(iout,6907)(eci(i),i=1,nroot)
 6907 format(5(1x,f13.8))
c
      del=eci(2)-eci(1)      
      if(abs(del).lt.small) then
         write(iout,*)' '
         write(iout,*)' error: ci roots are degenerate'
         call lnkerr(' stop in m1402 ')
      end if
      if(abs(del).lt.smltst) then
         write(iout,8907)
8907     format(5x,'warning: ci roots are nearly degenerate',/,
     $            '  check mcscf, ci and cphf convergence thresholds ')
      end if
c
      fact=1.d0/(eci(2)-eci(1))
      do 6910 i=1,ndf
         grad(i)=grad(i)*fact
 6910 continue
c
      write(iout,7906)
      write(iout,1201)(grad(ii),ii=1,ndf)
c
      call iosys('write real ci_nacme to rwf',
     $           ndf,grad,0,' ')      
c
 7906 format(5x,'ci-nacme ')
      do 6911 i=1,ndf
         grad(i)=grad(i)+ugrad(i)
 6911 continue
      write(iout,7907)
      write(iout,1201)(grad(ii),ii=1,ndf)
 7907 format(5x,'total ci-nacme ')
c
      call iosys(' write real "total ci nacme" to rwf',
     $           ndf,grad,0,' ')      
c
c
      write(iout,1203)
 1203 format(5x,'halfs completed')
c
      return
c
c
 1000 continue
      write(iout,11)ineed,ncore
   11 format(5x,'increase core: ineed,ncore ',2i8)
      call lnkerr(' m1402: halfs ')
c
c
      stop
      end

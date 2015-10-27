*deck %W%  %G%
      subroutine mchsdr(xjbuf,nbufj,nbufk,xj,txk,xk,xkt,den,hess,hmo,
     $     dab,ldab,f,lfm,fab,lfab,grad,lgrad, noba,lnoba, c,isymc,
     $     aii,aij,bii,bij,bijt, nco,nao,nvo,nocc,nob,nbf, nsym, deg,
     $     temp,tv, iscr, nfiv,nfav,nfhm,ndab,
     $     noci, itran, thrsh,
     $     imix,ilen,mixhes,mixinv,imixo, lok,locsym, locp,
     $     iocore,itflag, icas,nmix,lbufh,isymm,linear,nphes,cz,
     $     it46, igeom , iprtg, lhbuf,lhmo, icphf, nblock,
     $     nhd,ldar, fstrt,jstrt,kstrt, jnij,jmij,jnkl,jmkl,
     $     knij,kmij,knkl,kmkl, jkcore,bufax,bufix,lbufso,ilast,
     $     icc,iac,incorh,
     $     hbuf,shmo,lshmo,lbhmo,asort,ncors,mblkd,lblkd)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      900417   (yymmdd)
c
c      17 april, 1990  rlm at lanl
c          passing mblkd,lblkd to remove mclden entry point.
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
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
cc
      implicit real*8(a-h,o-z)
cc
      dimension xj(2),xk(2),xkt(2),lok(2),txk(2),
     $     den(2),hess(2), dab(2),ldab(2),fab(2),lfab(2),grad(2),
     $     f(2), hmo(2), icas(2), lfm(2),lgrad(2),noba(2),lnoba(2),
     $     nco(2),nao(2),nvo(2),nocc(2),nob(2),
     $     deg(2), temp(2), tv(2), isymc(2), c(2)
      dimension aii(2),aij(2),bii(2),bij(2),bijt(2),nbf(2)
      dimension imix(2),ilen(2),mixhes(2),imixo(2),locsym(2),
     $     isymm(2),mixinv(2),locp(2)
      real*8 bufax(lbufso),bufix(lbufso)
cbigmc
      dimension hbuf(*),shmo(*),lshmo(*),lbhmo(*),asort(*)
cbigmc
      dimension cz(2)
      dimension ihd(80)
      integer mblkd(51,2), lblkd(2)
      integer fstrt
      character*3 ians
c
c
      common /io/ inp,iout
c
      common / number / zero,pt5,one,two,four,eight
      common / nhexpk / nhexab
c
c      call jtime(it1)
c      nhexab=-100
c
ccccccccccccccccccccccccccccccccccccccccccc
c
cps      write (iout,9001)
cps 9001 format(/'0************ orbital hessian section ************')
c
c
c-----------------------------------c
c  build occ. orbital pointer array
c-----------------------------------c
c
      ix=0
      do 4 i=1,nsym
         lnoba(i)=ix
         ncoi=nco(i)
         if(ncoi.ne.0) then
            do 1 j=1,ncoi
               ix=ix+1
               noba(ix)=0
 1          continue
         end if
         naoi=nao(i)
         if(naoi.ne.0) then
            do 3 j=1,naoi
               ix=ix+1
               noba(ix)=j
 3          continue
         end if
    4 continue
      noct=ix
c=====
c     compute starting address for y-matrix on tape46
c     y-matrix is need for mccphf calculation
c     ...i2sec is a utility function which computes
c     ...the number of sectors on disk needed to store
c     ...a given number of i*4 words
c     ....the fock matrix is square at this point
c=====
cbl      call bufset(it46)
      lenlag=lgrad(nsym+1)-1
      lnfock=lfab(nsym+1)
cos      lenli=intowp(lenlag)
cps      lnfi=intowp(lnfock)
      lnblck=14*nblock
c
c-------------------------------------------c
c    setup vector coupling coefficients     c
c-------------------------------------------c
c
      call mcvcc(aii,bii,aij,bij,bijt,deg,nsym,linear)
c
c------------------------------------------------c
c  read the header record from the d.a. data set
c------------------------------------------------c
c
      ipflag = 0
c
c------------------------------
c     build the orbital hessian
c------------------------------
c
      call iosys('does mcscf_hessian exist on rwf',0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create real mcscf_hessian on rwf',lhmo,0,0,' ')
      endif
c
      if(iocore.eq.0) then
         call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
         nbb=min(nbufk,lhmo)
         do 90909 i=1,lhmo,nbb
            nbx=min(nbb,lhmo-i+1)
            call rzero(xjbuf,nbx)
            call iosys('write real mcscf_hessian to rwf '//
     $           'without rewinding',nbx,xjbuf,0,' ')
90909    continue
      endif
c
c-1     write(iout,*)'  before mchess lbufh ',lbufh
c
      lhbuf=lbufh
c
      call mchess(xjbuf,nbufj,nbufk,xj,txk,xk,xkt,den,hess,hmo,
     $     dab,ldab,f,lfm,fab,lfab,grad,lgrad, noba,lnoba, c,isymc,
     $     aii,aij,bii,bij,bijt, nco,nao,nvo,nocc,nob,nbf, nsym, deg,
     $     temp,tv,iscr, nfiv,nfav,nfhm,ndab,
     $     noci, itran, thrsh, mixinv,isymm,
     $     imix,ilen,imixo,lok,locsym,iocore,itflag,nmix,ldahab,ldafab,
     $     linear, lhbuf,lhmo,icphf,
     $     fstrt,jstrt,kstrt,jnij,jmij,jnkl,jmkl,knij,kmij,knkl,kmkl,
     $     jkcore ,bufix,bufax,lbufso,ilast,icc,iac,incorh,
     $     hbuf,shmo,lshmo,lbhmo,asort,ncors,mblkd,lblkd)
c
c-1     write(iout,*)'   after mchess lbufh ',lbufh
c-2      write(iout,*)'  after mchess nphes  ',nphes
c
c.io      call srew(nfhm)
c
c-------------------------------------------------------
c    add the fock operator contributions to the gradient
c-------------------------------------------------------
c===========================================c
c     first save core fock operator in xj   c
c===========================================c
      do 9 i=1,lnfock
         xj(i)=f(i)
    9 continue
c
      do 10 i=1,nsym
         if(nco(i).eq.0)go to 10
         call mcfgrd(grad(lgrad(i)),f(lfab(i)),fab(lfab(i)),
     $        nco(i),nob(i),nbf(i),c(isymc(i)),tv,temp,deg(i))
 10   continue
c
cps      write(iout,*) ' lagrangian '
cps      nnnn=nco(1)+nao(1)
cps      call vecout(grad,nob(1),nnnn)
c
c=====
c        ** orbital hessian without lagrangian **
c=====
c     store lagrangian for geometry optimization
c     store fock operator for final orbital rotations
c     store y-matrix for mccphf problem
c=====
c     fetch the y-matrix buffer size
c=====
cc
cc
      do 11 i=1,nsym
         ihd(i)=nco(i)
         ihd(nsym+i)=nao(i)
         ihd(2*nsym+i)=nvo(i)
         ihd(3*nsym+i)=nbf(i)
 11   continue
      lenlag=lgrad(nsym+1)-1
      lnfock=lfab(nsym+1)
      ihd(33)=nmix
      ihd(34)=lbufh
      ihd(35)=nblock
      ihd(36)=icphf
c     ihd(37)=lbcphf
c     ihd(38)=nsec46
      ihd(39)=lenlag
      ihd(40)=lnfock
      ihd(41)=nhd
c     ihd(42)=lda1
      ihd(43)=noci
      ihd(44)=ldar
c
      call iosys('write real mcscf_ao_total_fock to rwf',
     $     lnfock,f,0,' ')
      call iosys('write real mcscf_ao_core_fock to rwf',
     $     lnfock,xj,0,' ')
c
      if(ilast.ne.0) then
c
         call iosys('write real mcscf_ao_lagrangian to rwf',
     $        lenlag,grad,0,' ')
         call iosys('write integer mcscf_mixhes to rwf',
     $        lenlag,mixhes,0,' ')
c
      end if
c
 15   continue
c
c
      do 20 i=1,nsym
         call mcflag(grad(lgrad(i)),temp,c(isymc(i)),nbf(i),nob(i),
     $        nco(i),nao(i))
c
         if(iocore.eq.1) then
            call mch2h (grad(lgrad(i)), nbf(i),nob(i),nco(i),nao(i),
     $           nvo(i),hmo, mixhes(lgrad(i)), icas(i), nmix)
         else
            call sch2h (grad(lgrad(i)), nbf(i),nob(i),nco(i),nao(i),
     $           nvo(i),hmo, mixhes(lgrad(i)), icas(i), nmix,
     $           shmo,lshmo,lbhmo,asort,ncors)
         endif
 20   continue
c
      if(iocore.eq.0) then
         call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.false.)
ct      call mysort('end',asort,asort,0,0,0,0,0,0,0,0,0)
      endif
c
      if(ilast.ne.0) then
         call iosys('write real mcscf_mo_lagrangian to rwf',
     $        lenlag,grad,0,' ')
      end if
c
c.2
      if (nphes .eq. 0) go to 99
c
      if(iocore.eq.1) then
         write (iout,24)
 24      format('0 ***** sqr. orbital hessian **test***')
         call printm(hmo,nmix,2)
      else
         write (iout,2424)
 2424    format('0 ** new sqr. print *** orbital hessian **test***')
         call prthss(xjbuf,nbufk,xk,nmix)
      endif
 99   continue
c.2
c
      nmixtt=nmix*nmix
c
      if(iocore.eq.1) then
         call iosys('write real mcscf_hessian to rwf',
     $        nmixtt,hmo,0,' ')
      end if
c
      lm=1
      lmixt=0
      do 500 i=1,nsym
         call mcgrd(hmo(lm),grad(lgrad(i)),nbf(i),nco(i),nao(i),nvo(i),
     $        lmix,lmixt,mixhes(lgrad(i)),nob(i))
         lm=lm+lmix
 500  continue
c
c     print the gradient?
      if (iprtg .ne. 0) then
         write(iout,501) lmixt
         write(iout,502)(hmo(mm),mm=1,lmixt)
 501     format(/,'  gradient    lmixt = ',i6)
 502     format(5(2x,f18.10))
      end if
c
c
      call iosys('does mcscf_gradient exist on rwf',0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create real mcscf_gradient on rwf',lmixt,0,0,' ')
      endif
      call iosys('rewind mcscf_gradient on rwf',0,0,0,' ')
      call iosys('write real mcscf_gradient on rwf',lmixt,hmo,0,' ')
c
c
      return
      end

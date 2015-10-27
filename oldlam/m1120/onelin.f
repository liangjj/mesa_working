      subroutine onelin (ops)
      implicit integer(a-z)
      common z(1)
      common /io/ inp, iout
      real *8 z, eryd
      real *8 fpkey
      real *8 rgs, wtt,enbnd, convg, ovtol, charge, chrg
      real *8 rfs, del1, r11, wx, radii, rk, rmat
      real *8 rstr, rfnr, gam1, gam2, scale, r1, r2, r3, r4
      character *(*) ops
      character *3 itoc, ans, guess
      character *80 title, cpass
      character *8 chrkey, mshtyp, vtyp
      character *1600 card
      character *10 type
      logical logkey, varit, prnt, nowgt
      dimension a(1)
      dimension rgs(1000,2), wtt(1000,2), radii(2), npnts(2)
      dimension npnr(20), rstr(20), rfnr(20)
      dimension eryd(100)
      equivalence (a,z)
*
c                  solve one channel scattering
c                  problem by linear algebraic method
c                  two regions are used
*
*
*
*          all lengths are in bohr radii
*
*
*
*
      write (iout,270)
      write (iout,250)
      read (inp,180) title
      write(iout,290) title
c
      call cardin (card)
      charge=fpkey(card,'charge',0.d+00,' ')
      radii(1)=fpkey(card,'left-r-matrix-boundary',10.,' ')
      radii(2)=fpkey(card,'right-r-matrix-boundary',100.,' ')
      angmom=intkey(card,'angular-momentum',0,' ')
      nreg=intkey(card,'no-regions',1,' ')
      call posinp ('$energy',cpass)
      call cardin (card)
      nen=intkey(card,'no-energies',1,' ')
      call fparr (card,'energies',eryd,nen,' ')
c
c     ----- iter=max. no. iterations , convg=convergence criterion, ----
c     ----- ovtol=tolerence for linear dependence of iterates -----
c     ----- maxvec=maximum no. vectors -----
      iter=intkey(ops,'lam=no-iterations',50,' ')
      maxvec=intkey(ops,'lam=max-no-vecs',50,' ')
      convg=fpkey(ops,'lam=conv-criterion',1.d-09,' ')
      ovtol=fpkey(ops,'lam=overlap-tol',1.d-08,' ')
      prnt=logkey(ops,'lam=print',.false.,' ')
      varit=logkey(ops,'variation-iteration',.false.,' ')
      guess=chrkey(ops,'lam=guess','rhs',' ')
      nowgt=.false.
c
*
*
      write (iout,260)
      write (iout,270)
*
*
*
      write (iout,240) iter,maxvec,convg,ovtol
      call rzero(rgs,2000)
      call rzero(wtt,2000)
      rgs(1,1)=1.e-10
*
      do 5000 ireg=1,nreg
      call posinp ('$msh-'//itoc(ireg),cpass)
      call cardin (card)
      nrg=intkey(card,'no-regions',1,' ')
      mshtyp=chrkey(card,'type-quad','legendre',' ')
      call intarr (card,'no-points/reg',npnr,nrg,' ')
      call fparr (card,'start-r/reg',rstr,nrg,' ')
      call fparr (card,'final-r/reg',rfnr,nrg,' ')
c
c     ----- nrg=no. regions of integration , nps=no. points in a region
c     ----- rstart=starting value for a region , rfin=final value -----
c     ----- mshtyp=0 is gauss-legendre quadrature =1 trapezoidal rule --
c     ----- rgs and wtt hold the point and integration weights -----
c
      write (iout,300) nrg

c     first and last point special
      ic=1
      do 40 iz=1,nrg
      write (iout,310) npnr(iz),rstr(iz),rfnr(iz)
      rfs=rfnr(iz)-rstr(iz)
      nps=npnr(iz)
      if (mshtyp.eq.'legendre') then
      do 20 jz=1,nps
      ic=ic+1
      call lgndrx (nps,jz,wtt(ic,ireg),rgs(ic,ireg))
      wtt(ic,ireg)=wtt(ic,ireg)*rfs
      rgs(ic,ireg)=rgs(ic,ireg)*rfs+rstr(iz)
   20 continue
      else
*
      del1=rfs/nps
      r11=rstr(iz)-del1
      nps1=nps+1
      do 30 jv=1,nps1
      ic=ic+1
      r11=r11+del1
      wx=1.
      if (jv.eq.1.or.jv.eq.nps1) wx=.5
      wtt(ic,ireg)=wx*del1
   30 rgs(ic,ireg)=r11
      endif
*
   40 continue
*
      nptmx=ic+1
      npnts(ireg)=nptmx
      rgs(nptmx,ireg)=radii(ireg)
      wtt(nptmx,ireg)=0.d+00
 5000 continue
      rgs(1,2)=radii(1)
      wtt(1,2)=0.e+00
*
*
      npnt=max(npnts(1),npnts(2))
      chrg=0.
      v1=1
      v2=v1+npnt
      sf1=v2+npnt
      sf2=sf1+npnt
      sb1=sf1
      sb2=sf2
      rhsf=sb2+npnt
      rhsb1=rhsf
      rhsb2=rhsb1+npnt
      if (varit) then
          exvc=rhsb2+npnt
          exiter=exvc+iter*npnt
          aold=exiter+iter*npnt
          anew=aold+iter*iter
          bold=anew+iter*iter
          bnew=bold+iter
          ipvt=bnew+iter
          scr=ipvt+iter
          words=scr+npnt
      else
          mat=rhsb2+npnt
          ipvt=mat+npnt*npnt
          words=ipvt+npnt
      endif
      call getscm(words,z(1),ngot,'m1120',0)
      vloc=v1
      do 6000 ireg=1,nreg
         call modpot(z(vloc),rgs(1,ireg),npnts(ireg),charge,vtyp)
         vloc=vloc+npnt
 6000 continue
c
      do 7000 ien=1,nen
         rk=sqrt(eryd(ien))
c             calculate solution from zero to radii(1)
         type='function'
         call grnset(z(sf1),gam1,eryd(ien),rk,chrg,rgs(1,1),radii(1),
     1               lval,npnts(1),type)
         type='derivative'
         call grnset(z(sf2),gam1,eryd(ien),rk,chrg,rgs(1,1),radii(1),
     1               lval,npnts(1),type)
         call vscale(z(sf2),-1.e+00,npnts(1))
         call rhscal(z(sf1),z(sf2+npnts(1)-1),z(rhsf),npnts(1))
         call vscale(z(rhsf),-1.e+00,npnts(1))
         vloc=v1
         if (varit) then
             if (guess.eq.'coulomb') then
                 type='function'
                 call grnset(z(scr),gam1,eryd(ien),rk,charge,rgs(1,1),
     1                       radii(1),lval,npnts(1),type)
             elseif(guess.eq.'old') then
                 if (ien.eq.1) then
                     type='function'
                     call grnset(z(scr),gam1,eryd(ien),rk,chrg,rgs(1,1),
     1                           radii(1),lval,npnts(1),type)
                 else
                     call copy(z(rhsf),z(scr))
                 endif
             endif
             call slvit (z(sf1),z(sf2),z(vloc),z(rhsf),z(scr),z(exvc),
     1                   z(exiter),wtt(1,1),z(aold),z(anew),
     2                   z(bold),z(bnew),a(ipvt),z(scr),npnts(1),
     3                   iter,convg,ovtol,ops,nowgt,guess)
         else
             call conmat(z(sf1),z(sf2),z(vloc),z(mat),wtt(1,1),npnts(1))
             call sgefa(z(mat),npnts(1),npnts(1),a(ipvt),info)
             call sgesl(z(mat),npnts(1),npnts(1),a(ipvt),z(rhsf),0)
         endif
         title='forward solution'
      if (prnt) then
         call prntrm(title,z(rhsf),npnts(1),1,npnts(1),1,iout)
      endif
         rmat=z(rhsf+npnts(1)-1)
         call phase(rmat,rk,radii(1),charge,vtyp)
         if (nreg.ne.1) then
             vloc=vloc+npnt
             type='derivative'
             call grnset(z(sb1),gam1,eryd(ien),rk,chrg,rgs(1,2),
     1                   radii(1),lval,npnts(2),type)
             type='derivative'
             call grnset(z(sb2),gam2,eryd(ien),rk,chrg,rgs(1,2),
     1                   radii(2),lval,npnts(2),type)
             scale=1.e+00/(gam2-gam1)
             call vscale(z(sb2),scale,npnts(2))
             call rhscal(z(sb2),z(sb1),z(rhsb1),npnts(2))
             call rhscal(z(sb1),z(sb2+npnts(2)-1),z(rhsb2),npnts(2))
             call vscale(z(rhsb2),-1.e+00,npnts(2))
             if (varit) then
                 if (guess.eq.'coulomb') then
                     type='function'
                     call grnset(z(scr),gam1,eryd(ien),rk,charge,
     1                           rgs(1,2),radii(1),lval,npnts(2),type)
                 endif
                 call slvit (z(sb1),z(sb2),z(vloc),z(rhsb1),z(scr),
     1                       z(exvc),z(exiter),wtt(1,2),z(aold),z(anew),
     2                       z(bold),z(bnew),a(ipvt),z(scr),npnts(2),
     3                       iter,convg,ovtol,ops,nowgt,guess)
                 if (guess.eq.'coulomb') then
                     type='function'
                     call grnset(z(scr),gam1,eryd(ien),rk,charge,
     1                           rgs(1,2),radii(2),lval,npnts(2),type)
                 endif
                 call slvit (z(sb1),z(sb2),z(vloc),z(rhsb2),z(scr),
     1                       z(exvc),z(exiter),wtt(1,2),z(aold),z(anew),
     2                       z(bold),z(bnew),a(ipvt),z(scr),npnts(2),
     3                       iter,convg,ovtol,ops,nowgt,guess)
             else
                 call conmat(z(sb1),z(sb2),z(vloc),z(mat),wtt(1,2),
     1                       npnts(2))
                 call sgefa(z(mat),npnts(2),npnts(2),a(ipvt),info)
                 call sgesl(z(mat),npnts(2),npnts(2),a(ipvt),z(rhsb1),0)
                 call sgesl(z(mat),npnts(2),npnts(2),a(ipvt),z(rhsb2),0)
             endif
             if (prnt) then
                 title='backward solution left'
                 call prntrm(title,z(rhsb1),npnts(2),1,npnts(2),1,
     1                       iout)
                 title='backward solution right'
                 call prntrm(title,z(rhsb2),npnts(2),1,npnts(2),1,
     1                       iout)
             endif
            call prrmat(r1,r2,r3,r4,z(rhsb1),z(rhsb2),npnts(2))
            rmat=r1+r2*r4/(rmat-r3)
            call phase(rmat,rk,radii(2),charge,vtyp)
         endif
 7000 continue
      return
c
  180 format (a80)
  240 format (//,20x,'vector generation information',//,5x,'no. of itera
     1tions',13x,i4,//,5x,'max. no. vectors',14x,i4,//,5x,'convergence',
     2 20x,d15.8,//,5x,'overlap tolerence',14x,d15.8)
  250 format (///,20x,'***** linear algebraic scattering program *****')
  260 format (///,20x,'***** initial input data *****')
  270 format (//,1x,'***************************************************
     1***************************************************')
  290 format (//,20x,'***** title of run *****',//,20x,a80)
  300 format (///,5x,'gaussian mesh',//,5x,'regions=',1x,i3,///,26x,
     1'npts',10x,'r st',10x,'r fin')
  310 format (/,27x,i3,6x,f10.5,4x,f10.5)
  320 format (//,5x,'total no. points                          ',2x,i5)
  330 format (//,5x,'reformat local potential file')
      end

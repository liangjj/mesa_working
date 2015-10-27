      subroutine driver (ops)
      implicit integer(a-z)
      common /io/ inp, iout
      dimension eng(100), itp(4)
      logical logkey, statdm, nowgt
      character *8 orbctl, optctl
      character *8 vloc, orbs, cntpsi, optfle
      character *8 itp
      character *8 molsym
      character *(*) ops
      character *8 chrkey
      character *3 itoc, ans, soldir
      character *80 title
      character *800 card
      real *8 convg, ovtol, charge, rmtrad, eng, fpkey
      data itp/ 'vloc','moproj','optfle','cntpsi'/
*
*     ----- basic routine to recover the data to set up the  -----
*     -----           scattering calculation           -----
*
*
*     ----- units are bohr radii and rydbergs -----
*
*
*
      write (iout,90)
      read (inp,10) title
      call iosys ('does "opt pot file" exist on rwf',0,0,0,ans)
      if (ans.ne.'no') then
      call iosys ('read character "opt pot file" from rwf',0,0,0,itp(3))
      endif
      call iosys ('does "static pot file" exist on rwf',0,0,0,ans)
      if (ans.ne.'no') then
      call iosys ('read character "static pot file" from rwf',0,0,0,
     1   itp(1))
      endif
      call iosys ('does "orb proj file" exist on rwf',0,0,0,ans)
      if (ans.ne.'no') then
      call iosys ('read "orb proj file" from rwf',0,0,0,itp(2))
      endif
      call cardin (card)
      vloc=chrkey(card,'sc-mat',itp(1),' ')
      orbs=chrkey(card,'sc-mo',itp(2),' ')
      optfle=chrkey(card,'opt-mat',itp(3),' ')
      cntpsi=chrkey(card,'cnt-fun',itp(4),' ')
      memmax=intkey(card,'max-mem',0,' ')
      call iosys ('write character "soln file" to rwf',0,0,0,cntpsi)
      soldir=chrkey(card,'solver-dir','all',' ')
      write (iout,100) title
      call iosys ('read integer "no. energies" from rwf',1,neng,0,0)
      call iosys ('read real energies from rwf',neng,eng,0,0)
*
      call iosys ('read integer "no. mesh points" from rwf',1,nptmx,0,0)
      call iosys ('read integer "no. trgt states" from rwf',1,nsts,0,0)
      call iosys ('read integer "total no. chan" from rwf',1,ntchn,0,0)
      call iosys ('read integer "max. chn/st" from rwf',1,ncmax,0,0)
      call iosys ('read integer "max. l/chn" from rwf',1,lplsmx,0,0)
      call iosys ('write real "scat. eng." to rwf',neng,eng,0,0)
      call iosys ('read character "inversion sym" from rwf',0,0,0,molsym
     1 )
      call iosys ('read real "nuclear charge" from rwf',1,charge,0,0)
      call iosys ('read real "r matrix radius" from rwf',1,rmtrad,0,0)
      call iosys ('read integer "total m" from rwf',1,msmtot,0,0)
      call iosys ('read integer "total spin" from rwf',1,spntot,0,0)
      call iosys ('read integer "total parity" from rwf',1,partot,0,0)
      statdm=logkey(ops,'data',.false.,' ')
      nowgt=logkey(ops,'no-weights',.false.,' ')
      call iosys ('read character orbctl from rwf',0,0,0,orbctl)
      call iosys ('read character optctl from rwf',0,0,0,optctl)
      if (statdm) then
      orbctl=chrkey(ops,'orbctl','orb',' ')
      optctl=chrkey(ops,'optctl','opt',' ')
      endif
c     ----- write out some basic information -----
c
      write (iout,20) nsts,ntchn,nptmx,neng
      write (iout,30) molsym,charge,rmtrad
      write (iout,40) msmtot,spntot,partot
      if (orbctl.eq.'orb') then
      call iosys ('read integer "max. l bound" from rwf',1,maxlex,0,0)
      call iosys ('read integer "total l bound" from rwf',1,totlex,0,0)
      call iosys ('read integer "exch lam max" from rwf',1,lammax,0,0)
      call iosys ('read integer "total opt orbs" from rwf',1,nfopt,0,0)
      call iosys ('read integer "total lag orbs" from rwf',1,nflag,0,0)
      write (iout,50) nfopt,nflag,maxlex,lammax
      call iosys ('open mofile as old',0,0,0,orbs)
      call iosys ('read integer "total mos" from mofile',1,nmotot,0,0)
      call iosys ('read integer "bnd upper l" from mofile',1,uperl,0,0)
      write (iout,60) nmotot,uperl
      endif
      call iosys ('read integer "no. iterations" from rwf',1,iter,0,0)
      iter=intkey(ops,'iter',iter,' ')
      call iosys ('read real "convg. criterion" from rwf',1,convg,0,0)
      call iosys ('read real "overlap tol" from rwf',1,ovtol,0,0)
      convg=fpkey(ops,'convg',convg,' ')
      ovtol=fpkey(ops,'ovtol',ovtol,' ')
      call iosys ('read integer "max. no. vectors" from rwf',1,maxvec,0,
     1 0)
c
c     ----- maxvec capability not yet implemented -----
      maxvec=iter
      write (iout,70) iter,convg,ovtol,maxvec
*
      call iosys ('read integer "no. la solns" from rwf',1,nsol,0,0)
      if (optctl.ne.'opt') nsol=ntchn
      call iosys ('read integer "no. solutions" from rwf',1,numdo,0,0)
      call iosys ('read integer "size la matrix" from rwf',1,msize,0,0)
      write (iout,80) nsol,numdo,msize
*     ----- outlay memory and call routines to perform -----
*     -----              calculation                   -----
      call expand (ops,nsts,ntchn,nptmx,neng,nmotot,nfopt,nflag,iter
     1 ,maxvec,maxlex,totlex,lammax,ncmax,lplsmx,uperl,orbctl,optctl
     2 ,charge,rmtrad,convg,ovtol,nsol,numdo,msize,memmax,vloc,orbs,
     3 cntpsi,optfle,soldir,nowgt)
      return
c
   10 format (a80)
   20 format (//,5x,'no. states        ',2x,i4,/,5x,'total no. channels'
     1 ,2x,i4,/,5x,'no. mesh points   ',2x,i4,/,5x,'no. of energies   ',
     2 2x,i4)
   30 format (//,5x,'molecular symmetry:',2x,a8,/,5x,'residual charge',2
     1 x,f12.6,/,5x,'r-matrix radius',2x,f12.6)
   40 format (//,20x,'symmetry information of scattering state',/,5x,'ms
     1ym  ',2x,i3,/,5x,'spin  ',2x,i3,/,5x,'parity',2x,i3)
   50 format (//,5x,'total no. optical potential orbitals',2x,i4,/,5x,'t
     1otal no. laguerre orbitals         ',2x,i4,/,5x,'l max. optical
     2                   ',2x,i4,/,5x,'lambda max. exchange
     3   ',2x,i4)
   60 format (//,5x,'no. m.o. on mofile',2x,i4,/,5x,'largest l value   '
     1 ,2x,i4)
   70 format (//,5x,'max. no. iterations  ',2x,i4,/,5x,'convergence crit
     1erion',2x,d15.8,/,5x,'overlap tolerence    ',2x,d15.8,/,5x,'max. n
     2o. vectors     ',2x,i4)
   80 format (//,5x,'no. linearly independent solutions',2x,i5,/,5x,'no.
     1 solved',2x,i5,/,5x,'size of linear algebraic matrix   ',2x,i5)
   90 format (///,20x,'***** linear algebraic scattering code *****')
  100 format (//,20x,'title of run',//10x,10a8,//)
      end

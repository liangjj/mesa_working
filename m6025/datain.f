*deck datain
      subroutine datain      
      implicit real *8 (a-h,o-z)
c
      character*4096 ops
      character*1600 card
      character*80 cpass
      character *8 ityp(2), name(3), ichan(3)
      character *10 rsn
      character*80 chrkey
      integer chan,type,initial,final
      logical fit, debug, dollar, logkey
      common /doyle/ npad, nit, apad(4), npad1(7), apad1(8), npad2(2)
      common /kar/ apad2(150), bl(30), bu(30), v(30), apad3(60), 
     1             nuflag(30), nbd(30)
      common / kchg  / apad4(300,15), bsav(300), nersdm
      common /kvar/ npar, nprob, ners, npad3(2), apad5, solves,
     1               npad4(2), apad6(2)
      common / pdata / np, npbck, elist(256), plist(256), sreal(256),
     >                 simag(256), phase(256), pabs(256), pres(256),
     >                 pbck(256), edup(256), ilist(256), del(256)
      common /io/ inp, iout
      namelist/resfit/ emin,emax,chan,type,initial,final,fit,debug,
     #                 nit,npbck,escal,eref,rsn,elist,del,name,
     #                 np,xnpi
      emin = -999d0
      emax = -999d0
      chan = 1
      type = 2
      initial = 0
      final = -1
      fit = .true.
      debug = .false.
      de = 0d0
      ne = -1
      escal = 1d0
      eref = 0d0
      nbar = 0
      npar = 3
      npbck = 2
      xnpi = 2.d0
      np = 0
      name(1) = ' '
      name(2) = ' '
      name(3) = ' '
      nit = 10
      rsn = ' '
      read(inp,resfit,end=100)
      do 3 i=1,npar
         read(inp,*) nuflag(i), v(i), nbd(i), bl(i), bu(i)
 3    continue   
c      IF ( dollar('$resfit',card,cpass,inp) ) then 
c           emin=fpkey(card,'emin',emin,' ')
c           emax=fpkey(card,'emin',emax,' ')
c           chan=intkey(card,'no.-chan',chan,' ')
c           type=intkey(card,,'type-scattering',type,' ')
c           initial=intkey(card,'initial-channel',initial,' ')
c           final=intkey(card,'final-channel',final,' ')
c           fit=logkey(card,'fit',fit,' ')
c           debug=logkey(card,'debug',debug,' ')
c           nit=intkey(card,'no.-its',nit,' ')
c           np=intkey(card,'no.-energies',np,' ')
c           npbck=intkey(card,'no.-background-terms',npbck,' ')
c           escal=fpkey(card,'energy-scale-factor',escal,' ')
c           eref=fpkey(card,'reference-energy',eref,' ')
c           rsn=chrkey(card,'label','resonance',' ')
c           call fparr(card,'energies',elist,np,' ')
c           call fparr(card,'phase-shifts',del,np,' ')
c           xnpi=fpkey(card,'xnpi',xnpi,' ')
c           call intarr(card,'vary-flag',nuflag,npar,' ')
c           call intarr(card,'bound-flag',nbd,npar,' ')
c           call fparr(card,'initial-guesses',v,npar,' ')
c           call fparr(card,'lower-bounds',bl,npar,' ')
c           call fparr(card,'upper-bounds',bu,npar,' ')
c      endif
 100  continue
      return
      end








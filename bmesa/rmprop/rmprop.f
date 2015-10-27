*deck rmprop.f 
c***begin prologue     m6236
c***date written       941019   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6236, link 6236, propagation
c***author             schneider, b. i.(nsf)
c***source             m6236
c***purpose            coupled channel r-matrix propagation
c***description        propagation of a set of coupled second order
c***                   differential equations is accomplished using 
c***                   the r-matrix propagation approach of light-walker.
c***references         papers by light-walker j. chem. phys.
c
c***routines called    iosys, util and mdutil
c***end prologue       m6236
      program rmprop
      implicit integer(a-z)
      common z(1)
      real*8 z
      dimension a(1)
      equivalence (a,z)
      real*8 charge, rbeg, rfinal, stpmin, stpmax, fparr
      real*8 beta, fpkey, tmp, psesum, range
      logical prntv, prntr, debug, logkey, iflgl, bypass
      logical posinp, prntmt, prntcp, prntk, prntep, iefl
      character*4096 ops
      character*32 cpass
      character*1600 card
      character*32 chrkey, pottyp
      character*128 fillam
      common /io/ inp, iout
      common/memory/ioff
c
c
c*****************
c
c this program takes the r-matrix calculated at a left hand boundary and
c propagates it to a right hand boundary.  scattering parameters are 
c determined at the two points.
c
c*****************
c
c
c rbeg   = leftmost starting point
c rfinal = asymptotic matching radius (bohr)
c stpmin = min. step size
c stpmax = max. allowed step size
c beta   = step-size condition parameter
c
c
c prntv = .true. - potential every npv points
c prntr = .true. - r matrix every npv points
c debug = .true. - debug print
c nc    = no. of channels
c nen   = no. of energies
c charge = residual ion charge
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      bypass=logkey(ops,'data=input',.false.,' ')
      lmax=intkey(ops,'maximum-l-value',0,' ')
      charge=fpkey(ops,'residual-charge',0.d+00,' ')
      if (charge.ne.0.d0) then
          call lnkerr ('quit. residual charge must be zero')
      endif 
      pottyp=chrkey(ops,'potential-type','square-well',' ')
      prntr=logkey(ops,'print=m6236=r-matrix',.false.,' ')
      prntk=logkey(ops,'print=m6236=k-matrix',.false.,' ')
      prntv=logkey(ops,'print=m6236=potential',.false.,' ')
      prntmt=logkey(ops,'print=m6236=scattering-matrices',.false.,' ')
      prntcp=logkey(ops,'print=m6236=coupling-constants',.false.,' ')
      prntep=logkey(ops,'print=m6236=eigenphase-vectors',.false.,' ')
      debug=logkey(ops,'print=m6236=debug',.false.,' ')
      npv=intkey(ops,'print=m6236=potential=points',1,' ')
      npr=intkey(ops,'print=m6236=r-matrix=points',1,' ')
      npr=intkey(ops,'m6236=points-rmat',1,' ')
      if (posinp('$integration-parameters',cpass) ) then
          call cardin(card)
      endif                    
      stpmax=fpkey(card,'maximum-integration-step',10.d0,' ')
      stpmin=fpkey(card,'minimum-integration-step',1.d-04,' ')
      beta=fpkey(card,'step-size-condition',1.d-05,' ')
      iflgl=logkey(card,'zero-off-diagonal',.false.,' ')
      iefl=logkey (card,'prop=diagonalize-each-energy',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1             -1,0,0,fillam)
      if (bypass) then
          if (posinp('$grid',cpass) ) then
              call cardin(card)
          endif              
          rbeg=fpkey(card,'r-matrix-radius',10.d0,' ')
          rfinal=fpkey(card,'final-r',100.d+00,' ')
          range=fpkey(card,'potential-range',rbeg,' ')   
          if ( posinp('$channels',cpass) ) then
               call cardin(card)
          endif
          nc=intkey(ops,'number-of-channels',1,' ')
          call iosys ('open lamdat as new',0,0,0,fillam)
      else           
          call iosys ('open lamdat as old',0,0,0,fillam)
          call iosys ('read integer "number of channels" from lamdat',
     1                 1,nc,0,' ')
          call iosys('read real "r-matrix radius" from lamdat',1,
     1                rbeg,0,' ')
          call iosys('read real "final radius" from lamdat',1,
     1                rfinal,0,' ')
          call iosys ('read character "potential type" from lamdat',
     1                 -1,0,0,pottyp)
          call iosys('read real "potential range" from lamdat',1,
     1                range,0,' ')
          call iosys ('read integer "maximum l value" from lamdat',1,
     1                 lmax,0,' ') 
          call iosys ('read integer "size of hamiltonian matrix" '//
     1                'from lamdat',1,neig,0,' ')     
      endif
      if ( posinp('$energy',cpass) ) then
           call cardin(card)
      endif           
      nen=intkey(card,'number-of-energies',1,' ')
      write (iout,*)
      write (iout,*)           
      write(iout,*) '                         m6236:coupled '//
     1                         'channel r-matrix propagation'
      write (iout,*)'                               link'
      tmp=lmax*30.d0
      ltop=lmax+sqrt(tmp)
      ltop=max(ltop,lmax)
      lp=ltop+1     
      call iosys ('open scratch as scratch on ssd',1000000,0,0,
     1            'scrtch')
      call iosys ('create real tev on scratch',nc*nc,0,0,0)
      call iosys ('create real vdiag on scratch',-1,0,0,0)
      iflg=1
      if (iflgl) then
          iflg=0
      endif
      ief=0
      if (iefl) then
          ief=1
      endif           
      write (iout,1) nc, rbeg, rfinal, nen, charge
*
*     ----- lay out memory -----
*
      ioff=1
      do 10 i=1,2
         vcij=ioff
         eig=vcij+nc*nc
         srfprj=eig+neig
         ec=srfprj+neig*nc
         lval=wpadti(ec+nc)
         energy=iadtwp(lval+nc)
         j=energy+nen
         jp=j+lp
         y=jp+lp
         yp=y+lp
         wron=yp+lp
         y0=wron+lp
         dy0=y0+nc
         y1=dy0+nc
         dy1=y1+nc
         rmat=dy1+nc
         coef=rmat+nc*nc
         kmat=coef+2*nc*nc
         kmatc=coef
         tmtrx=kmat+nc*nc
         cross=tmtrx+2*nc*nc
         scr=cross+nc*nc
         vmat=coef
         tev=coef+nc*nc
         tp=tmtrx
         tsave=tmtrx+nc*nc
         pdiag=scr+2
         enrg=pdiag+4*nc
         scr1=enrg+nc
         ipvt=wpadti(scr1++3*nc*nc)
         words=ipvt+nc
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'rmprop',0)
         endif
 10   continue
      if (bypass) then
          if ( posinp('$channels',cpass) ) then
               call cardin(card)
          endif
          call intarr(card,'channel-l-values',a(lval),nc,' ')
          call fparr(card,'channel-energies',z(ec),nc,' ')     
      else
          call iosys ('read integer "channel l values" from lamdat',nc,
     1                 a(lval),0,' ')              
          call iosys ('read real "channel energies" from lamdat',nc,
     1                 z(ec),0,' ')      
          call iosys ('read real "hamiltonian eigenvalues" from lamdat',
     1                 neig,z(eig),0,' ')
          call iosys ('read real "hamiltonian eigenfunction surface '//
     1                'projections" from lamdat',neig*nc,z(srfprj),
     2                 0,' ')               
      endif
      if ( posinp('$energy',cpass) ) then
           call cardin(card)
      endif
      call fparr (card,'energies',z(energy),nen,' ')    
      call vcpl(card,cpass,z(vcij),nc,prntcp)
c
c       loop on energies
c
      eloc=energy
      do 20 i=1,nen
c
         call openc (z(ec),z(eloc),nc,nopen,nclosd)
c
c
c       compute the r-matrix at rbeg
c
         call rmtrx (z(rmat),z(srfprj),z(eig),z(eloc),nc,neig,
     1               bypass,prntr)
c
c      calculate the k-matrix, cross section and eigenphases at the
c      r-matrix boundary.  since no propagation has been done we are
c      in the asymptotic, not adiabatic basis, and no transformation
c      is needed.
         write( iout,2) rbeg
c
         call kmtrx(z(rmat),z(kmat),z(coef),z(coef),z(tmtrx),z(cross),
     1              z(y0),z(dy0),z(y1),z(dy1),z(j),z(jp),z(y),z(yp),
     2              z(wron),z(scr),a(ipvt),a(lval),z(ec),z(eloc),
     3              rbeg,nc,nopen,prntk)

         call egnpse(z(kmat),z(coef),psesum,z(scr1),nc,nopen,
     1               z(eloc),prntep)
c
c       propagate the r-matrix from r=rbeg to r=rfinal
c
         if (rfinal.gt.rbeg) then
             if (ief.eq.0) then
                 iex=1
             endif   
             if (ief.eq.1) then
                 iex=i
             endif
             write(iout,3) rfinal
             call rmtxp (z(vmat),z(vcij),range,z(tev),z(ec),z(enrg),
     1                   a(ipvt),z(tp),z(tsave),z(pdiag),z(rmat),
     2                   z(scr1),z(eloc),nc,itri,iex,prntr,
     3                   rbeg,rfinal,stpmin,stpmin,stpmax,beta,prntv,
     4                   npv,prntr,npr,debug,pottyp)
c
             call kmtrx(z(rmat),z(kmat),z(coef),z(coef),z(tmtrx),
     1                  z(cross),z(y0),z(dy0),z(y1),z(dy1),z(j),z(jp),
     2                  z(y),z(yp),z(wron),z(scr),a(ipvt),a(lval),z(ec),
     3                  z(eloc),rfinal,nc,nopen,prntk)
             call egnpse(z(kmat),z(coef),psesum,z(scr1),nc,nopen,
     1                   z(eloc),prntep)
         endif
         eloc=eloc+1
   20 continue
      call iosys ('endfile vdiag on scratch',0,0,0,' ')
      call iosys ('rewind all on scratch read-and-write',0,0,0,' ')
      call iosys ('close scratch',0,0,0,' ')
      call iosys ('destroy scratch',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
c
 1    format(/,1x,'number of channels = ',i3,1x,'r-matrix radius = ',
     1         e15.8,/,1x,'final r value =      ',e15.8,1x,
     2                    'number of energies = ',i3,/,1x,
     3                    'residual atomic charge = ',e15.8)
    2 format (/20x,'initial scattering information at r = ',e15.8) 
    3 format (/,20x,'r-matrix propagation to r = ',e15.8,/) 
      end

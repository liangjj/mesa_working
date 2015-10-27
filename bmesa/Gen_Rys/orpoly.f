*deck orpoly.f 
c***begin prologue     orpoly
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             orpoly
c***purpose            driver for orthogonal functions and generalized gauss
c***                   quadrature.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       orpoly
      program orpoly
c
      implicit integer (a-z)
      character*4096 ops
      character*2 itoc
      character*80 cpass, chrkey, phrse, typwt, refwt, title
      character*16 typint
      character*800 card
      character*128 filkohn
      logical dollar, logkey
      logical prn, switch
      real*8 edge, alpha, beta, fpkey, grid, mu, delalf
      real*8 alf, bet, refalf, refbet, alfmax, fac, tmpalf
      dimension prn(10), edge(2)
      common/io/inp, iout      
      data on /0/
      pointer (pquad,grid(1))
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "kohn filename" from rwf',-1,0,
     1             0,filkohn)
      call iosys ('open kohn as unknown',0,0,0,filkohn)
      prn(1)=logkey(ops,'print=m7500=reference-points-and-weights',
     1              .false.,' ')
      prn(2)=logkey(ops,'print=m7500=ratio-weight',.false.,' ')
      prn(3)=logkey(ops,'print=m7500=lanczos-coefficients',.false.,' ')
      prn(4)=logkey(ops,'print=m7500=points-and-weights',.false.,' ')
      prn(5)=logkey(ops,'print=m7500=summary',.false.,' ')
      switch=logkey(ops,'switch=on',.false.,' ')
      alfmax=fpkey(ops,'rys-max-alpha',50.d0,' ')
      typint=chrkey(ops,'interval-type','finite',' ')
      typwt=chrkey(ops,'weight-type','legendre',' ')
      nrys=intkey(ops,'number-of-runs',1,' ')
      write(iout,2) nrys, typint, typwt
      if(typint.eq.'finite') then
         nfix=intkey(ops,'number-of-fixed-points',0,' ')
         call fparr(ops,'end-points',edge,2,' ')
         write(iout,3) (edge(i),i=1,2), nfix
      endif	 
      if(typwt.eq.'jacobi'.or.typwt.eq.'laguerre') then
         alpha=fpkey(ops,'alpha',0.d0,' ')
         beta=fpkey(ops,'beta',0.d0,' ')
         write(iout,4) alpha, beta
      elseif(typwt.eq.'rys') then
         alpha=fpkey(ops,'initial-rys-alpha',0.d0,' ')      
         delalf=fpkey(ops,'rys-step-in-alpha',1.d0,' ')      
         write(iout,5) alpha, delalf
      endif
      alf=alpha
      bet=beta
      mxnord=0
      do 10 i=1,nrys
         if ( dollar('$run-'//itoc(i),card,cpass,inp) ) then
              nord=intkey(card,'polynomial-n',10,' ')	      
              mxnord=max(mxnord,nord)
              refwt=chrkey(card,'reference-weight','legendre',' ')
              nquad=intkey(card,'reference-quadrature-size',2*nord,' ')
              if(refwt.eq.'jacobi'.or.refwt.eq.'laguerre') then
                 refalf=fpkey(card,'reference-alpha',0.d0,' ')
                 refbet=fpkey(card,'reference-beta',0.d0,' ')
              endif	                   
              tmpalf=alf
c
c             if alpha is large use a scaled integral and hermite quadrature
c             as a reference
c
              if(switch) then
                 write(iout,*) 'switch is on'
                 if(typwt.eq.'rys') then
                    if(alf.ge.alfmax) then
                       refwt='hermite'
                       nquad=3*nord
                       fac=.5d0/sqrt(alf)
                       alf=1.d0
                       refalf=alf
                    endif
                 endif
              endif
	      write(iout,6) nord, nquad, refwt, refalf, refbet, tmpalf     
         endif
         pt=1
         wt=pt+nquad
         ply=wt+nquad
         arg=ply+nquad*(nord+1)
         a=arg+nquad
         b=a+nord+1
         rwt=b+nord+1
         scr=rwt+nquad
         vec=scr+nquad
         need=wpadti(vec+nquad*(nord+1))
         if(on.eq.0) then
            call memory(need,pquad,ngot,'quad',0)
            on=1
         else
            if(need.gt.ngot) then
               call memory(-ngot,pquad,idum,'quad',idum)
               call memory(need,pquad,ngot,'quad',0)
            endif
         endif
         call gauwpt(refwt,nquad,refalf,refbet,0,edge,grid(scr),
     1               grid(pt),grid(wt),prn(1))
         if(typwt.eq.'rys'.and.refwt.eq.'hermite') then
            call vscale(grid(wt),grid(wt),fac,nquad)
         endif
         if(prn(1)) then
            title='reference points'
            call prntfm(title,grid(pt),nquad,1,nquad,1,iout)
            title='reference weights'
            call prntfm(title,grid(wt),nquad,1,nquad,1,iout)
         endif
         call genrwt(grid(rwt),grid(pt),typwt,alf,bet,nquad)
         call genrwt(grid(scr),grid(pt),refwt,refalf,refbet,nquad)
         call vdiv(grid(rwt),grid(rwt),grid(scr),nquad)
         if(prn(2)) then
            title='ratio weight factor'
            call prntfm(title,grid(rwt),nquad,1,nquad,1,iout)
         endif
c
c        generate the recursion coefficients
c
         if(typwt.eq.'rys') then
            call vmul(grid(arg),grid(pt),grid(pt),nquad)
         else
            call copy(grid(pt),grid(arg),nquad)
         endif
         call lancz(grid(vec),grid(arg),grid(a),grid(b),grid(wt),
     1              grid(rwt),mu,grid(scr),nquad,nord)
         call iosys('write integer "quad size-'//itoc(i)
     1               //'" to kohn',1,nord,0,' ')
         call iosys('write real "alpha-'//itoc(i)
     1              //'" to kohn',nord,grid(a),0,' ')
         call iosys('write real "beta-'//itoc(i)
     1              //'" to kohn',nord,grid(b),0,' ')
         if(typwt.eq.'rys') then
            call iosys('write real "X-'//itoc(i)
     1                 //'" to kohn',1,tmpalf,0,' ')
         endif
         if(prn(3)) then
            title='lanczos a coefficients'
            call prntfm(title,grid(a),nord,1,nord,1,iout)
            title='lanczos b coefficients'
            call prntfm(title,grid(b),nord-1,1,nord-1,1,iout)
         endif
         call copy(grid(a),grid(arg),nord)       
         call copy(grid(b),grid(scr),nord)       
         call genq(grid(arg),grid(wt),grid(scr),mu,nfix,edge,nord)
         if(typwt.eq.'rys') then
            call vsqrt(grid(arg),grid(arg),nord)
	    if(refwt.eq.'hermite') then
               call vscale(grid(arg),grid(arg),2.d0*fac,nquad)
	    endif   
         endif	 
         if(prn(4)) then
            title='final points'
            call prntfm(title,grid(arg),nord,1,nord,1,iout)
            title='final weights'
            call prntfm(title,grid(wt),nord,1,nord,1,iout)
         endif
         call iosys('write real "points-'//itoc(i)
     1              //'" to kohn',nord,grid(arg),0,' ')
         call iosys('write real "weights-'//itoc(i)
     1              //'" to kohn',nord,grid(wt),0,' ')
	 alf=tmpalf+delalf
 10   continue   
      call iosys('write integer "max rys quad size" to kohn',
     1            1,mxnord,0,' ')      
      call memory(-ngot,pquad,idum,'quad',idum)
      if(prn(5)) then
         call sumary(nrys,mxnord,typwt)
      endif
      call chainx(0)               
      stop
 1    format(/,20x,'orthogonal polynomial basis function code')
 2    format(/,1x,'number of runs     = ',i3,
     1       /,1x,'interval type      = ',a16,
     1       /,1x,'weight function    = ',a16)
 3    format(/,1x,'left boundary          = ',e15.8,
     1       /,1x,'right boundary         = ',e15.8,
     2       /,1x,'number of fixed points = ',i1)
 4    format(/,1x,'alpha = ',e15.8,
     1       /,1x,'beta = ',e15.8) 
 5    format(/,1x,'initial rys parameter = ',e15.8,
     2       /,1x,'rys stepsize          = ',e15.8)
 6    format(/,1x,'polynomial n                 = ',i3,
     1       /,1x,'size of reference quadrature = ',i3,
     2       /,1x,'reference weight function    = ',a32,
     3       /,1x,'reference alpha              = ',e15.8,
     4       /,1x,'reference beta               = ',e15.8,
     5       /,1x,'rys alpha                    = ',e15.8)
      end



















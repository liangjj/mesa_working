*deck rcbsrt.f 
c***begin prologue     m6234
c***date written       930627   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6234, link 6234, bessel, roots
c***author             schneider, b. i.(nsf)
c***source             m6234
c***purpose            driver for ricatti-bessel root finder
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6234
      program rcbsrt
c
      implicit integer (a-z)
      character*4096 ops
      character*8 cpass
      character*320 card
      character*128 fillam
      character*3 itoc
      logical posinp, logkey, prnt, prntit, prntnt
      common z(1)
      dimension ia(1)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, left, right, del, convg, root, tmp, fpkey
c
      call drum
      write(iout,*)
      write(iout,*) '                         m6234:ricatti-bessel'//
     1                                        ' root finder'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prntit=logkey(ops,'print=m6234=iterations',.false.,' ')
      prntnt=logkey(ops,'print=m6234=zeros',.false.,' ')
      prnt=logkey(ops,'print=m6234',.false.,' ')
      if (prnt) then
          prntnt=.true.
          prntit=.true.
      endif
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$bsroot',cpass) ) then
           call cardin(card)
           left=fpkey(card,'left-boundary',1.d-06,' ')
           if (left.le.1.d-06) then
               left=1.d-06
           endif
           right=fpkey(card,'right-boundary',100.d0,' ')
           step=intkey(card,'number-of-steps',1000,' ')
           mzro=intkey(card,'number-of-zeros',50,' ')
           lmax=intkey(card,'maximum-l-value',2,' ')
           niter=intkey(card,'maximum-number-of-iterations',100,' ')
           convg=fpkey(card,'convergence-criterion',1.d-8,' ')
           del=(right-left)/step
           write(iout,1) mzro, left, right, del, lmax   
      endif
      tmp=lmax*30.d0
      ltop=lmax+sqrt(tmp)
      ltop=max(ltop,lmax)
      lp=ltop+1     
      ioff=1
      do 10 i=1,2
         x=ioff
         j=x+step
         jp=j+lp*step
         y=jp+lp*step
         yp=y+lp*step
         wron=yp+lp*step
         scr=wron+lp
         jdum=scr+2*step
         jpdum=jdum+lp
         ydum=jpdum+lp
         ypdum=ydum+lp
         zval=ypdum+lp
         zdval=zval+mzro+1
         root=zdval+mzro+1
         droot=root+mzro
         words=droot+mzro
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'rcbsrt',0)
         endif
   10 continue
c           make the x grid and compute bessel and derivatives.
      call mkx(z(x),left,del,step)
      call rcbes(z(x),z(j),z(jp),z(y),z(yp),z(wron),z(scr),step,lmax,
     1           ltop,'derivatives',.false.)
      locf=j
      locdf=jp
      do 20 l=0,lmax        
c              locate zeros of bessel
         call intrv(z(locf),z(zval),left,del,step,mzro,cntzro,prntnt)
         lroot=root
         locz=zval
         if (prntit) then         
             write(iout,2)
             write(iout,6)
         endif    
         do 30 nroot=1,cntzro
              call nwtrap(z(locz),z(locz+1),z(lroot),convg,niter,number,
     1                    z(jdum),z(jpdum),z(ydum),z(ypdum),l,ltop,'j')
              if (prntit) then
                  write (iout,3) number, z(lroot)
              endif
              lroot=lroot+1
              locz=locz+1
   30    continue
         call iosys ('write integer "number of j roots '//
     1               'for l='//itoc(l)//'" to lamdat',1,
     2                cntzro,0,' ')
         call iosys ('write real "j roots for '//
     1               'l='//itoc(l)//'" to lamdat',cntzro,
     2                z(root),0,' ')

c              locate zeros of derivative of bessel
         call intrv(z(locdf),z(zdval),left,del,step,mzro,cntzro,prntnt)
         locz=zdval
         lroot=droot
         if (prntit) then         
             write(iout,4)
             write(iout,6)
         endif    
         do 40 nroot=1,cntzro
              call nwtrap(z(locz),z(locz+1),z(lroot),convg,niter,number,
     1                    z(jdum),z(jpdum),z(ydum),z(ypdum),l,ltop,'jp')
              if (prntit) then
                  write (iout,3) number, z(lroot)
              endif
              lroot=lroot+1
              locz=locz+1
   40    continue
         call iosys ('write integer "number of jp roots '//
     1               'for l='//itoc(l)//'" to lamdat',1,
     2                cntzro,0,' ')
         call iosys ('write real "jp roots for '//
     1               'l='//itoc(l)//'" to lamdat',cntzro,
     2                z(droot),0,' ')
         write (iout,5) l
         locf=locf+step
         locdf=locdf+step
   20 continue       
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
    1 format(/,5x,'search for ',i3,' zeros in range [',e15.8,',',
     1            e15.8' ]',/,5x,'function evaluation every ',
     2            e15.8,' bohr',/5x,'maximum l value ',i3)
    2 format(/,5x,'searching on ricatti-bessel')
    3 format (15x,i3,22x,e15.8)
    4 format(/,5x,'searching on derivative of ricatti-bessel')
    5 format(/,25x,'finished l = ',i3)
    6 format(/,5x,'no. newton-raphson iterations',8x,'value of root')      
      end


*deck @(#)cbasis.f
c***begin prologue     cbasis
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           cbasis, link 7000
c***author             schneider, barry (nsf)
c***source             m7000
c***description        fills up arrays with numerical values of
c***                   basis functions and the effect of the
c***                   laplacian on the complex gaussian basis.
c***                   note that angular momentum is included.
c
c***references       
c
c***routines called
c***end prologue       cbasis
      subroutine cbasis(fns,ddfns,lval,mval,alpha,power,pt,nsym,npts,
     1                  ntot,dimsym,prnt,type)
      implicit integer (a-z)
      complex *16 fns, ddfns, alpha, afac, afac1, fac 
      real *8 pt, fac1
      character *24 cpass
      character *1600 card
      logical prnt, posinp
      character *8 type
      dimension fns(npts,ntot), ddfns(npts,ntot), pt(npts)
      dimension lval(ntot), mval(ntot), power(ntot), alpha(ntot)
      dimension nsym(dimsym)
      common/io/inp,iout
      write(iout,1) type
      call iosys ('write integer "number complex '//
     1            'functions" to atomci',1,ntot,0,' ')   
      call izero(lval,ntot)
      call izero(mval,ntot)
      call izero(power,ntot)
      call czero(alpha,ntot)
      call izero(nsym,dimsym)
      icnt=0
      if ( posinp('$s-type-c',cpass) ) then
           call cardin(card)
           nsym(1)=intkey(card,'number',1,' ')
           call iosys ('write integer "number s type '//
     1                 'complex functions" to atomci',1,nsym(1),0,' ')
           ibeg=icnt
           do 10 i=1,nsym(1)
              icnt=icnt+1
              lval(icnt)=0
              mval(icnt)=0
   10      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(1),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(1),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$p(-1)-type-c',cpass) ) then
           call cardin(card)
           nsym(2)=intkey(card,'number',1,' ')
           call iosys ('write integer "number p(-1) type '//
     1                 'complex functions" to atomci',1,nsym(2),0,' ')
           ibeg=icnt
           do 20 i=1,nsym(2)
              icnt=icnt+1
              lval(icnt)=1
              mval(icnt)=-1
   20      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(2),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(2),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$p(0)-type-c',cpass) ) then
           call cardin(card)
           nsym(3)=intkey(card,'number',1,' ')
           call iosys ('write integer "number p(0) type '//
     1                 'complex functions" to atomci',1,nsym(3),0,' ')   
           ibeg=icnt
           do 30 i=1,nsym(3)
              icnt=icnt+1
              lval(icnt)=1
              mval(icnt)=0
   30      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(3),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(3),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$p(1)-type-c',cpass) ) then
           call cardin(card)
           nsym(4)=intkey(card,'number',1,' ')
           call iosys ('write integer "number p(1) type '//
     1                 'complex functions" to atomci',1,nsym(4),0,' ')   
           ibeg=icnt
           do 40 i=1,nsym(4)
              icnt=icnt+1
              lval(icnt)=1
              mval(icnt)=1
   40      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(4),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(4),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$d(-2)-type-c',cpass) ) then
           call cardin(card)
           nsym(5)=intkey(card,'number',1,' ')
           call iosys ('write integer "number d(-2) type '//
     1                 'complex functions" to atomci',1,nsym(5),0,' ')   
           ibeg=icnt
           do 50 i=1,nsym(5)
              icnt=icnt+1
              lval(icnt)=2
              mval(icnt)=-2
   50      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(5),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(5),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$d(-1)-type-c',cpass) ) then
           call cardin(card)
           nsym(6)=intkey(card,'number',1,' ')
           call iosys ('write integer "number d(-1) type '//
     1                 'complex functions" to atomci',1,nsym(6),0,' ')   
           ibeg=icnt
           do 60 i=1,nsym(6)
              icnt=icnt+1
              lval(icnt)=2
              mval(icnt)=-1
   60      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(6),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(6),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$d(0)-type-c',cpass) ) then
           call cardin(card)
           nsym(7)=intkey(card,'number',1,' ')
           call iosys ('write integer "number d(0) type '//
     1                 'complex functions" to atomci',1,nsym(7),0,' ')   
           ibeg=icnt
           do 70 i=1,nsym(7)
              icnt=icnt+1
              lval(icnt)=2
              mval(icnt)=0
   70      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(7),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(7),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$d(1)-type-c',cpass) ) then
           call cardin(card)
           nsym(8)=intkey(card,'number',1,' ')
           call iosys ('write integer "number d(1) type '//
     1                 'complex functions" to atomci',1,nsym(8),0,' ')   
           ibeg=icnt
           do 80 i=1,nsym(8)
              icnt=icnt+1
              lval(icnt)=2
              mval(icnt)=1
   80      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(8),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(8),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$d(2)-type-c',cpass) ) then
           call cardin(card)
           nsym(9)=intkey(card,'number',1,' ')
           call iosys ('write integer "number d(2) type '//
     1                 'complex functions" to atomci',1,nsym(9),0,' ')   
           ibeg=icnt
           do 90 i=1,nsym(9)
              icnt=icnt+1
              lval(icnt)=2
              mval(icnt)=2
   90      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(9),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(9),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$f(-3)-type-c',cpass) ) then
           call cardin(card)
           nsym(10)=intkey(card,'number',1,' ')
           call iosys ('write integer "number f(-3) type '//
     1                 'complex functions" to atomci',1,nsym(10),0,' ')   
           ibeg=icnt
           do 100 i=1,nsym(10)
              icnt=icnt+1
              lval(icnt)=3
              mval(icnt)=-3
  100      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(10),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(10),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$f(-2)-type-c',cpass) ) then
           call cardin(card)
           nsym(11)=intkey(card,'number',1,' ')
           call iosys ('write integer "number f(-2) type '//
     1                 'complex functions" to atomci',1,nsym(11),0,' ')   
           ibeg=icnt
           do 200 i=1,nsym(11)
              icnt=icnt+1
              lval(icnt)=3
              mval(icnt)=-2
  200      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(11),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(11),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$f(-1)-type-c',cpass) ) then
           call cardin(card)
           nsym(12)=intkey(card,'number',1,' ')
           call iosys ('write integer "number f(-1) type '//
     1                 'complex functions" to atomci',1,nsym(12),0,' ')   
           ibeg=icnt
           do 300 i=1,nsym(12)
              icnt=icnt+1
              lval(icnt)=3
              mval(icnt)=-1
  300      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(12),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(12),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$f(0)-type-c',cpass) ) then
           call cardin(card)
           nsym(13)=intkey(card,'number',1,' ')
           call iosys ('write integer "number f(0) type '//
     1                 'complex functions" to atomci',1,nsym(13),0,' ')   
           ibeg=icnt
           do 400 i=1,nsym(13)
              icnt=icnt+1
              lval(icnt)=3
              mval(icnt)=0
  400      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(13),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(13),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$f(1)-type-c',cpass) ) then
           call cardin(card)
           nsym(14)=intkey(card,'number',1,' ')
           call iosys ('write integer "number f(1) type '//
     1                 'complex functions" to atomci',1,nsym(14),0,' ')   
           ibeg=icnt
           do 500 i=1,nsym(14)
              icnt=icnt+1
              lval(icnt)=3
              mval(icnt)=1
  500      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(14),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(14),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$f(2)-type-c',cpass) ) then
           call cardin(card)
           nsym(15)=intkey(card,'number',1,' ')
           call iosys ('write integer "number f(2) type '//
     1                 'complex functions" to atomci',1,nsym(15),0,' ')   
           ibeg=icnt
           do 600 i=1,nsym(15)
              icnt=icnt+1
              lval(icnt)=3
              mval(icnt)=2
  600      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(15),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(15),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$f(3)-type-c',cpass) ) then
           call cardin(card)
           nsym(16)=intkey(card,'number',1,' ')
           call iosys ('write integer "number f(3) type '//
     1                 'complex functions" to atomci',1,nsym(16),0,' ')   
           ibeg=icnt
           do 700 i=1,nsym(16)
              icnt=icnt+1
              lval(icnt)=3
              mval(icnt)=3
  700      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(16),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(16),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(-4)-type-c',cpass) ) then
           call cardin(card)
           nsym(17)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(-4) type '//
     1                 'complex functions" to atomci',1,nsym(17),0,' ')   
           ibeg=icnt
           do 800 i=1,nsym(17)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=-4
  800      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(17),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(17),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(-3)-type-c',cpass) ) then
           call cardin(card)
           nsym(18)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(-3) type '//
     1                 'complex functions" to atomci',1,nsym(18),0,' ')   
           ibeg=icnt
           do 900 i=1,nsym(18)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=-3
  900      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(18),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(18),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(-2)-type-c',cpass) ) then
           call cardin(card)
           nsym(19)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(-2) type '//
     1                 'complex functions" to atomci',1,nsym(19),0,' ')   
           ibeg=icnt
           do 1000 i=1,nsym(19)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=-2
 1000      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(19),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(19),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(-1)-type-c',cpass) ) then
           call cardin(card)
           nsym(20)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(-1) type '//
     1                 'complex functions" to atomci',1,nsym(20),0,' ')   
           ibeg=icnt
           do 1100 i=1,nsym(20)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=-1
 1100      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(20),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(20),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(0)-type-c',cpass) ) then
           call cardin(card)
           nsym(21)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(0) type '//
     1                 'complex functions" to atomci',1,nsym(21),0,' ')   
           ibeg=icnt
           do 1200 i=1,nsym(21)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=0
 1200      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(21),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(21),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(1)-type-c',cpass) ) then
           call cardin(card)
           nsym(22)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(1) type '//
     1                 'complex functions" to atomci',1,nsym(22),0,' ')   
           ibeg=icnt
           do 1300 i=1,nsym(22)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=1
 1300      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(22),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(22),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(2)-type-c',cpass) ) then
           call cardin(card)
           nsym(23)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(2) type '//
     1                 'complex functions" to atomci',1,nsym(23),0,' ')   
           ibeg=icnt
           do 1400 i=1,nsym(23)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=2
 1400      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(23),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(23),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(3)-type-c',cpass) ) then
           call cardin(card)
           nsym(24)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(3) type '//
     1                 'complex functions" to atomci',1,nsym(24),0,' ')   
           ibeg=icnt
           do 1500 i=1,nsym(24)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=3
 1500      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(24),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(24),' ')
      endif
      if (icnt.eq.ntot) go to 5000
      if ( posinp('$g(4)-type-c',cpass) ) then
           call cardin(card)
           nsym(25)=intkey(card,'number',1,' ')
           call iosys ('write integer "number g(4) type '//
     1                 'complex functions" to atomci',1,nsym(25),0,' ')   
           ibeg=icnt
           do 1600 i=1,nsym(25)
              icnt=icnt+1
              lval(icnt)=4
              mval(icnt)=4
 1600      continue     
           call intarr(card,'powers',power(ibeg+1),nsym(25),' ')
           call fparr(card,'exponents',alpha(ibeg+1),2*nsym(25),' ')
      endif
 5000 continue
      call iosys ('write integer "complex symmetry list" to atomci',
     1             dimsym,nsym,0,' ')
      if (prnt) then
          write(iout,2)
          do 2000 i=1,ntot
             write(iout,3) lval(i), mval(i), power(i), alpha(i)
 2000     continue
      endif
      call iosys('write integer "l values of complex orbitals" to'
     1           //' atomci',ntot,lval,0,' ')
      call iosys('write integer "m values of complex orbitals" to'
     1           //' atomci',ntot,mval,0,' ')
      call iosys('write integer "powers of complex orbitals" to'
     1           //' atomci',ntot,power,0,' ')
      call iosys('write real "exponents of complex orbitals" to'
     1           //' atomci',2*ntot,alpha,0,' ')
c----------------------------------------------------------------------c
c           compute the functions on the choosen grid and the value    c
c                                   of                                 c
c                      ( -1/2 Del**2 ) on the function                 c
c           this includes the l*(l+1)/r*r angular momentum due to      c
c           the implicit spherical harmonic multiplying the radial     c
c           function.                                                  c
c----------------------------------------------------------------------c   
      if (type.eq.'gaussian') then
          do 3000 i=1,ntot
             afac=4.d0*alpha(i)*alpha(i)
             nfac=power(i)*(power(i)+1)
             afac1=-2.d0*alpha(i)*(2*power(i)+3)
             do 4000 j=1,npts
                fac=exp(-alpha(i)*pt(j)*pt(j))
                fac1=pt(j)**power(i)
                fns(j,i)=fac*fac1
                ddfns(j,i)=afac*pt(j)*pt(j)+afac1
                if(power(i).gt.0) then
                   ddfns(j,i)=ddfns(j,i)+nfac/(pt(j)*pt(j))
                endif
                ddfns(j,i)=ddfns(j,i)*fac*fac1
                ddfns(j,i)=-.5d0 * ( ddfns(j,i) -
     1                       lval(i)*(lval(i)+1)*fns(j,i)/(pt(j)*pt(j)))
 4000        continue
 3000     continue   
      elseif (type.eq.'slater') then
          do 6000 i=1,ntot
             afac=alpha(i)*alpha(i)
             nfac=power(i)*(power(i)+1)
             afac1=-alpha(i)*(2*power(i)+2)
             do 7000 j=1,npts
                fac=exp(-alpha(i)*pt(j))
                fac1=pt(j)**power(i)
                fns(j,i)=fac*fac1
                ddfns(j,i)=afac+afac1/pt(j)
                if(power(i).gt.0) then
                   ddfns(j,i)=ddfns(j,i)+nfac/(pt(j)*pt(j))
                endif
                ddfns(j,i)=ddfns(j,i)*fac*fac1
                ddfns(j,i)=-.5d0 * ( ddfns(j,i) -
     1                       lval(i)*(lval(i)+1)*fns(j,i)/(pt(j)*pt(j)))
 7000        continue
 6000     continue   
      endif
      return
    1 format(//,20x,'using',1x,a8,1x,'type orbitals')
    2 format(//,10x,'complex orbitals',//,6x,' l ',5x,' m ',3x,
     1       ' power ',20x,'   exponent   ',/)
    3 format(5x,i3,5x,i3,5x,i3,10x,f15.8,2x,f15.8)   
      end

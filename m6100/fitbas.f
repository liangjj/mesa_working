*deck @(#)fitbas.f
c***begin prologue     fitbas
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, link 6080
c***author             schneider, barry (lanl)
c***source             m6080
c***purpose            internal fitting functions on grid
c***description        fills up arrays with numerical values of
c***                   functions, second derivatives and matching information.
c*** 
c
c***references       
c
c***routines called
c***end prologue       fitbas
      subroutine fitbas(fns,ddfns,alpha,power,pt,mtch,rmatch,
     1                  npts,ntot,type)
      implicit integer (a-z)
      real *8 fns, ddfns, pt, alpha, afac, afac1, fac, fac1
      real *8 mtch, rmatch
      character *1600 card
      character *16 type
      dimension fns(npts,ntot), ddfns(npts,ntot), pt(npts)
      dimension power(ntot), alpha(ntot), mtch(3,ntot)
      common/io/inp,iout
      call cardin(card)
      call intarr(card,'powers',power,ntot,' ')
      call fparr(card,'exponents',alpha,ntot,' ')
      if (type.eq.'gaussian') then
          do 10 i=1,ntot
             afac=4.d0*alpha(i)*alpha(i)
             nfac=power(i)*(power(i)-1)
             afac1=-2.d0*alpha(i)*(2*power(i)+1)
             do 20 j=1,npts
                fac=exp(-alpha(i)*pt(j)*pt(j))
                fac1=pt(j)**power(i)
                fns(j,i)=fac*fac1
                ddfns(j,i)=afac*pt(j)*pt(j)+afac1
                if(power(i).gt.1) then
                   ddfns(j,i)=ddfns(j,i)+nfac/(pt(j)*pt(j))
               endif
               ddfns(j,i)=ddfns(j,i)*fac*fac1
   20        continue
             fac=exp(-alpha(i)*rmatch*rmatch)
             fac1=rmatch**power(i)
             mtch(1,i)=fac*fac1
             mtch(2,i)=-2.d0*alpha(i)*rmatch
             mtch(3,i)=afac*rmatch*rmatch+afac1
             if (power(i).ne.0) then
                 mtch(2,i)=mtch(2,i)+power(i)/rmatch
             endif
             mtch(2,i)=mtch(2,i)*fac*fac1
             if(power(i).gt.1) then
                mtch(3,i)=mtch(3,i)+nfac/(rmatch*rmatch)
             endif
             mtch(3,i)=mtch(3,i)*fac*fac1
   10     continue    
      elseif(type.eq.'polynomial') then
          do 30 i=1,ntot
             do 40 j=1,npts
                fns(j,i)=pt(j)**power(i)
                ddfns(j,i)=0.d0
                if (power(i).gt.1) then
                    ddfns(j,i)=power(i)*(power(i)-1)*
     1                         fns(j,i)/(pt(j)*pt(j))
                endif
   40        continue     
             mtch(1,i)=rmatch**power(i)
             mtch(2,i)=0.d0
             if(power(i).gt.0) then
                mtch(2,i)=power(i)*rmatch**(power(i)-1)
             endif
             mtch(3,i)=0.d0
             if(power(i).gt.1) then
                mtch(3,1)=power(i)*(power(i)-1)*rmatch**(power(i)-2)
             endif
   30     continue         

      endif
      return
      end


















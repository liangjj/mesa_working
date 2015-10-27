*deck @(#)basis.f
c***begin prologue     basis
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, link 6080
c***author             schneider, barry (lanl)
c***source             m6080
c***purpose            internal basis functions on grid
c***description        fills up arrays with numerical values of
c***                   basis functions and second derivatives.
c*** 
c
c***references       
c
c***routines called
c***end prologue       basis
      subroutine basis(fns,ddfns,alpha,power,pt,npts,ntot)
      implicit integer (a-z)
      real *8 fns, ddfns, pt, alpha, afac, afac1, fac, fac1
      character *1600 card
      dimension fns(npts,ntot), ddfns(npts,ntot), pt(npts)
      dimension power(ntot), alpha(ntot)
      common/io/inp,iout
      call cardin(card)
      call intarr(card,'powers',power,ntot,' ')
      call fparr(card,'exponents',alpha,ntot,' ')
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
   20    continue
   10 continue    
      return
      end

*deck %W%  %G%
      subroutine mcflag(xlag,temp,c,nbf,nob,nco,nao)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
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
c
      implicit integer(a-z)
      real*8 xlag(*),c(*),temp(*)
c
      noc=nco+nao
      ntot=noc*nob
c
      call ebtc(temp,c,xlag,nob,nbf,noc)
c
      do 10 i=1,ntot
         xlag(i)=temp(i)
 10   continue
c
      return
      end

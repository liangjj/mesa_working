*deck %W%  %G%
      subroutine mcabix(c,xjk,nbf,nco,navo,temp,tv,nfiv,
     $     factiv,ikt,tfile,buf)
c
c***begin prologue     mcabix
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
c----------------------------------------------------------------
c
c     this program performs a limited transformation
c     of the an exchange integral block
c       one of the indices is transformed over the core orbitals.
c       the other index is transformed over the active + virtual
c       orbitals. the result is stored on unit nfiv
c
c       c     molecular orbitals
c       xjk    exchange integral block
c       temp  temporary storage for the transformed fock operator
c       tv    temporary storage .. vector of length  nbf
c       nco   the number of core orbitals
ccc....nope     navo  the number of active + virtual orbitals
c       navo  the total number of orbitals  nob
c       nfiv  output unit
c
c----------------------------------------------------------------
c
c***references
c
c***routines called    (none)
c
c***end prologue       mcabix
c
      implicit real*8 (a-h,o-z)
c
      real*8 buf(*)
      character*(*) tfile
      dimension c(nbf,2),xjk(nbf,2),temp(navo,2),tv(2)
c
      common /io/ inp,iout
      common / number / zero,pt5,one,two,four,eight
      common / nhexpk / nhexab
c
      call ebc(temp,xjk,c,nbf,nbf,nco)
      call ebtc(buf,c,temp,navo,nbf,nco)
c
      if(ikt.ne.0)go to 100
      ntot=navo*nco
      do 11 i=1,ntot
         buf(i)=buf(i)*factiv
 11   continue
      return
c
 100  continue
      call ebtc(temp,xjk,c,nbf,nbf,nco)
      call apbtc(buf,c,temp,navo,nbf,nco)
      ntot=navo*nco
      do 10 i=1,ntot
         buf(i)=buf(i)*factiv
 10   continue
c
c
      return
      end

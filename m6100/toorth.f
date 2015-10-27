*deck @(#)toorth.f
c***begin prologue     toorth
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           toorth, link 6080
c***author             schneider, barry (lanl)
c***source             m6080
c***purpose            free functions to set orthogonal to bound
c***                   functions
c*** 
c
c***references       
c
c***routines called
c***end prologue       toorth
      subroutine toorth(fns,ddfns,cfn,ddcfn,wt,sbbp,sfbp,sfbt,
     1                  scrc,tempc,nbf,count,npts)
      implicit integer (a-z)
      real *8  fns, ddfns,sbbp, wt
      complex *16 cfn, ddcfn, scrc, sfbp, sfbt, tempc
      dimension fns(npts,nbf), ddfns(npts,nbf), wt(npts)
      dimension cfn(npts,2), ddcfn(npts,2)
      dimension scrc(npts,2), sfbp(2,nbf), sbbp(nbf,count)
      dimension sfbt(2,count), tempc(2,nbf)
      common/io/inp,iout
      do 10 i=1,2
         do 20 j=1,npts
            scrc(j,i)=wt(j)*cfn(j,i)
   20    continue     
   10 continue
      call ecbtc(sfbp,scrc,fns,2,npts,nbf)
      call ecbc(sfbt,sfbp,sbbp,2,nbf,count)
      call ecbct(tempc,sfbt,sbbp,2,count,nbf)
      call ambcct(cfn,fns,tempc,npts,nbf,2)
      call ambcct(ddcfn,ddfns,tempc,npts,nbf,2)
      return
      end


















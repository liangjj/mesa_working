*deck @(#)match.f
c***begin prologue     match
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           match, link 6080
c***author             schneider, barry (lanl)
c***source             m6080
c***purpose            match function, first and second derivative
c***                   to outgoing wave
c*** 
c
c***references       
c
c***routines called
c***end prologue       match
      subroutine match(fns,ddfns,cfn,ddcfn,diff,pt,mtch,mtrx,rhs,energy,
     1                 rmatch,last,ipvt,npts,nbfn)
      implicit integer (a-z)
      real *8 fns, ddfns, mtch, energy, rmatch, k, pt, sqk
      complex *16 cfn, ddcfn, ck, mtrx, rhs, diff
      dimension fns(npts,nbfn), ddfns(npts,nbfn)
      dimension cfn(npts,2), ddcfn(npts,2)
      dimension mtch(3,nbfn), mtrx(3,*), rhs(*), pt(npts), ipvt(*)
      dimension diff(npts)
      common/io/inp,iout
      k=sqrt(2.d0*energy)
      sqk=1.d0/sqrt(k)
      ck=cmplx(0.d0,k)
      rhs(1)=exp(ck*rmatch)
      rhs(2)=ck*rhs(1)
      if (nbfn.gt.2) then
          rhs(3)=ck*rhs(2)
      endif
      limit=nbfn
      call czero(cfn,2*npts)
      call czero(ddcfn,2*npts)
      do 10 i=1,limit
         do 20 j=1,limit
            mtrx(i,j)=mtch(i,j)
   20    continue
   10 continue
      call cgefa(mtrx,3,limit,ipvt,info)
      call cgesl(mtrx,3,limit,ipvt,rhs,0)
      do 30 i=1,limit
         do 40 j=1,last
            cfn(j,2)=cfn(j,2)+rhs(i)*fns(j,i)
            ddcfn(j,2)=ddcfn(j,2)+rhs(i)*ddfns(j,i)
            cfn(j,1)=imag(cfn(j,2))
            ddcfn(j,1)=imag(ddcfn(j,2))
            diff(j)=ddcfn(j,2)-ck*ck*cfn(j,2)
   40    continue
   30 continue
      do 50 i=1,last
         cfn(i,1)=cfn(i,1)*sqk
         ddcfn(i,1)=ddcfn(i,1)*sqk 
         cfn(i,2)=cfn(i,2)*sqk
         ddcfn(i,2)=ddcfn(i,2)*sqk 
   50 continue     
      do 60 i=last+1,npts
         cfn(i,1)=sqk*cmplx(sin(k*pt(i)),0.d0)
         cfn(i,2)=sqk*exp(ck*pt(i))
         ddcfn(i,1)=-k*k*cfn(i,1)
         ddcfn(i,2)=ck*ck*cfn(i,2)
   60 continue
      return
      end


















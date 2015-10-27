*deck @(#)hsmat.f
c***begin prologue     hsmat
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, link 6080
c***author             schneider, barry (lanl)
c***source             m6090
c***purpose            kohn matrix elements for numerical functions
c***description        matrix elements of ( h - E ) are computed for
c***                   bound-bound, bound-free and free-free functions.
c***                   they are computed by radial quadrature
c*** 
c
c***references       
c
c***routines called
c***end prologue       hsmat
      subroutine hsmat(fns,ddfns,cfn,ddcfn,v,energy,wt,hbbp,hfbp,hffp,
     1                 sbbp,sfbp,hbbt,hfbt,scrb,scrc,temp,tempc,nbfn,
     2                 npts,count,prntfb,prntff)
      implicit integer (a-z)
      real *8 fns, ddfns, v, energy, wt, hbbp, sbbp, hbbt
      real *8 scrb, temp
      complex *16 cfn, ddcfn, hfbp, hffp, sfbp, hfbt, sfbt, scrc, tempc 
      character *80 title
      logical prntfb, prntff
      dimension fns(npts,nbfn), ddfns(npts,nbfn), cfn(npts,2)
      dimension ddcfn(npts,2), v(npts), wt(npts), hbbp(nbfn,nbfn)
      dimension hfbp(2,nbfn), hffp(2,2), sbbp(nbfn,nbfn)
      dimension sfbp(2,nbfn), hbbt(*), hfbt(*)
      dimension scrb(npts,nbfn), scrc(npts,2), temp(nbf,nbf)
      dimension tempc(nbf,2)
      common /io/ inp, iout
c
c               bound-bound matrix elements ( e - h )
c
      do 10 i=1,nbfn
         hbbt(i,i)=energy-hbbp(i,i)
   10 continue     
c
c               free-bound matrix elements  ( e - h )
c
      do 20 i=1,nbfn
         do 30 j=1,npts
            scrb(j,i)=(.5d0*ddfns(j,i)+(energy-v(j))*fns(j,i))
     1                            *wt(j)
   30    continue
   20 continue
      call ecbtc(hfbp,cfn,scrb,2,npts,nbfn)
      if(prntfb) then
         title='free-bound primitive hamiltonian'
         call prntcmn(title,hfbp,2,nbfn,2,nbfn,iout,'e')
      endif
c
c               free-bound overlaps
c
      call ecbtc(sfbp,cfn,scrb,2,npts,nbfn)
      if(prntfb) then
         title='free-bound primitive overlap'
         call prntcmn(title,sfbp,2,nbfn,2,nbfn,iout,'e')
      endif
c
c               free-free matrix elements  ( e - h )   
c
      do 50 i=1,2
         do 60 j=1,npts
            scrc(j,i)=(.5d0*ddcfn(j,i)+(energy-v(j))*cfn(j,i))
     1                               *wt(j)
   60    continue
   50 continue
      call cebtc(hffp,cfn,scrc,2,npts,2)
c
      if (prntff) then
          title='free-free hamiltonian'
          call prntcmn(title,hffp,2,2,2,2,iout,'e')
      endif
c
c               done with raw matrix elements      
c               transform to orthogonal bound basis
c
c               free-bound
c
c           first transform bound set
      call ecbc(hfbt,hfbp,sbbp,2,nbfn,count)
      call ecbc(sfbt,sfbp,sbbp,2,nbfn,count)
c           now transform to orthogonal free set    
      call amcbc(hfbt,sfbt,hbbt,2,count,count)
      if(prntfb) then
         title='free-bound final hamiltonian'
         call prntcmn(title,hbbt,2,count,2,count,iout,'e')
      endif
c
c               now transform free functions to space
c               orthogonal to bound functions
      call cambct(hffp,sfbt,hfbt,2,count,2)
      call cambct(hffp,hfbt,sfbt,2,count,2)
      call ecbc(tempc,sfbt,hbbt,2,count,count)
      call capbct(hffp,tempc,sfbt,2,count,2)
      if(prntff) then
         title='free-free final hamiltonian'
         call prntcmn(title,hffp,2,2,2,2,iout,'e')
      endif
      return
      end







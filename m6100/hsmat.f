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
     1                 sbbp,sfbp,hbbt,hfbt,sfbt,scrb,scrc,tempc,nbf,
     2                 npts,count,prntfb,prntff,north)
      implicit integer (a-z)
      real *8 fns, ddfns, v, energy, wt, hbbp, sbbp, hbbt
      real *8 scrb
      complex *16 cfn, ddcfn, hfbp, hffp, sfbp, hfbt, sfbt, scrc, tempc 
      complex *16 test
      character *80 title
      logical prntfb, prntff, north
      dimension fns(npts,nbf), ddfns(npts,nbf), cfn(npts,2)
      dimension ddcfn(npts,2), v(npts), wt(npts), hbbp(nbf,nbf)
      dimension hfbp(2,nbf), hffp(2,2), sbbp(nbf,count)
      dimension sfbp(2,nbf), hbbt(count,count), hfbt(2,count)
      dimension sfbt(2,count), scrb(npts,nbf), scrc(npts,2)
      dimension tempc(2,nbf)
      common /io/ inp, iout
c
c               bound-bound matrix elements ( h - e )
c
      do 10 i=1,nbf
         do 11 j=1,nbf
            hbbt(i,j)=hbbp(i,j)
   11    continue
         hbbt(i,i)=hbbt(i,i)-energy   
   10 continue     
c
c               free-bound matrix elements  ( h - e )
c
      do 20 i=1,nbf
         do 30 j=1,npts
            scrb(j,i)=(-.5d0*ddfns(j,i)+(v(j)-energy)*fns(j,i))
     1                            *wt(j)
   30    continue
   20 continue
      call ecbtc(hfbp,cfn,scrb,2,npts,nbf)
      if(prntfb) then
         title='free-bound primitive hamiltonian'
         call prntcmn(title,hfbp,2,nbf,2,nbf,iout,'e')
      endif
c
c               free-bound overlaps
c
      do 40 i=1,nbf
         do 50 j=1,npts
            scrb(j,i)=fns(j,i)*wt(j)
   50    continue
   40 continue  
      call ecbtc(sfbp,cfn,scrb,2,npts,nbf)
      if(prntfb) then
         title='free-bound primitive overlap'
         call prntcmn(title,sfbp,2,nbf,2,nbf,iout,'e')
      endif
c
c               free-free matrix elements  ( h - e )   
c
      do 60 i=1,2
         do 70 j=1,npts
            scrc(j,i)=(-.5d0*ddcfn(j,i)+(v(j)-energy)*cfn(j,i))
     1                               *wt(j)
   70    continue
   60 continue
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
c           first transform bound primitive to bound orthonormal set
      call ecbc(hfbt,hfbp,sbbp,2,nbf,count)
      call ecbc(sfbt,sfbp,sbbp,2,nbf,count)
      if(prntfb) then
         title='free-bound final hamiltonian'
         call prntcmn(title,hfbt,2,count,2,count,iout,'e')
         title='free-bound overlap before orthogonalization'
         call prntcmn(title,sfbt,2,count,2,count,iout,'e')
      endif
c
c               now transform all free functions to space
c               orthogonal to bound functions if desired
      if (.not.north) then
          call cambct(hffp,sfbt,hfbt,2,count,2)
          call cambct(hffp,hfbt,sfbt,2,count,2)
          call ecbc(tempc,sfbt,hbbt,2,count,count)
          call capbct(hffp,tempc,sfbt,2,count,2)
          call amcbc(hfbt,sfbt,hbbt,2,count,count)
      endif
      if(prntff) then
         title='free-free final hamiltonian'
         call prntcmn(title,hffp,2,2,2,2,iout,'e')
      endif
      return
      end







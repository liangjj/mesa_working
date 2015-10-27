      subroutine rhscal (cfn,ddcfn,v,rhs,energy,npts,prnt)
      implicit integer(a-z)
      common /io/ inp, iout
      real *8 energy, v 
      complex *16 rhs, cfn, ddcfn
      logical prnt
      character *80 title
      dimension rhs(npts), v(npts), cfn(npts), ddcfn(npts)
      do 10 j=1,npts
         rhs(j)=.5d0*ddcfn(j) + ( energy - v(j) ) * cfn(j)
   10 continue
      if (prnt) then
          title='matrix rhs'
          call prntcmn(title,rhs,npts,1,npts,1,iout,'e')
      endif 
      return
c
      end

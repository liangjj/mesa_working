!===============================================================================!
! This subroutine calculates the electron-electron Coulomb potential of 1/r_12. !
! Note that NO matrix elements of 1/r_12 are computed in the present subroutine.!
! It only calculates the function 1/r_12 itself directly in Cartesian           !
! coordinates and also by using the von Neumann expansion in terms of the P and !
! Q functions. And then compare these two results with each other. The purpose  !
! of calling this subroutine is to get a feeling about how large lmax is        !
! required even only for 1/r_12 itself. This could be useful for computation of !
! matrix element of 1/r_12.                                                     !
!-------------------------------------------------------------------------------!
! Input 'lmax' is the maximum of l, at which the summation over l is truncated. !
! Input coordinates of two electrons are given in prolate spheoidal coordinates.!
!-------------------------------------------------------------------------------!
! Output 'couinv_r12_direct' is obtained by directly calculating it in          !
! coordinates. The output 'couinv_r12_neumann' is obtained by Neumann expansion.!
!===============================================================================!
      subroutine clmb_r12_itself(lmax,xi1,eta1,varphi1,xi2,eta2,varphi2,&
                                 radius_moeq,couinv_r12_direct,couinv_r12_neumann)
      use accuracy
      implicit none
      integer            intent(in) :: lmax
      real(KIND=idp),    intent(in) :: xi1,eta1,varphi1,xi2,eta2,varphi2,radius_moeq
      real(KIND=idp),   intent(out) :: couinv_r12_direct,couinv_r12_neumann
      integer                       :: l,m,nmax,mn
      real(KIND=idp)                :: a,b,c,x,rho1,rho2
      real(KIND=idp)                :: x1,x2,y1,y2,z1,z2,r12,x12,y12,z12
      real(KIND=idp)                :: xi_small,xi_large,p1,p2,rsum,factor,&
                                       frac1,frac2,psmall,qlarge,temp,varphi,varp,vam
      character(LEN=24),            :: recur
      character(LEN=16),            :: directive

      if(lmax <= 0) then
         write(idwrite,*) 'In subroutine : clmb_r12_itself'
         write(idwrite,*) 'lmax <=0 !!!'
         write(idwrite,*) 'lmax =',lmax
         write(idwrite,*) 'It must be > 0 to make a meaningfull &
     & summation in Neumann expansion !!!'
         write(idwrite,*) 'I stopped !!!'
         return
      end if

      recur = 'Miller'
      directive = 'irregular'



      a = 0.5_idp*radius_moeq

! the electron 1:
      rho1 = a*SQRT((xi1*xi1-1.0_idp)*(1.0_idp-eta1*eta1))
      x1 = rho1 * COS(varphi1)
      y1 = rho1 * SIN(varphi1)
      z1 = a*xi1*eta1

! the electron 2:
      rho2 = a*SQRT((xi2*xi2-1.0_idp)*(1.0_idp-eta2*eta2))
      x2 = rho2 * COS(varphi2)
      y2 = rho2 * SIN(varphi2)
      z2 = a*xi2*eta2

      x12 = x1-x2
      y12 = y1-y2
      z12 = z1-z2
      r12 = x12*x12 + y12*y12 + z12*z12

      couinv_r12_direct = 1.0_idp / SQRT(r12)

      xi_small = MIN(xi1,xi2)
      xi_large = MAX(xi1,xi2)


! P_LM(eta):
      directive = 'regular'
      normalize = .false. 
      mmax = lmax
      call legendre_driver(mmax,lmax,recur,directive,normalize,  &
                                 x,fun_legen,error_messg)

! Neumann expansion:
         rsum = 0.0_idp
      do l = 0, lmax
      do m = 0, l

                  vam = 2.0_idp
         if(m==0) vam = 1.0_idp
         varphi = varphi1 - varphi2
         nmax = l - m
         call factorial_min_max(1,nmax,frac1)
         nmax = l + m
         call factorial_min_max(1,nmax,frac2)
         mn = (-1)**m
         temp = dble(mn)
         factor = frac1/frac2
         factor = factor*factor*dble(l+l+1)*temp
         varp = vam * COS(m*varphi)
         call associate_legendre_xip(l,m,xi_small,psmall)
         call associate_legendre_xiq(l,m,xi_large,qlarge)
         call associate_legendre_eta(l,m,eta1,p1)
         call associate_legendre_eta(l,m,eta2,p2)

         temp = factor * psmall * qlarge * p1 * p2 * varp
         rsum = rsum + temp
      end do
      end do

      couinv_r12_neumann = rsum / a

      return
      end subroutine clmb_r12_itself

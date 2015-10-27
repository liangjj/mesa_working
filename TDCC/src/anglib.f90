module anglib
  USE nrtype
  ! Library of angular momentum coupling coefficient routines in fortran 90
  ! Paul Stevenson, Oxford University/Oak Ridge National Laboratory.
  ! spaul@mail.phy.ornl.gov
contains

!******************************************************
  function cleb(j1,m1,j2,m2,j,m)
    implicit none
! calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
! arguments are integer and twice the true value.

    REAL(DP)    :: cleb,factor,sum
    INTEGER(I4B) :: j1,m1,j2,m2,j,m,par,z,zmin,zmax

! some checks for validity (let's just return zero for bogus arguments)

    if (2*(j1/2)-int(2*(j1/2.0)) /= 2*(abs(m1)/2)-int(2*(abs(m1)/2.0)) .or.&
         2*(j2/2)-int(2*(j2/2.0)) /= 2*(abs(m2)/2)-int(2*(abs(m2)/2.0)) .or.&
         2*(j/2)-int(2*(j/2.0)) /= 2*(abs(m)/2)-int(2*(abs(m)/2.0)) .or. &
         j1<0 .or. j2<0 .or. j<0 .or. abs(m1)>j1 .or. abs(m2)>j2 .or.&
         abs(m)>j .or. j1+j2<j .or. abs(j1-j2)>j .or. m1+m2/=m) then
       cleb= 0.0
    else

       factor = 0.0
       factor = binom(j1,(j1+j2-j)/2) / binom((j1+j2+j+2)/2,(j1+j2-j)/2)
       factor = factor * binom(j2,(j1+j2-j)/2) / binom(j1,(j1-m1)/2)
       factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2)
       factor = sqrt(factor)

       zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
       zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)

       sum=0.0
       do z = zmin,zmax
          par=1
          if(2*(z/2)-int(2*(z/2.0)) /= 0) par=-1
          sum=sum+par*binom((j1+j2-j)/2,z)*binom((j1-j2+j)/2,(j1-m1)/2-z)*&
               binom((-j1+j2+j)/2,(j2+m2)/2-z)
       end do

       cleb = factor*sum
    end if
  end function cleb

!***********************************************************
!  Wigner3J000 from kenichi ishikawas code (in wigner.f90)
!***********************************************************
  FUNCTION Wigner3J000(l1,l2,l3)
    ! Wigner 3J symbol, in which all m's are zero.
    ! Wigner3J000(l1,l2,l3) = Wigner3J(l1,l2,l3,0,0,0)
    IMPLICIT none
    REAL(DP) :: Wigner3J000
    INTEGER(I4B), INTENT(IN) :: l1,l2,l3
    INTEGER(I4B) :: J
    if(mod(l1+l2+l3,2)==1)then
       Wigner3J000 = 0.d0
    else   
       J = (l1+l2+l3)/2
       if(J-l1<0.or.J-l2<0.or.J-l3<0)then
          Wigner3J000 = 0.d0
       else   
          Wigner3J000 = (-1)**J * &
               sqrt(factorial(2*J-2*l1) * factorial(2*J-2*l2) * factorial(2*J-2*l3) / factorial(2*J+1)) * &
               factorial(J) / ( factorial(J-l1) * factorial(J-l2) * factorial(J-l3) )
       endif
    endif
  END FUNCTION Wigner3J000

!******************************************************
  function sixj(a,b,c,d,e,f)
    implicit none
    INTEGER(I4B), intent(in) :: a,b,c,d,e,f
    REAL(DP) :: sixj
    INTEGER(I4B) :: phase, nlo, nhi, n
    REAL(DP) :: outfactors, sum, sumterm
! calculates a Wigner 6-j symbol. Argument a-f are integer and are
! twice the true value of the 6-j's arguments, in the form
! { a b c }
! { d e f }
! Calculated using binomial coefficients to allow for (reasonably) high
! arguments.

! First check for consistency of arguments:
    sixj=0.0
    if(mod(a+b,2)/=mod(c,2)) return
    if(mod(c+d,2)/=mod(e,2)) return
    if(mod(a+e,2)/=mod(f,2)) return
    if(mod(b+d,2)/=mod(f,2)) return
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(c-d)>e .or. c+d<e) return
    if(abs(a-e)>f .or. a+e<f) return
    if(abs(b-d)>f .or. b+d<f) return

    phase=(-1)**((a+c+d+f)/2)

    outfactors = angdelta(a,e,f)/angdelta(a,b,c)
    outfactors = outfactors * angdelta(b,d,f)*angdelta(c,d,e)

!    write(6,*) outfactors

    nlo = max( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
    nhi = min( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)

    sum=0.0
    do n=nlo,nhi
       sumterm = (-1)**n
       sumterm = sumterm * binom(n+1,n-(a+b+c)/2)
       sumterm = sumterm * binom((a+b-c)/2,n-(c+d+e)/2)
       sumterm = sumterm * binom((a-b+c)/2,n-(b+d+f)/2)
       sumterm = sumterm * binom((b-a+c)/2,n-(a+e+f)/2)
!       write(6,*) ',sumterm: ',sumterm
       sum=sum+sumterm
    end do

    sixj = phase * sum * outfactors

  end function sixj


!******************************************************
  function angdelta(a,b,c)
    implicit none
    INTEGER(I4B) :: a,b,c
    REAL(DP)    :: angdelta, scr1
! calculate the function delta as defined in varshalovich et al. for
! use in 6-j symbol:
    scr1= factorial((a+b-c)/2)
    scr1=scr1/factorial((a+b+c)/2+1)
    scr1=scr1*factorial((a-b+c)/2)
    scr1=scr1*factorial((-a+b+c)/2)
    angdelta=sqrt(scr1)
  end function angdelta


!******************************************************
  function ninej(a,b,c,d,e,f,g,h,i)
    implicit none
    INTEGER(I4B) :: a,b,c,d,e,f,g,h,i
    REAL(DP)    :: ninej, sum
    INTEGER(I4B) :: xlo, xhi
    INTEGER(I4B) :: x
! calculate a 9-j symbol. The arguments are given as integers twice the
! value of the true arguments in the form
! { a b c }
! { d e f }
! { g h i }

    ninej=0.0
! first check for bogus arguments (and return zero if so)
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(d-e)>f .or. d+e<f) return
    if(abs(g-h)>i .or. g+h<i) return
    if(abs(a-d)>g .or. a+d<g) return
    if(abs(b-e)>h .or. b+e<h) return
    if(abs(c-f)>i .or. c+f<i) return

    xlo = max(abs(b-f),abs(a-i),abs(h-d))
    xhi = min(b+f,a+i,h+d)

    sum=0.0
    do x=xlo,xhi,2
       sum=sum+(-1)**x*(x+1)*sixj(a,b,c,f,i,x)*sixj(d,e,f,b,x,h)*&
            sixj(g,h,i,x,a,d)
    end do
    ninej=sum

  end function ninej


!******************************************************
  recursive function factorial(n) result(res)
    implicit none
    INTEGER(I4B) :: n
    REAL(DP) :: res

    if (n==0 .or. n==1) then
       res=1.d0
    else
       res=n*factorial(n-1)
    end if
  end function factorial


!******************************************************
  recursive function binom(n,r) result(res)
    implicit none
    INTEGER(I4B) :: n,r
    REAL(DP) :: res
 
    if(n==r .or. r==0) then
       res = 1.0
    else if (r==1) then
       res = real(n,DP)
    else
       res = real(n,DP)/real(n-r,DP)*binom(n-1,r)
    end if
  end function binom

end module anglib

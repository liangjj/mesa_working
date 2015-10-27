!======================================================================
    FUNCTION mk2(i1,j1,i2,j2,k)
!======================================================================
!              k                      k
!   Returns ( N (i1, j1; i2, j2)  +  N (j1, i1; j2, i2) ) / 2   
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_momentc
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, iv, ii, ik, ivi,ivj
    REAL(KIND=8), DIMENSION(nv) :: v1, v2, v3, v4
    REAL(KIND=8), DIMENSION(ks*(ks+1)/2) :: a,b
    REAL(KIND=8) :: mk2, s1, s2

    !    moments calculations

    Call mmk2(k)

    ik= ks*(ks+1)/2

    mk2 = 0.d0

    Do iv = 1,nv

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=i,ks
         ivj=iv+j-1
         ii = ii+1
         if(i.eq.j) then
          a(ii) = p(ivi,i1)*p(ivj,i2)
          b(ii) = p(ivi,j1)*p(ivj,j2)
         else
          a(ii) = p(ivi,i1)*p(ivj,i2) + p(ivi,i2)*p(ivj,i1)
          b(ii) = p(ivi,j1)*p(ivj,j2) + p(ivi,j2)*p(ivj,j1)
         end if
        End do
       End do

       v1(iv) = SUM(a*rkd2(1:ik,iv))
       v2(iv) = SUM(b*rkd1(1:ik,iv))
       v3(iv) = SUM(a*rkd1(1:ik,iv))
       v4(iv) = SUM(b*rkd2(1:ik,iv))

       ! .. diagonal contributions

       Do i = 1,ik
        Do j = 1,ik
         mk2 = mk2 + a(i)*rkd(i,j,iv)*b(j)
        End do
       End do

    End do

    ! the upper and lower regions

    s1 = 0.d0
    s2 = 0.d0
    Do iv = 2,nv
      s1 = s1 + v2(iv-1)
      mk2 = mk2 + s1*v1(iv)
      s2 = s2 + v3(iv-1)
      mk2 = mk2 + s2*v4(iv)
    End do

    mk2 = mk2/2 * fine
    
    End function mk2


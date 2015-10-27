!======================================================================
    REAL(8) FUNCTION tk2 (i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  T (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_momentc
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, iv, ii, ivi,ivj
    REAL(KIND=8), DIMENSION(nv) :: v1, v2, v3, v4
    REAL(KIND=8), DIMENSION(ks*ks) :: a,b
    REAL(KIND=8) :: s1, s2

    !    moments calculations

    Call Mtk2(k)

    tk2 = 0.d0

    Do iv = 1,nv

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=1,ks
         ivj=iv+j-1
         ii = ii+1
         a(ii) = p(ivi,i1)*p(ivj,i2)
         b(ii) = p(ivi,j1)*p(ivj,j2)
        End do
       End do

       v1 (iv) = SUM(a*rkd1(:,iv))
       v2 (iv) = SUM(b*rkd1(:,iv))
       v3 (iv) = SUM(a*rkd2(:,iv))
       v4 (iv) = SUM(b*rkd2(:,iv))

       ! .. diagonal contributions

       Do i = 1,ks*ks
        Do j = 1,ks*ks
          tk2 = tk2 + a(i)*rkd(i,j,iv)*b(j) 
        End do
       End do

    End do

    ! the upper and lower regions

    s1 = 0.d0
    s2 = 0.d0
    Do iv =  2,nv
      s1 = s1 + v1(iv-1)
      tk2 = tk2 + s1*v4(iv)
      s2 = s2 + v2(iv-1)
      tk2 = tk2 + s2*v3(iv)
    End do

    tk2 = tk2 * fine / (k + k + 1)

    End FUNCTION tk2


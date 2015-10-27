!======================================================================
    FUNCTION vk2(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  V (i1, j1; i2, j2) using  the  cell-summation algorith
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_momentc
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER(4) :: i,j, iv, ivi, ivj, ii, jk 
    REAL(8), DIMENSION(nv) :: v1, v2, v3, v4
    REAL(8), DIMENSION(ks*ks) ::  a, b
    REAL(8) :: vk2
    REAL(8) :: s1, s2

    ! .. check the need of calculations
                
    Call Mvk2(k)

    jk = ks*ks

    vk2 = 0.d0

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

     v1(iv) = SUM(rkd3(:,iv)*a)
     v2(iv) = SUM(rkd2(:,iv)*b)
     v3(iv) = SUM(rkd4(:,iv)*a)
     v4(iv) = SUM(rkd1(:,iv)*b)
     
     ! the diagonal cell contribution
       
     Do i = 1,jk
      Do j = 1,jk
       vk2 = vk2 + a(i)*rkd(i,j,iv)*b(j)
      End do
     End do
     
    End do 
     
    ! the upper and lower regions

    s1 = 0.d0
    s2 = 0.d0
    Do iv =  2,nv
      s1 = s1 + v1(iv-1)
      vk2 = vk2 + s1*v2(iv)
      s2 = s2 + v4(iv-1)
      vk2 = vk2 + s2*v3(iv)
    End do     
     
    vk2 = vk2 *fine

    EndFUNCTION vk2


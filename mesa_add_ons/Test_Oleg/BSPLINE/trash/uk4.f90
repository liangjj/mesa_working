!======================================================================
    FUNCTION uk4(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  U (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_moments
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, iv, ii, ik,jk, ivi,ivj
    REAL(KIND=8), DIMENSION(nv) :: v1, v2, v3, v4
    REAL(KIND=8), DIMENSION(ks*ks) :: a,b,c,d
    REAL(KIND=8) :: uk4, s1, s2, u1,u2

    !    moments calculations

    Call muk4(k)

    ik = ks*(ks+1)/2
    jk = ks * ks
    
    uk4 = 0.d0

    Do iv = 1,nv

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=i,ks
         ivj=iv+j-1
         ii = ii+1
         if(i.eq.j) then
          a(ii) = p(ivi,i1)*p(ivj,i2)
         else
          a(ii) = p(ivi,i1)*p(ivj,i2) + p(ivi,i2)*p(ivj,i1)
         end if
        End do
       End do

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=1,ks
         ivj=iv+j-1
         ii = ii+1
         b(ii) = p(ivi,j1)*p(ivj,j2)
        End do
       End do

       c(1:ik) = a(1:ik)*rkm1(1:ik,iv)
       v1(iv) = SUM(c(1:ik))
       c = b*rkm2(:,iv)
       v2(iv) = SUM(c)
       c = b*rkm3(:,iv)
       v3(iv) = SUM(c)
       c(1:ik) = a(1:ik)*rkm4(1:ik,iv)
       v4(iv) = SUM(c(1:ik))

       ! .. diagonal contributions

       Do j = 1,jk
         c = rkd(1:ik,j,iv)*a(1:ik)
         d(j) =  SUM(c)
       End do
       c = d*b
       uk4 = uk4 + SUM(c)

    End do

    ! the lower and  upper regions

    s1 = 0.d0
    s2 = 0.d0
    u1 = 0.d0
    u2 = 0.d0
    
    Do iv =  2,nv
      s1 = s1 + v2(iv-1)
      u1 = u1 + s1*v1(iv)
      s2 = s2 + v4(iv-1)
      u2 = u2 + s2*v3(iv)
    End do

    write(*,'(i3,3E12.3)') k,uk4, u1,u2
    
    uk4 = uk4 + (k-1)*u1 - (k+2)*u2

    uk4 = uk4 * fine  / (2*k+1)

    End FUNCTION uk4


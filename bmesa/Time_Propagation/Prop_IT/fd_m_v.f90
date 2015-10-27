!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE fd_m_v
                         USE dvrprop_global_it
                         USE dvrprop_global_rt
                         USE dvr_shared
                         USE dvr_global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck ds3_bmm_d
!**begin prologue     ds3_bmm_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        the two non-zero elements of the symmetric
!**                   tridiagonal matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmm_d
  SUBROUTINE ds3_bmm_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n4-1
      
!\begin{eqnarray}
! v_{out}(i,j,k) &=&  v_{out}(i,j,k) + m(i,i+1) v_{in}(i+1,j,k)
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!            &=&  v_{out}(i,j,k)     + m(i+1,i) v_{in}(i+1,j,k) \nonumber
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                            &
                                   +                             &
                       grid(dim)%ke(2,i) * v_in(i+1,:,:,:)       &
                                   +                             &
                       grid(dim)%ke(1,i) * v_in(i,:,:,:)         &
                                   +                             &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)
  END DO

!\center Do the special case of the first and last row.

      v_out(1,:,:,:)  = v_out(1,:,:,:)                           & 
                                   +                             &
                        grid(dim)%ke(2,1) * v_in(2,:,:,:)        &
                                   +                             &
                        grid(dim)%ke(1,1) * v_in(1,:,:,:)
      v_out(n4,:,:,:) = v_out(n4,:,:,:)                          &
                                   +                             &
                        grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)  & 
                                   +                             &
                        grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)
END SUBROUTINE ds3_bmm_d
!deck ds3_bmm_z
!**begin prologue     ds3_bmm_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        the two non-zero elements of the symmetric
!**                   tridiagonal matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmm_z
  SUBROUTINE ds3_bmm_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n4-1
      
!\begin{eqnarray}
! v_{out}(i,j,k) &=&  v_{out}(i,j,k) + m(i,i+1) v_{in}(i+1,j,k)
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!            &=&  v_{out}(i,j,k)     + m(i+1,i) v_{in}(i+1,j,k) \nonumber
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                            &
                                   +                             &
                       grid(dim)%ke(2,i) * v_in(i+1,:,:,:)       &
                                   +                             &
                       grid(dim)%ke(1,i) * v_in(i,:,:,:)         &
                                   +                             &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)
  END DO

!\center Do the special case of the first and last row.

      v_out(1,:,:,:)  = v_out(1,:,:,:)                           & 
                                   +                             &
                        grid(dim)%ke(2,1) * v_in(2,:,:,:)        &
                                   +                             &
                        grid(dim)%ke(1,1) * v_in(1,:,:,:)
      v_out(n4,:,:,:) = v_out(n4,:,:,:)                          &
                                   +                             &
                        grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)  & 
                                   +                             &
                        grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)
END SUBROUTINE ds3_bmm_z
!deck ds3_bmmt_d.f
!**begin prologue     ds3_bmmt_d
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        the two non-zero elements of the tridiagonal matrix
!**                   are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmmt_d
  SUBROUTINE ds3_bmmt_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n3-1
      
!\begin{eqnarray}
! v_{out}(j,i,k) &=&  v_{out}(j,i,k) + m(i,i+1) v_{in}(j,i+1,k)
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!              &=&  v_{out}(j,i,k)   + m(i+1,i) v_{in}(j,i+1,k) \nonumber
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(:,i,:,:) = v_out(:,i,:,:)                              &
                                         +                         &
                       grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)       &
                                         +                         &
                       grid(dim)%ke(1,i)   * v_in(:,i,:,:)         &
                                         +                         &
                       grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special case of the first and last row.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                 &
                                         +                         &
                    grid(dim)%ke(2,1)    * v_in(:,2,:,:)           &
                                         +                         &
                    grid(dim)%ke(1,1)    * v_in(:,1,:,:)
  v_out(:,n3,:,:) = v_out(:,n3,:,:)                                &
                                         +                         &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)        &
                                         +                         &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)
END SUBROUTINE ds3_bmmt_d
!deck ds3_bmmt_z.f
!**begin prologue     ds3_bmmt_z
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        the two non-zero elements of the tridiagonal matrix
!**                   are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmmt_z
  SUBROUTINE ds3_bmmt_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n3-1
      
!\begin{eqnarray}
! v_{out}(j,i,k) &=&  v_{out}(j,i,k) + m(i,i+1) v_{in}(j,i+1,k)
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!              &=&  v_{out}(j,i,k)   + m(i+1,i) v_{in}(j,i+1,k) \nonumber
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(:,i,:,:) = v_out(:,i,:,:)                              &
                                         +                         &
                       grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)       &
                                         +                         &
                       grid(dim)%ke(1,i)   * v_in(:,i,:,:)         &
                                         +                         &
                       grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special case of the first and last row.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                 &
                                         +                         &
                    grid(dim)%ke(2,1)    * v_in(:,2,:,:)           &
                                         +                         &
                    grid(dim)%ke(1,1)    * v_in(:,1,:,:)
  v_out(:,n3,:,:) = v_out(:,n3,:,:)                                &
                                         +                         &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)        &
                                         +                         &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)
END SUBROUTINE ds3_bmmt_z
!deck ds5_bmm_d
!**begin prologue     ds5_bmm_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmm_d
  SUBROUTINE ds5_bmm_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                            :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                            :: i
!
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO  i=3,n4-2
!\begin{eqnarray}
! v_{out}(i,j,k) =  v_{out}(i,j,k) + m(i,i+2) v_{in}(i+2,j,k)
!                                  + m(i,i+1) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!                                  + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!                =  v_{out}(i,j,k) + m(i+2,i) v_{in}(i+2,j,k)
!                                  + m(i+1,i) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(i-1,j,k)
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)
  END DO
!
!\center Do the special cases of the first two and last two rows.
!
  v_out(1,:,:,:)   =  v_out(1,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(1,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,1) * v_in(3,:,:,:)
  v_out(2,:,:,:)   =  v_out(2,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(1,2) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,2) * v_in(3,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,2) * v_in(4,:,:,:)
  v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                 &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)           &
                                      +                                 & 
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)           &
                                      +                                 &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)
  v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                   &
                                      +                                 &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)           & 
                                      +                                 &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)
END SUBROUTINE ds5_bmm_d
!deck ds5_bmm_z
!**begin prologue     ds5_bmm_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmm_z
  SUBROUTINE ds5_bmm_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO  i=3,n4-2
!\begin{eqnarray}
! v_{out}(i,j,k) =  v_{out}(i,j,k) + m(i,i+2) v_{in}(i+2,j,k)
!                                  + m(i,i+1) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!                                  + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!                =  v_{out}(i,j,k) + m(i+2,i) v_{in}(i+2,j,k)
!                                  + m(i+1,i) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(i-1,j,k)
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)
  END DO
!
!\center Do the special cases of the first two and last two rows.
!
  v_out(1,:,:,:)   =  v_out(1,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(1,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,1) * v_in(3,:,:,:)
  v_out(2,:,:,:)   =  v_out(2,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(1,2) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,2) * v_in(3,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,2) * v_in(4,:,:,:)
  v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                 &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)           &
                                      +                                 & 
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)           &
                                      +                                 &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)
  v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                   &
                                      +                                 &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)           & 
                                      +                                 &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)
END SUBROUTINE ds5_bmm_z
!deck ds5_bmmt_d.f
!**begin prologue     ds5_bmmt_d
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmmt_d
  SUBROUTINE ds5_bmmt_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8,     DIMENSION(n4,n3,n2,n1)        :: v_in, v_out
  INTEGER                                :: i
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO i=3,n3-2
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                              &  
                                  +                               &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)           &
                                  +                               &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)         &
                                  +                               &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)
  END DO

!\center Do the special cases of the first two and last two rows.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                &
                                  +                               &
                  grid(dim)%ke(1,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,1) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,1) * v_in(:,3,:,:)
  v_out(:,2,:,:)  = v_out(:,2,:,:)                                &
                                  +                               & 
                  grid(dim)%ke(2,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(1,2) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,2) * v_in(:,3,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,2) * v_in(:,4,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                           &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                             &
                                  +                               &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)
END SUBROUTINE ds5_bmmt_d
!deck ds5_bmmt_z.f
!**begin prologue     ds5_bmmt_z
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmmt_z
  SUBROUTINE ds5_bmmt_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                 :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)      :: v_in, v_out
  INTEGER                                 :: i
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO i=3,n3-2
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                              &  
                                  +                               &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)           &
                                  +                               &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)         &
                                  +                               &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)
  END DO

!\center Do the special cases of the first two and last two rows.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                &
                                  +                               &
                  grid(dim)%ke(1,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,1) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,1) * v_in(:,3,:,:)
  v_out(:,2,:,:)  = v_out(:,2,:,:)                                &
                                  +                               & 
                  grid(dim)%ke(2,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(1,2) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,2) * v_in(:,3,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,2) * v_in(:,4,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                           &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                             &
                                  +                               &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)
END SUBROUTINE ds5_bmmt_z
!deck ds7_bmm_d.f
!**begin prologue     ds7_bmm_d
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%(4,i) = m(i,i+3)
!\end{eqnarray}

!**references
!**routines called
!**end prologue       ds7_bmm_d
  SUBROUTINE ds7_bmm_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix
!
  DO  i=4,n4-3
      
!\begin{eqnarray}
! v_{out}(i,:,:) =  v_{out}(i,:,:) + m(i,i+3) v_{in}(i+3,:,:)
!                                  + m(i,i+2) v_{in}(i+2,:,:)
!                                  + m(i,i+1) v_{in}(i+1,:,:)  \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:) \nonumber \\
!                =  v_{out}(i,:,:) + m(i+3,i) v_{in}(i+3,:,:)
!                                  + m(i+2,i) v_{in}(i+2,:,:)
!                                  + m(i+1,i) v_{in}(i+1,:,:) \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)

      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(4,i)   * v_in(i+3,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)             &
                                      +                                  &   
                       grid(dim)%ke(4,i-3) * v_in(i-3,:,:,:)               
  END DO

!\center Do the special cases of the first three and last three rows.

    v_out(1,:,:,:) = v_out(1,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(1,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,1) * v_in(4,:,:,:)
    v_out(2,:,:,:) = v_out(2,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(1,2) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,2) * v_in(4,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(4,2) * v_in(5,:,:,:)
    v_out(3,:,:,:) = v_out(3,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(2,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(1,3) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,3) * v_in(4,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,3) * v_in(5,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,3) * v_in(6,:,:,:)
    v_out(n4-2,:,:,:) = v_out(n4-2,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(1,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-3) * v_in(n4-3,:,:,:)            & 
                                      +                                  &
                      grid(dim)%ke(3,n4-4) * v_in(n4-4,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-5) * v_in(n4-5,:,:,:)
    v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-4) * v_in(n4-4,:,:,:)
    v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                  &
                                      +                                  &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-3) * v_in(n4-3,:,:,:)
END SUBROUTINE ds7_bmm_d
!!deck ds7_bmm_z.f
!**begin prologue     ds7_bmm_z
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%(4,i) = m(i,i+3)
!\end{eqnarray}

!**references
!**routines called
!**end prologue       ds7_bmm_z
  SUBROUTINE ds7_bmm_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix
!
  DO  i=4,n4-3
      
!\begin{eqnarray}
! v_{out}(i,:,:) =  v_{out}(i,:,:) + m(i,i+3) v_{in}(i+3,:,:)
!                                  + m(i,i+2) v_{in}(i+2,:,:)
!                                  + m(i,i+1) v_{in}(i+1,:,:)  \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:) \nonumber \\
!                =  v_{out}(i,:,:) + m(i+3,i) v_{in}(i+3,:,:)
!                                  + m(i+2,i) v_{in}(i+2,:,:)
!                                  + m(i+1,i) v_{in}(i+1,:,:) \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)

      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(4,i)   * v_in(i+3,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)             &
                                      +                                  &   
                       grid(dim)%ke(4,i-3) * v_in(i-3,:,:,:)               
  END DO

!\center Do the special cases of the first three and last three rows.

    v_out(1,:,:,:) = v_out(1,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(1,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,1) * v_in(4,:,:,:)
    v_out(2,:,:,:) = v_out(2,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(1,2) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,2) * v_in(4,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(4,2) * v_in(5,:,:,:)
    v_out(3,:,:,:) = v_out(3,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(2,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(1,3) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,3) * v_in(4,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,3) * v_in(5,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,3) * v_in(6,:,:,:)
    v_out(n4-2,:,:,:) = v_out(n4-2,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(1,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-3) * v_in(n4-3,:,:,:)            & 
                                      +                                  &
                      grid(dim)%ke(3,n4-4) * v_in(n4-4,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-5) * v_in(n4-5,:,:,:)
    v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-4) * v_in(n4-4,:,:,:)
    v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                  &
                                      +                                  &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-3) * v_in(n4-3,:,:,:)
END SUBROUTINE ds7_bmm_z
!deck ds7_bmmt_d.f
!**begin prologue     ds7_bmmt_d
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%ke(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%ke(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%ke(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%ke(4,i) = m(i,i+3)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds7_bmmt_d
  SUBROUTINE ds7_bmmt_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix

  DO i=4,n3-3
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+3) v_{in}(j,i+3,k)
!                                  + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+3,i) v_{in}(j,i+3,k)
!                                  + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                                &
                                      +                             &
                    grid(dim)%ke(4,i)   * v_in(:,i+3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)           & 
                                      +                             &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)           &
                                      +                             &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)             &
                                      +                             &
                    grid(dim)%ke(4,i-3) * v_in(:,i-3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)           &
                                      +                             &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special cases of the first three and last three rows.

  v_out(:,1,:,:) = v_out(:,1,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(1,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,1) * v_in(:,4,:,:)
  v_out(:,2,:,:) = v_out(:,2,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,1,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(1,2) * v_in(:,2,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,2) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,2) * v_in(:,5,:,:)
  v_out(:,3,:,:) = v_out(:,3,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(1,3) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,3) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,3) * v_in(:,5,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,3) * v_in(:,6,:,:)
  v_out(:,n3-2,:,:) = v_out(:,n3-2,:,:)                             &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3,:,:)           &
                                   +                                &
                   grid(dim)%ke(2,n3-2) * v_in(:,n3-1,:,:)          &
                                   +                                &
                   grid(dim)%ke(1,n3-2) * v_in(:,n3-2,:,:)          &
                                   +                                &
                   grid(dim)%ke(2,n3-3) * v_in(:,n3-3,:,:)          &
                                   +                                &
                   grid(dim)%ke(3,n3-4) * v_in(:,n3-4,:,:)          &
                                   +                                &
                   grid(dim)%ke(4,n3-5) * v_in(:,n3-5,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                             &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-4) * v_in(:,n3-4,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                               &
                                   +                                &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-3) * v_in(:,n3-3,:,:)
END SUBROUTINE ds7_bmmt_d
!deck ds7_bmmt_z.f
!**begin prologue     ds7_bmmt_z
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%ke(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%ke(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%ke(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%ke(4,i) = m(i,i+3)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds7_bmmt_z
  SUBROUTINE ds7_bmmt_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                  :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)       :: v_in, v_out
  INTEGER                                  :: i
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix

  DO i=4,n3-3
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+3) v_{in}(j,i+3,k)
!                                  + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+3,i) v_{in}(j,i+3,k)
!                                  + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                                &
                                      +                             &
                    grid(dim)%ke(4,i)   * v_in(:,i+3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)           & 
                                      +                             &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)           &
                                      +                             &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)             &
                                      +                             &
                    grid(dim)%ke(4,i-3) * v_in(:,i-3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)           &
                                      +                             &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special cases of the first three and last three rows.

  v_out(:,1,:,:) = v_out(:,1,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(1,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,1) * v_in(:,4,:,:)
  v_out(:,2,:,:) = v_out(:,2,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,1,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(1,2) * v_in(:,2,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,2) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,2) * v_in(:,5,:,:)
  v_out(:,3,:,:) = v_out(:,3,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(1,3) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,3) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,3) * v_in(:,5,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,3) * v_in(:,6,:,:)
  v_out(:,n3-2,:,:) = v_out(:,n3-2,:,:)                             &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3,:,:)           &
                                   +                                &
                   grid(dim)%ke(2,n3-2) * v_in(:,n3-1,:,:)          &
                                   +                                &
                   grid(dim)%ke(1,n3-2) * v_in(:,n3-2,:,:)          &
                                   +                                &
                   grid(dim)%ke(2,n3-3) * v_in(:,n3-3,:,:)          &
                                   +                                &
                   grid(dim)%ke(3,n3-4) * v_in(:,n3-4,:,:)          &
                                   +                                &
                   grid(dim)%ke(4,n3-5) * v_in(:,n3-5,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                             &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-4) * v_in(:,n3-4,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                               &
                                   +                                &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-3) * v_in(:,n3-3,:,:)
END SUBROUTINE ds7_bmmt_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          END MODULE fd_m_v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

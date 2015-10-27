*deck cua3d.f
c***begin prologue     cua3d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           unitary transformation, 3d
c***author             schneider, barry (nsf)
c***source             trap3d              
c***purpose            davidson preconditioner 
c***description        perform the transformation on a vector
c***                   from the dvr to the representation defined by a
c***                   separable zeroth order hamiltonian.  this is used
c***                   as a preconditioner in the davidson routine.                  
c***references         
c
c                      the arrays resid and vec may be the same but must
c                      be distinct from t1 and t2.
c                      note that the vector and residual are assumed to
c                      be calculated in packed form and as if they were
c                      matrices a(n3,n2,n1).
c                      
c***routines called    
c***end prologue       cua3d

      subroutine cua3d(vec,resid,u01,u02,u03,t1,t2,n1,n2,n3,m)
      implicit integer (a-z)
      complex*16 resid, vec, u01, u02, u03, t1, t2
      character*80 title
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension vec(n3,n2,n1,m), resid(n3,n2,n1,m), t1(n3,*)
      dimension t2(n3,n2,n1)
      common/io/inp, iout
c
c         transform residual from dvr to h0 representation
c 
      do 10 i=1,m
         do 20 j=1,n1
c
c           transform on the n3 and n2 variables for fixed n1
c         
            call cebc(t1,u03,vec(1,1,j,i),n3,n3,n2)
            call cebct(t2(1,1,j),t1,u02,n3,n2,n2)
 20      continue
c
c        finish the transformation on n1 treating the vector as an
c        (n3*n2,n1) array.
c
         call cebct(resid(1,1,1,i),t2,u01,n3*n2,n1,n1)          
 10   continue   
      return
      end       




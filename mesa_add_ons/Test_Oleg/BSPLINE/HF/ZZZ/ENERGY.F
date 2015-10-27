!=======================================================================
   SUBROUTINE energy(etotal) 
!=======================================================================
!
!       determines the total energy
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: energy

      PARAMETER(nt=100,kx=8,nwd=10)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
      DOUBLE PRECISION ek(nwd)
*
*  *****  compute kinetic energy if necessary
*
      do i = 1,nwf
        ek(i) = -0.5d0*hl(el,i,i,rel)
      end do
*
      etotal = 0.d0
      do i = 1,nwf
         etotal = etotal + sum(i)*(ek(i))
	 print *, 'etotal', ek(i), etotal
         do j = 1,i
            do k = 0,2*min0(l(i),l(j)),2
               c = a(i,j,k)*sum(i)
               if (i .eq. j) c = c/2.d0
               if (abs(c).ne.0.d0) etotal = etotal + c*fky(i,j,k,rel)
	       print *, i,j,c,fky(i,j,k,rel),etotal
            end do
         end do
         do j = 1,i-1
           do k = abs(l(i)-l(j)),l(i)+l(j),2
             c = b(i,j,k)*sum(i)
             if (abs(c).ne.0.d0) etotal=etotal+c*gky(i,j,k,rel)
	       print *, i,j,c,gky(i,j,k,rel),etotal
           end do
         end do
       end do
     END SUBROUTINE energy

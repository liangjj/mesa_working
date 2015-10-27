

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  1D RSPF Program: Real time(complex) Propagation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 MODULE gen_data
!
!   Set the dimensions for the general arrays throughout
!   the program
!
! np8 = Maximum number of radial points
! nprmx = Maximum number of potential parameters
!
   REAL*8, SAVE :: zero=0.e0, one=1.e0, two=2.e0
   REAL*8, SAVE :: convfs=2.419e-02, convh=27.2114e0
   REAL*8, SAVE :: bohr=0.529177e0, conv1=2.19473e+5
!
!  Global arrays for trapezoidal points and weights
!     and initial wavefunction
!
   REAL*8, DIMENSION(:), AllOCATABLE, SAVE :: xn,wtn,wfbnd,wfn_acc, &
           wfr,wfi,vpt,upr,upi,vij,vijt,eig_prop,upr_k,upi_k,      &
           work,tm1,tm2,tm3,tm4,prb1,pbr,pbi,pb,prb2r,prb2i,pb1
   REAL*8, DIMENSION(:,:),ALLOCATABLE, SAVE :: U_prop,U_propT


END MODULE gen_data
MODULE splint_arrays
!
!  Arrays for B-spline fit of initial wavefunction
!
!
!
   REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: rs,wftmp,dy2,u
!
END MODULE splint_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       program pdep_ms
       USE gen_data
       USE splint_arrays
!
       IMPLICIT NONE
!
       REAL*4 :: ttx1,ttx2,tarry(2),etime
       REAL*8, PARAMETER :: tol=1.e-08
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: bb,dd,ee
       REAL*8 :: x_shift,phase,phase_coeff,t_restart,&
                 alpha,scale_ho,velocity,tole,&
                 scl_de,d_0,d_1,d_2,d_3,pi,rstart,rstop,tstart,&
                 tstop,delt,delr,xxx,atoms_N,a_d0,capdl2,ddd,uuu,&
                 ust,uu1,tt2,tt4,t,ttx3,ttx4,ttt,tav,pbxx,term_nl,&
                 yp1,ypn
       INTEGER :: n_shift,mask,n_work,info,order_de,nrp,nrad,&
                  ntop,ntp,ntime,ir,it,nrp1,i,k,iwftyp,itxx,nrp5
!
       open (5,file='rs2ms.dat',status='old')
       open (6,file='rs2ms.out',status='unknown')
       open(22,file='wfn1.fle',status='unknown')
       open(23,file='wfn2.fle',status='unknown')
       open(31,file='pot.fle',status='unknown')
       open(32,file='prj.fle',status='unknown')
       open(33,file='prb.fle',status='unknown')
       open(34,file='rexp.fle',status='unknown')
       open(35,file='pbe.fle',status='unknown')
       open(50,file='phase.fle',status='unknown')
       open(60,file='wav.fle',status='unknown')
!
       pi = acos(-one)
!
!
!
!***************************************************************
! DESCRIPTION
!***************************************************************
!
! This program solves a coupled-state time dependent
! Schrodinger Equation in a single spatial variable
!
!****
!
! method of solution: real-space product formula(2)
!
!
!
!** all integrations by trapezoidal rule since the radial
!       points are equally spaced
!
!
!**********************************************************
!**********************************************************
!
!
!**
!** INPUT: radial mesh parameters
!**
!
! nrad = no. of points bebtween successive printings
!         of the solution
! nrp = no. of points in x mesh( must be EVEN)
! rstart = starting value of x
! rstop = final value of x
!
!
        read(5,*) nrad,nrp,rstart,rstop
!
!
!**
!** INPUT: temporal mesh parameters
!**
!
! ntop = input option
!        0 - read in ntp,tstart,tstop - calculate delt
!        1 - read in ntp,tstart,delt - determine tstop
! ntp = no. of temporal points
! tstart = starting value of t
! tstop = final value of t
! ntime = increment brtween time printings
!
       read(5,*) ntop,ntp,tstart,tstop,ntime,delt
!
!
!**
!** INPUT: potential parameters
!**
!
!**** Diagonal potential
!
!             atoms_N = number of atoms
!             a_d0 = nonlinear scaling term
!
         read(5,*) atoms_N,a_d0,mask,phase,phase_coeff,x_shift,n_shift
         term_nl = atoms_N*a_d0
         phase = phase*pi
!
!
!
!**
!** INPUT: additional parameters
!**
!
! iwftyp = type of initial wavefunction
!          1 - read IT solution(real part only) from unit 44
!          2 - gaussian - reall part only of width alpha and velocity
!          3 - restart RT solution (complex)
! scale_ho = scale of the HO potential for NLS
! tole = tolerence in the check of the eigenvalues of T
! order_de = order of the difference scheme
!            3 - three-point
!            5 - five-point
!
       read(5,*)  iwftyp,alpha,velocity,scale_ho,tole,order_de
!
!
!
!
!
!***************************************************
!
!
!
       write(6,100)
!
! form of temporal mesh
!
       delr = (rstop - rstart)/float(nrp)
       if(ntop .eq. 0) then
          delt = (tstop - tstart)/float(ntp)
       else if(ntop .eq. 1) then
          tstop = float(ntp)*delt
       endif
!
       write(6,101) nrp,rstart,rstop,delr
       write(6,102) ntp,tstart,tstop,delt
       write(6,122)
!       x_shift = n_shift*delr
       write(6,124) atoms_N,a_d0,phase,phase_coeff,x_shift
!
!
! set the radial mesh
!
!
       capdl2 = 1.e0/(2.e0*(delr**2))
!
       nrp1 = nrp + 1
       ddd = delt*capdl2
       term_nl = term_nl*(delt/2.e0)
!
       write(6,755) delt,delr,ddd
!
!
!! Allocate arrays
!
      nrp5 = nrp1+5
       ALLOCATE(xn(1:nrp1),wtn(1:nrp1),wfr(1:nrp1),wfi(1:nrp1),upr(1:nrp1),&
                upi(1:nrp1),vpt(1:nrp1),vij(1:nrp1),vijt(1:nrp1), &
                eig_prop(1:nrp5),upr_k(1:nrp1),upi_k(1:nrp1),work(1:5*nrp1), &
                U_prop(1:nrp1,1:nrp1),U_propT(1:nrp5,1:nrp5),rs(1:nrp1), &
                wftmp(1:nrp1),dy2(1:nrp1),wfbnd(1:nrp1),wfn_acc(1:nrp1), &
                
bb(1:nrp1,1:nrp1),dd(1:nrp1,1:nrp1),ee(1:nrp1,1:nrp1),tm1(1:nrp1),&
                tm2(1:nrp1),tm3(1:nrp1),tm4(1:nrp1),prb1(1:nrp1),pbr(1:nrp1),&
                pbi(1:nrp1),pb(1:nrp1),prb2r(1:nrp1),prb2i(1:nrp1),u(1:nrp1), &
                pb1(1:nrp1))
!
! set_up a pairwise radial array - xn
!       and the weights for spatial integration - wtn
!
       xxx = rstart - delr
       do ir = 1,nrp1
          xxx = xxx + delr
          xn(ir) = xxx
          wtn(ir) = delr
       enddo
!
!
       tt2 = delt/2.
       tt4 = delt/4.
!
! DIAGONAL, TIME-INDEPENDENT POTENTIAL and PROPAGATOR
!
        do i = 1,nrp1
           call potvdi(scale_ho,xn(i),ust)
            uuu = (ust)*tt2
            vpt(i) = uuu
        enddo
!
       upr(1:nrp1) = cos(vpt(1:nrp1))
       upi(1:nrp1) = sin(vpt(1:nrp1))
!
!
! WRITE: POTENTIALs to unit 31
!
          write(31,1331)
          do i = 1,nrp1
             uu1 = (vpt(i)/tt4)
             write(31,211) xn(i),uu1
          enddo
!
!          write(31,1331)
!          do i = 1,nrp1
!               uu1 = vij(i)/tt4
!               write(31,211) xn(i),uu1
!          enddo
!
!
!
!
! set the initial conditions at the first time point t0
!
!
       write(23,106) nrp1
       call funint(nrp1,mask,phase,phase_coeff,x_shift,n_shift,   &
                   iwftyp,t_restart,alpha,velocity)
!
!
!
!
! KINETIC ENERGY OPERATOR: PARAMETERS B
!
!
!
!
!
       write(32,1331)
       write(33,1331)
       write(34,1331)
!
!
! propagation in time
!
!
       if(iwftyp <= 2) then
          t = tstart
       elseif(iwftyp == 3) then
          t = t_restart
       endif
       ttx1 = etime(tarry)
!
! Set the diagonal KE propagator
!
       U_propT = 0.e0
       if(order_de == 3) then
          scl_de = -1.e0/(2.e0*delr**2)
          d_0 = -2.e0*scl_de
          d_1 = scl_de
          d_2 = 0.e0
          d_3 = 0.e0
       elseif(order_de == 5) then
          scl_de = -1.e0/(24.e0*delr**2)
          d_3 = 0.e0
          d_2 = -1.e0*scl_de
          d_1 = +16.e0*scl_de
          d_0 = -30.e0*scl_de
       elseif(order_de == 7) then
          scl_de = -1.e0/(720.e0*delr**2)
          d_3 = 4.e0*scl_de
          d_2 = -54.e0*scl_de
          d_1 = +540.e0*scl_de
          d_0 = -980.e0*scl_de
       endif
       do i = 1,nrp1
          U_propT(i,i) = d_0
          U_propT(i,i+1) = d_1
          U_propT(i+1,i) = d_1
          U_propT(i,i+2) = d_2
          U_propT(i+2,i) = d_2
          U_propT(i,i+3) = d_3
          U_propT(i+3,i) = d_3
       enddo
       ee = U_propT
       n_work = 5*nrp1
       call dsyev('v','l',nrp1,U_propT,nrp5,eig_prop,work,n_work,info)
       U_prop = U_propT
       vijt(1:nrp1) = eig_prop(1:nrp1)*delt
       upr_k(1:nrp1) = cos(vijt(1:nrp1))
       upi_k(1:nrp1) = sin(vijt(1:nrp1))
       do i = 1,nrp1
          do k = 1,nrp1
             U_propT(i,k) = U_prop(k,i)
          enddo
       enddo
       bb = 0.e0
       do i = 1,nrp1
          bb(i,i) = eig_prop(i)
       enddo
       dd(1:nrp1,1:nrp1) = matmul(ee(1:nrp1,1:nrp1),U_prop(1:nrp1,1:nrp1))
       bb(1:nrp1,1:nrp1) = matmul(U_propT(1:nrp1,1:nrp1),dd(1:nrp1,1:nrp1))
       write(6,212)
       do i = 1,nrp1
          if(abs(bb(i,i) - eig_prop(i)) > tole) then
             write(6,211) bb(i,i),eig_prop(i)
          endif
       enddo
       ttx2 = etime(tarry) - ttx1
!
!********************************************
!
! TEMPORAL PROPAGATION
!
!********************************************
!
       write(6,105)
       ttx3 = etime(tarry)
       do it = 1,ntp
          t = t + delt
!
! add in time-dependent and/or nonlinear  potential
!
          call timfn(t,delt,nrp1,term_nl)
!
          call expu(nrp1)
          call expk(nrp1)
          call expu(nrp1)
!
!
!
          if(mod(it,ntime) .eq. 0) then
              call prbfn(it,nrp1,iwftyp,t,delr,pbxx)
          endif
!
!
       enddo
!
       ttx4 = etime(tarry) - ttx3
       ttt = ttx4
       tav = ttt/ntp
       write(6,811) ntp,ttt,tav,ttx2
!
!
 100   format(//' pde solution by product formula'///)
 101   format(// &
     ' radial mesh'//' no. of pomits=',i10/ &
     ' rstart=',e15.6,5x,' rstop=',e15.6/' delr=',e15.6)
 102   format(//' temporal mesh'//' no. of ponts=',i10/ &
     ' tstart=',f10.5,' tstop=',f15.5/' delt=',e15.6/)
 105   format(///' propagation in r and t'//)
 106   format(i10)
 122   format(//' Potential Type'//)
 123   format(' 1-D well'///)
 124   format(' GP nonlinear'//' Number of atoms =',1pe15.6/      &
              ' a/d0 =',1pe15.6/' phase =',1pe15.6,5x, &
              ' coefficient =',1pe15.6/' shift =',f12.5// )
 131   format(/5f12.5)
 133   format(5f12.5)
 134   format(5e15.6)
 210   format(/' time=',f10.5,' it=',i10)
 211   format(6e15.6)
 212   format(//' check of the eigevalues'/)
 231   format(//' radial mesh'/)
 232   format(5e15.8)
 711   format(' norm=',e15.6/)
 755   format(///' summary of mesh parameters'// &
     ' delt =',e15.6/' delr =',e15.6/' tv =',e15.6///)
 811   format(/////' timimg information '// &
     ' number of time steps =',i10/ &
     ' total time(sec) =',f20.6/ &
     ' average time/step =',f20.12//' digonalization =',f20.12/)
 1201  format(//' error1 - nrp must be even',i10//)
 1331  format('start')
 1332  format('c',e15.6,i10)
 1333  format('c',2i5)
 1466  format(///'error in continuum mesh'// &
     ' npair =',i10,5x,' npair(input_45)=',i10///)
 1666  format('c',6e15.6)
!
       stop
       end
       subroutine funint(nrp1,mask,phase,phase_coeff,x_shift,n_shift,  &
                          iwftyp,t_restart,alpha,velocity)
!
       USE gen_data
       USE splint_arrays
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nrp1,iwftyp,n_shift,mask
       REAL*8, INTENT(IN) :: phase,phase_coeff,velocity,alpha,x_shift
       REAL*8, INTENT(OUT) :: t_restart
       REAL*8 :: pbt,pbb,rexp,r2exp,ebnd,phase_x,repsi,impsi,phasew,&
                 xscl,pi,xxx,yp1,ypn
       INTEGER :: ix,nsrp,ir,itxx,i
!
! set the initial values of the wavefunction
! at time t0 for all x
!
       if(iwftyp == 1) then
          open(44,file='wfn.tmp',status='old')
          read(44,*) nsrp
          read(44,*) ebnd
          do ix = 1,nsrp
             read(44,*) rs(ix),wftmp(ix)
          enddo
!
! spine-fit the wavefunction from bndbs
!
             yp1 = 0.e0
             ypn = 0.e0
             call spline(rs,wftmp,nsrp,yp1,ypn,dy2)
!
! generate the wavefunction on the new radial mesh
!
             do ir = 1,nrp1
                if(xn(ir).ge.rs(1) .and. xn(ir).le.rs(nsrp)) then
                   call  splint(rs,wftmp,dy2,nsrp,xn(ir),wfn_acc(ir))
                else
                   wfn_acc(ir) = 0.
                endif
                write(60,*) xn(ir),wfn_acc(ir)
             enddo
             wfr(1:nrp1) = wfn_acc(1:nrp1)
             wfi(1:nrp1) = 0.e0
       elseif(iwftyp == 2) then
             wftmp(1:nrp1) = exp(-alpha*(xn(1:nrp1)**2))
             wfr(1:nrp1) = cos(velocity*xn(1:nrp1))*wftmp(1:nrp1)
             wfi(1:nrp1) = sin(velocity*xn(1:nrp1))*wftmp(1:nrp1)
       elseif(iwftyp == 3) then
          open(44,file='wfn.tmp',status='old')
          read(44,*) t_restart,itxx
          do ix = 1,nrp1
             read(44,*) rs(ix),wfr(ix),wfi(ix)
          enddo
       endif
!
!
!
!
! construct the pairwise initial wavefunction
!
!
!
! check the normalization of the initial function
!      and re-normalize if necessary
!
!
       if(iwftyp <= 2) then
          call wfnorm(nrp1,pbb,rexp,r2exp)
          write(6,334) pbb
          xscl = 1.e0/sqrt(pbb)
          wfr(1:nrp1) = xscl*wfr(1:nrp1)
          wfi(1:nrp1) = xscl*wfi(1:nrp1)
       endif
!
       call wfnorm(nrp1,pbb,rexp,r2exp)
       write(6,113) pbb
       write(6,779) rexp
       write(6,789) r2exp
!
!
! add in a velocity term
!
!
        if(n_shift > 0 .and. iwftyp <= 2) then
           wfr(1:nrp1) = 0.e0
           wfi(1:nrp1) = 0.e0
           wfn_acc(n_shift:nrp1+n_shift) = wfbnd(1:nrp1)
           wfi(1:nrp1-n_shift) = wfbnd(n_shift:nrp1)
           wfr(1:nrp1) = wfn_acc(1:nrp1) + wfi(1:nrp1) 
           wfi(1:nrp1) = 0.e0
        endif
        write(50,110)
        if(mask == 1) then
           do i = 1,nrp1
              phase_x = 0.5e0*phase*(1.e0 + tanh((xn(i))/phase_coeff))
              write(50,112) xn(i),phase_x
              wfn_acc(i) = sin(phase_x)*wfr(i) + cos(phase_x)*wfi(i)
              wfbnd(i) = cos(phase_x)*wfr(i) - sin(phase_x)*wfi(i)
              wfi(i) = wfn_acc(i)
              wfr(i) = wfbnd(i)
           enddo
        elseif(mask == 2) then
           do i = 1,nrp1
              phase_x = 0.5e0*phase*(1.e0 + tanh((xn(i)-x_shift)/phase_coeff))
              phase_x = phase_x - 0.5e0*phase*(-1.e0 +       &
                           tanh((xn(i)+x_shift)/phase_coeff))
              write(50,112) xn(i),phase_x
              wfn_acc(i) = sin(phase_x)*wfr(i) + cos(phase_x)*wfi(i)
              wfbnd(i) = cos(phase_x)*wfr(i) - sin(phase_x)*wfi(i)
              wfi(i) = wfn_acc(i)
              wfr(i) = wfbnd(i)
           enddo
        elseif(mask == 3) then
           do i = 1,nrp1
              phase_x = 0.5e0*phase*exp(-((xn(i)/phase_coeff)**2))
              write(50,112) xn(i),phase_x
              wfn_acc(i) = sin(phase_x)*wfr(i) + cos(phase_x)*wfi(i)
              wfbnd(i) = cos(phase_x)*wfr(i) - sin(phase_x)*wfi(i)
              wfi(i) = wfn_acc(i)
              wfr(i) = wfbnd(i)
           enddo
        endif
!
!
! print selected values of the wavefunction and probability
!
       write(22,110)
       write(22,111) 
       pi = acos(-one)
       write(6,*) pi, nrp1
       do ix = 1,nrp1
             pbt = (wfr(ix)**2) + (wfi(ix)**2)
             repsi = wfr(ix)
             impsi = wfi(ix)
!             phasew = datan(impsi/repsi)
               if(repsi .lt. 0.e0) then
                  phasew = phasew + pi
               elseif(repsi.gt.0.e0 .and. impsi.lt.0.e0) then
                  phasew = phasew + 2.e0*pi
               endif
             write(22,112) xn(ix),pbt,phasew
       enddo
!
!
 110   format('start')
 111   format('c  initial wf')
 112   format(6e15.6)
 113    format(///' norm at initial time='/3e15.6//)
 334   format(//' initial function overlap='/5e15.6/)
 771   format(////' discrete basis parameters'// &
     ' no. of states =',i10/ &
     ' no. of bound states =',i10/ &
     ' no. of continuum states =',i10///)
 779   format('   <r>=',5e15.6/)
 789   format('   <r2>=',5e15.6/)
!
       return
       end
       subroutine potvdi(scale_ho,rr,vvv)
!
       USE gen_data
!
       IMPLICIT NONE
!
       REAL*8, INTENT(IN) :: rr,scale_ho
       REAL*8, INTENT(OUT) :: vvv
!
!
           vvv = 0.5e0*scale_ho*rr*rr
!
!
      return
      end
       subroutine expk(np)
!
       USE gen_data
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: np
!
       INTEGER :: is,i,j
!
!
!
! Multiply by the real eigenvector matrix U_prop
!
      tm1(1:np) = matmul(U_propT(1:np,1:np),wfr(1:np))
      tm2(1:np) = matmul(U_propT(1:np,1:np),wfi(1:np))
      wfr(1:np) = tm1(1:np)
      wfi(1:np) = tm2(1:np)
!
! Exponential eigenvalue operation:
!
          tm1(1:np) = (wfr(1:np)*upr_k(1:np)) &
                    + (wfi(1:np)*upi_k(1:np))
          tm2(1:np) = (wfi(1:np)*upr_k(1:np)) &
                    - (wfr(1:np)*upi_k(1:np))
          wfr(1:np) = tm1(1:np)
          wfi(1:np) = tm2(1:np)
!
! Multiply by the real eigenvector matrix transpose U_prop
!
      tm1(1:np) = matmul(U_prop(1:np,1:np),wfr(1:np))
      tm2(1:np) = matmul(U_prop(1:np,1:np),wfi(1:np))
      wfr(1:np) = tm1(1:np)
      wfi(1:np) = tm2(1:np)
!
      return
      end
      subroutine expu(np)
!
      USE gen_data
!
      IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: np
!
       INTEGER :: is,i
!
!
! split propagator for the diagonal potential
!
!
! operations performed on pairs of radial points
!
!
!* operate using u on delt
!
          tm1(1:np) = (wfr(1:np)*upr(1:np)) &
                    + (wfi(1:np)*upi(1:np))
          tm2(1:np) = (wfi(1:np)*upr(1:np)) &
                    - (wfr(1:np)*upi(1:np))
          wfr(1:np) = tm1(1:np)
          wfi(1:np) = tm2(1:np)
!
       return
       end
       subroutine timfn(t,delt,np,term_nl)
!
       USE gen_data
!
       IMPLICIT NONE 
!
!
       REAL*8, INTENT(IN) :: t,delt,term_nl
       INTEGER, INTENT(IN) :: np
!
       REAL*8 :: ccc,ex1,ff
!
!
!  nonlinear term
!
         prb1(1:np) = term_nl*((wfr(1:np)**2) + (wfi(1:np)**2))
!
       vijt(1:np) = vpt(1:np) + prb1(1:np)
       upr(1:np) = cos(vijt(1:np))
       upi(1:np) = sin(vijt(1:np))
!
      return
      end
      subroutine prbfn(it,np,iwftyp,t,delr,pbb)
!
       USE gen_data
       USE splint_arrays
!
       IMPLICIT NONE
!
!
       INTEGER, INTENT(IN) :: it,np,iwftyp
       REAL*8, INTENT(IN) :: t,delr
!
       REAL*8 :: pbnd,pbt,rexp,r2exp,pbion,pbiont
       REAL*8 :: tpp,ptot,pbr1,pbi1,pbb,prbr,exp_nlr, &
                 exp_nli,exp_nl,repsi,impsi,phasew,pi,dev
       INTEGER :: ib
!
       pi = acos(-one)
!
           write(6,210) t,it
!
!
! determine the norm of the wavefunction
!
          call wfnorm(np,pbb,rexp,r2exp)
!
          write(6,711) pbb
          write(33,211) t,tpp,pbb
          dev = sqrt(r2exp - (rexp*rexp))
          ptot = 0.e0
          write(6,779) rexp
          write(6,789) r2exp,dev
          write(34,722) t,rexp
!
!
          write(22,1331)
          write(22,1332) t,it
          write(23,1332) t,it
!
!
! print selected values of the wavefunction and probability
!
       pbr(1:np) = (wfr(1:np)**2) + (wfi(1:np)**2)
       do ib = 1,np
          repsi = wfr(ib)
          impsi = wfi(ib)
!          phasew = datan(impsi/repsi)
          if(repsi .lt. 0.e0) then
               phasew = phasew + pi
          elseif(repsi.gt.0.e0 .and. impsi.lt.0.e0) then
               phasew = phasew + 2.e0*pi
          endif
          write(22,211) xn(ib),pbr(ib),phasew
       enddo
       do ib = 1,np
          write(23,*) xn(ib),wfr(ib),wfi(ib)
       enddo
       prb2r(1:np) = pbr(1:np)*wfr(1:np)*wtn(1:np)
       prb2i(1:np) = pbr(1:np)*wfi(1:np)*wtn(1:np)
       exp_nlr = dot_product(wfr(1:np),prb2r(1:np))
       exp_nli = dot_product(wfi(1:np),prb2i(1:np))
       exp_nl = exp_nlr + exp_nli
!       write(6,101) t,exp_nl

!
 101  format(' non-linear strength =',e15.6)
!
!
 210   format(/' time=',f15.5,5x,' it=',i10)
 211   format(6e15.6)
 711   format('   norm=',5e15.6)
 714   format('   total =',e20.12)
 779   format('   <r>=',5e15.6)
 789   format('  <r2>=',1pe15.6,5x,' deviation=',1pe15.6/)
 722   format(6e15.6)
 723   format('    bnd proj =',6e15.6)
 1331  format('start')
 1332  format('c',e15.6,i10)
 1333  format(3e15.6,i10)
!
       return
       end
       subroutine wfnorm(np,pbb,rexp,r2exp)
!
       USE gen_data
       USE splint_arrays
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: np
       REAL*8,  INTENT(OUT) :: pbb,rexp,r2exp
!
!
! determine the norm of the wavefunction
! rexp = expectation value of <x>
! r2exp = sqrt(<x2> - <x>*<x>)
!
             pbb = zero
             rexp = zero
             r2exp = zero
             pb1(1:np) = (wfr(1:np)**2) + (wfi(1:np)**2)
             pbb = dot_product(pb1(1:np),wtn(1:np))
             pb1(1:np) = pb1(1:np)*xn(1:np)
             rexp = dot_product(pb1(1:np),wtn(1:np))
             pb1(1:np) = pb1(1:np)*xn(1:np)
             r2exp = dot_product(pb1(1:np),wtn(1:np))
!
       return
       end
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
!
      USE splint_arrays
!
      INTEGER n
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      INTEGER i,k
      REAL*8 p,qn,sig,un
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
      1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) &
      -sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) write(6,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) &
      *(h** 2)/6.
      return
      END
      SUBROUTINE SPL1D2(N,X,F,W,IJ,Y,TAB)
      implicit real*8 (a-h,o-z)
      DIMENSION X(N),F(N),W(N),TAB(N)
      integer ilo
      save ilo
      data ilo / 1 /
      CALL FndInt ( N-1, X, Y, ilo, I )
      I = MAX0(I,1)
   30 MI=(I-1)*IJ+1
      K1=MI+IJ
      FLK=X(I+1)-X(I)
      A=(W(MI)*(X(I+1)-Y)**3+W(K1)*(Y-X(I))**3)/(6.d0*FLK)
      B=(F(K1)/FLK-W(K1)*FLK/6.d0)*(Y-X(I))
      C=(F(MI)/FLK-FLK*W(MI)/6.d0)*(X(I+1)-Y)
      TAB(1)=A+B+C
      A=(W(K1)*(Y-X(I))**2-W(MI)*(X(I+1)-Y)**2)/(2.d0*FLK)
      B=(F(K1)-F(MI))/FLK
      C=FLK*(W(MI)-W(K1))/6.d0
      TAB(2)=A+B+C
      TAB(3)=(W(MI)*(X(I+1)-Y)+W(K1)*(Y-X(I)))/FLK
      RETURN
      END

      subroutine FndInt ( n, array, value, start, find )
      real*8    array(n)
      real*8    value
      integer high, low, find, start
      high = n
      low = 1
      IF ( value .lt. array(1) ) THEN
         find = 0
         return
      ELSE IF ( value .ge. array(n) ) THEN
         find = n
         return
      END IF
      look = min ( n, max(1,start) )
      IF ( value .lt. array(look) ) THEN
         high = look
         IF ( value .lt. array(high-1) ) THEN
            high = high - 1
         ELSE
            find = high - 1
            return
         END IF
      ELSE
         low = look
         IF ( value .lt. array(low+1) ) THEN
            find = low
            return
         ELSE
            low = low + 1
            IF ( value .lt. array(low+1) ) THEN
               find = low
               return
            ELSE
               low = low + 1
            END IF
         END IF
      END IF
      look = (high + low)/2
      look = min(high-1, look)
      look = max(low +1, look)
    1 continue
      IF ( high - low .eq. 1 ) THEN
         find = low
         return
      ELSE
         IF ( value .lt. array(look) ) THEN
            high = look
         ELSE
            low = look
         END IF
         look = (high + low)/2
         look = min(high-1, look)
         look = max(low +1, look)
      END IF
      go to 1
      end


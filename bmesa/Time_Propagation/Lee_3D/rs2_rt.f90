!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 MODULE gen_data
!
!
   IMPLICIT NONE
!
!
   REAL*8, SAVE :: zero=0.e0, one=1.e0, two=2.e0
   REAL*8, SAVE :: convfs=2.419e-02, convh=27.2114e0
   REAL*8, SAVE :: bohr=0.529177e0, conv1=2.19473e+5
   REAL*8, SAVE :: pi=3.141592653589793
!
!
   REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: xn,yn,zn,wtnx,wtny,wtnz
   REAL*8, SAVE :: delrx,delry,delrz
   REAL*8, SAVE :: amp_opt,width_opt,t_int,t_fin,f_time
   INTEGER, SAVE :: n_shift

END MODULE gen_data
!
MODULE temp_mats
!
! Sets allocatable arrys
!

       REAL*8, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: wfr,wfi
       REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: tm1,tm2,tm3,tm4
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: phxy
       REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: vpt_nl,vpt,upr,  &
                                                  upi,wfnlgr,wfnlgi, &
                                                  vijt,vptr,tm5,tm6
!
END MODULE temp_mats
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       program pdep_ms
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
!
       REAL*8 :: twopi,xstart,xstop,ystart,ystop,zstart,zstop,  &
                 tstart,tstop,delt,atoms_N,a_d0,gauss,gamma,alpha, &
                 beta,chem_g,term_nl,capdl2x,capdl2y,capdl2z,dddx, &
                 dddy,dddz,xxx,tvvdx,tvvdy,tvvdz,tvvd,tt2,tt4,ust, &
                 uuu,cc2x,fc4x,fs4x,cc3x,fc2x,fs2x,cc2y,fc4y,fs4y, &
                 cc3y,fc2y,fs2y,cc2z,fc4z,fs4z,cc3z,fc2z,fs2z,t,   &
                 ttt,tav,v_0,phase,phase_frac,phase_coeff
       REAL*4 :: ttx1,ttx2,tarry(2),etime
       INTEGER :: iflag,nrad,nrpx,nrpy,nrpz,nrp1x,nrp1y,nrp1z,ntop, &
                  ntp,ntime,iwftyp,N_tot,i_count,ix,iy,iz,it,icord, &
                  i_con,ir,mask		  
!
       open (5,file='rs2nl.dat',status='old')
       open (6,file='rs2nl.out',status='unknown')
       open(23,file='wfn2.fle',status='unknown')
       open(31,file='pot.fle',status='unknown')
       open(32,file='prj.fle',status='unknown')
       open(33,file='prb.fle',status='unknown')
       open(34,file='rexp.fle',status='unknown')
       open(35,file='pbe.fle',status='unknown')
       open(42,file='wfnx.fle',status='unknown')
       open(43,file='wfny.fle',status='unknown')
       open(45,file='wfnz.fle',status='unknown')
       open(46,file='contr.fle',status='unknown')
       open(47,file='phase.fle',status='unknown')
       open(48,file='fleld.fle',status='unknown')
       open(49,file='wfnxy_z0.fle',status='unknown', &
               form='unformatted')
!
       twopi = 2.e0*pi
       iflag = 0
!
!
!
!***************************************************************
! DESCRIPTION
!***************************************************************
!
! This program solves a time dependent
! Schrodinger Equation in THREE  spatial dimensions
!
!*** Real Time
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
! nrad = no. of points between successive printings
!         of the solution
! nrpx = no. of points in x mesh( must be EVEN)
! xstart = starting value of x
! xstop = final value of x
!
!
        read(5,*) nrad
        read(5,*) nrpx,xstart,xstop
        read(5,*) nrpy,ystart,ystop
        read(5,*) nrpz,zstart,zstop
!
      if(mod(nrpx,2) .ne. 0 .or. mod(nrpy,2) .ne. 0 .or.   &
         mod(nrpz,2) .ne. 0) then
         write(6,1201) nrpx,nrpy,nrpz
         stop
      endif
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
!!  Trap units -
!
! atoms_N = number of atoms
! a_d0 = scattering length/d_0
! alpha = frequecy ratio(w_y/w_x) for trap potential
! beta  = frequecy ratio(w_z/w_x) for trap potential
! gauss = exponential coefficient of the initial Gaussian function
! chem_g = guess of chemical potential
! 
      read(5,*) atoms_N,a_d0,gauss,gamma,alpha,beta,chem_g,v_0
!
      term_nl = 4.e0*pi*atoms_N*a_d0
!
! phase_frac = imposed phase of the initial wavefunction(solitons)
!               fraction of pi
! phase_coef = "width" of FD-like function in x
! mask = 0  - no phase mask
!        1 - phase mask
!
       read(5,*) phase_frac,phase_coeff,mask,n_shift
       phase = phase_frac*pi
!
!*** Parameters for a radiation interaction term
!
!    amp_opt   = amplitude of the Gauusian laser pulse
!    width_opt = width of the pulse
!    linear ramp off from t_int to t_fin
!
     read(5,*) amp_opt,width_opt,t_int,t_fin
!
!**
!** INPUT: additional parameters
!**
!
! iwftyp = type of initial wavefunction
!          1 - Gaussian packet 
!          2 - Initial Bound wavefuntion(unit 44) real only
!          3 - Initial RT wavefuntion(unit 44) complex
!
       read(5,*)  iwftyp
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
       delrx = (xstop - xstart)/float(nrpx)
       delry = (ystop - ystart)/float(nrpy)
       delrz = (zstop - zstart)/float(nrpz)
!
       if(ntop .eq. 0) then
          delt = (tstop - tstart)/float(ntp)
       else if(ntop .eq. 1) then
          tstop = float(ntp)*delt
       endif
!
       write(6,101) nrpx,xstart,xstop,delrx,          &
        nrpy,ystart,ystop,delry,                      &
        nrpz,zstart,zstop,delrz                      
       write(6,102) ntp,tstart,tstop,delt
       write(6,103) atoms_N,a_d0,alpha,beta,term_nl,chem_g,phase,&
                    phase_coeff
       write(6,104) amp_opt,width_opt,t_int,t_fin
       write(6,122)
!
!
! set the radial mesh
!
!
       capdl2x = 1.e0/(2.e0*(delrx**2))
       capdl2y = 1.e0/(2.e0*(delry**2))
       capdl2z = 1.e0/(2.e0*(delrz**2))
!
       nrp1x = nrpx + 1
       nrp1y = nrpy + 1
       nrp1z = nrpz + 1
       if(nrpx .eq. 2) then
            nrpx = 1
            nrp1x = 1
       endif
       if(nrpy .eq. 2) then
            nrpy = 1
            nrp1y = 1
       endif
       if(nrpz .eq. 2) then
            nrpz = 1
            nrp1z = 1
       endif
       
       dddx = delt*capdl2x
       dddy = delt*capdl2y
       dddz = delt*capdl2z
       N_tot = nrp1x*nrp1y*nrp1z
!
       write(6,755) delt,delrx,dddx,delry,dddy,delrz,dddz
!
!
! set-up a pairwise radial array - xn
!       and the weights for spatial integration - wtn
!
       ALLOCATE(xn(nrp1x+1),wtnx(nrp1x+1),yn(nrp1y+1),   &
                wtny(nrp1y+1),zn(nrp1z+1),wtnz(nrp1z+1))
!
       xxx = xstart - delrx
       do ir = 1,nrp1x
          xxx = xxx + delrx
          xn(ir) = xxx
          wtnx(ir) = delrx
       enddo
       xxx = ystart - delry
       do ir = 1,nrp1y
          xxx = xxx + delry
          yn(ir) = xxx
          wtny(ir) = delry
       enddo
       xxx = zstart - delrz
       do ir = 1,nrp1z
          xxx = xxx + delrz
          zn(ir) = xxx
          wtnz(ir) = delrz
       enddo
!
       tvvdx = 2.*capdl2x
       tvvdy = 2.*capdl2y
       tvvdz = 2.*capdl2z
       if(nrp1x .eq. 1) then
          wtnx(1) = 1.e0
          tvvdx = 0.e0
       endif
       if(nrp1y .eq. 1) then
          wtny(1) = 1.e0
          tvvdy = 0.e0
       endif
       if(nrp1z .eq. 1) then
          wtnz(1) = 1.e0
          tvvdz = 0.e0
       endif
!
!
       tvvd = tvvdx + tvvdy + tvvdz
       tt2 = delt/2.
       tt4 = delt/4.
!
! DIAGONAL, TIME-INDEPENDENT POTENTIAL and PROPAGATOR
!
!
!! ALLOCATION: vpt,vptr,upr,wfr,wfnlg,tm1,tm2
!
        ALLOCATE(vpt(1:N_tot+1),vptr(1:N_tot+1),upr(1:N_tot+1),      &
                 wfnlgr(1:N_tot+1),vpt_nl(1:N_tot+1),vijt(1:N_tot+1), &
                 tm5(1:N_tot+1),upi(1:N_tot+1),tm6(1:N_tot+1),        &
                 wfnlgi(1:N_tot+1))
        ALLOCATE(wfr(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 wfi(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 tm1(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 tm2(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 tm3(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 tm4(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1))
        ALLOCATE(phxy(1:nrp1x+1,1:nrp1y+1))    
!!
!
        i_count = 0
!        write(6,*) N_tot,nrp1x
        do ix = 1,nrp1x
           do iy = 1,nrp1y
              do iz = 1,nrp1z
                 i_count = i_count + 1
                 call pot_trap(xn(ix),yn(iy),zn(iz),gamma,  &
                               alpha,beta,ust)
                 uuu = (ust + tvvd)*delt
                 vpt(i_count) = uuu
                 vptr(i_count) = ust
              enddo
           enddo
        enddo
!
       upr(1:N_tot) = cos(vpt(1:N_tot))
       upi(1:N_tot) = sin(vpt(1:N_tot))
!
!
! set the initial conditions at the first time point t0
!
!
       call funint(nrp1x,nrp1y,nrp1z,iwftyp,gauss,v_0,phase,  &
                   phase_coeff,mask)
!       write(6,*) 'wfr',wfr(1:5,1:5,1:5)
       i_con = 1
       call conv_3d1d(nrp1x,nrp1y,nrp1z,N_tot,i_con)
!       write(6,*) 'wfnlg',wfnlg(1:100)
!
!
!
!
! KINETIC ENERGY OPERATOR: PARAMETERS B
!
! set the cosine/sine coefficients
!
       cc2x = tt4*capdl2x
       fc4x = cos(cc2x)
       fs4x = sin(cc2x)
       cc3x = tt2*capdl2x
       fc2x = cos(cc3x)
       fs2x = sin(cc3x)
       cc2y = tt4*capdl2y
       fc4y = cos(cc2y)
       fs4y = sin(cc2y)
       cc3y = tt2*capdl2y
       fc2y = cos(cc3y)
       fs2y = sin(cc3y)
       cc2z = tt4*capdl2z
       fc4z = cos(cc2z)
       fs4z = sin(cc2z)
       cc3z = tt2*capdl2z
       fc2z = cos(cc3z)
       fs2z = sin(cc3z)
!
!
!
!
       write(32,1331)
       write(33,1331)
       write(34,1331)
       write(48,1331)
       write(6,105)
!
!
! propagation in time
!
!
       t = tstart
       ttx1 = etime(tarry)
!
!********************************************
!
! TEMPORAL PROPAGATION
!
!********************************************
!
       do it = 1,ntp
          t = t + delt
!
!
! Determine the Time-dependent Potential
!
          call timfn(t,delt,N_tot,nrp1x,nrp1y,nrp1z,   &
                     term_nl,chem_g)
!          write(6,*) 'timfn',vpt(1:10),vptr(1:10)
!          write(6,*) upr(1:10)
!          write(6,*) 'wfr',wfr(nrp1x/2:nrp1x/2+10,1,1)
!          write(6,*) 'wflgn',wfnlg(nrp1x/2:nrp1x/2+10)
!
!
!  Propogate the Kinetic Energy in x
!
          if(nrp1x .gt. 1) then
             icord = 1
             call expk(nrp1x,nrp1y,nrp1z,icord,fc2x,fs2x,fc4x,fs4x) 
          endif
!
!  Propogate the Kinetic Energy in y
!
          if(nrp1y .gt. 1) then
             icord = 2
             call expk(nrp1x,nrp1y,nrp1z,icord,fc2y,fs2y,fc4y,fs4y) 
          endif
!
!  Propogate the Kinetic Energy in z
!
          if(nrp1z .gt. 1) then
             icord = 3
             call expk(nrp1x,nrp1y,nrp1z,icord,fc2z,fs2z,fc4z,fs4z) 
          endif
!
          i_con = 1
          call conv_3d1d(nrp1x,nrp1y,nrp1z,N_tot,i_con)
!
!
! Propagate the diagonal potential operator
!
          call expu(N_tot)
          i_con = 2
          call conv_3d1d(nrp1x,nrp1y,nrp1z,N_tot,i_con)
!
!  Propogate the Kinetic Energy in x
!
          if(nrp1x .gt. 1) then
             icord = 1
             call expk(nrp1x,nrp1y,nrp1z,icord,fc2x,fs2x,fc4x,fs4x) 
          endif
!
!  Propogate the Kinetic Energy in y
!
          if(nrp1y .gt. 1) then
             icord = 2
             call expk(nrp1x,nrp1y,nrp1z,icord,fc2y,fs2y,fc4y,fs4y) 
          endif
!
!  Propogate the Kinetic Energy in z
!
          if(nrp1z .gt. 1) then
             icord = 3
             call expk(nrp1x,nrp1y,nrp1z,icord,fc2z,fs2z,fc4z,fs4z) 
          endif
!
!
          if(mod(it,ntime) .eq. 0) then
              write(46,1341) it,t
              write(49) it,t
              call prbfn(it,ntp,nrp1x,nrp1y,nrp1z,iwftyp,t,term_nl)  
          endif
!
!
       enddo
!
!
! Write to unit 23 the wavefunction at the final time step
!
       do ix = 1,nrp1x
          do iy = 1,nrp1y
             do iz = 1,nrp1z
                 write(23,*) xn(ix),yn(iy),zn(iz),wfr(ix,iy,iz),&
                             wfi(ix,iy,iz)
             enddo
          enddo
       enddo
!
       ttx2 = etime(tarry) - ttx1
       ttt = ttx2
       tav = ttt/ntp
       write(6,811) ntp,N_tot,ttt,tav
!
!
 100   format(//' pde solution by product formula'///)
 101   format(// &
     ' radial mesh'//' no. of pomits=',i10/ &
     ' rstart=',e15.6,5x,' rstop=',e15.6/' delr=',e15.6)
 102   format(//' temporal mesh'//' no. of ponts=',i10/ &
     ' tstart=',f10.5,' tstop=',f15.5/' delt=',e15.6/)
 103   format(///' POTENTIAL PARAMETERS'//' Number of atoms ='1pe15.6/  &
                   ' a/d0 =',1pe15.6,5x,'alpha =',1pe15.6,             &
                   5x,' beta =',1pe15.6/                               &
                   ' non-linear term =',1pe15.6//     &
                   ' Chemical Potential(Initial Guess) =',1pe15.6/     &
                   ' Phase(initial) =',1pe15.6/' Coefficient =',1pe15.6/)
 104   format(//' Laser field interaction'/    &
                  ' Amplitude =',1pe15.6 &
                  ' Width =',1pe15.6/    &
                  ' Initial ramp time  =',1pe15.6, &
                  ' Final ramp time  =',1pe15.6//)
 105   format(///' propagation in r and t'//)
 106   format(i10)
 122   format(//' Potential Type'/)
 711   format(' norm=',e15.6/)
 755   format(///' summary of mesh parameters'// &
     ' delt =',e15.6/' delrx =',e15.6,5x,' tvx =',e15.6/ &
     ' delry =',e15.6,5x,' tvy =',e15.6/     &
     ' delrz =',e15.6,5x,' tvz =',e15.6///)     
 811   format(/////' timimg information '// &
     ' number of time steps =',i10/   &
     ' number of radial pts. ='i10/ &
     ' total time(sec) =',f20.6/ &
     ' average time/step =',f20.12/)
 1201  format(//' error1 - nrp must be even',3i10//)
 1331  format('start')
 1332  format(2i10)
 1333  format(e15.6)
 1341  format('start'/'c ',i10,f12.5)
 1666  format('c',6e15.6)
!
       stop
       end
       subroutine funint(nrp1x,nrp1y,nrp1z,iwftyp,alpha,v_0,phase, &
                         phase_coeff,mask)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nrp1x,nrp1y,nrp1z,iwftyp,mask
       REAL*8, INTENT(IN) :: alpha,v_0,phase, phase_coeff
       REAL*8 ::  pbt,pbb,rexp,r2exp,ebnd,xexp,yexp,zexp,scl, &
                  csx,ssx,xnorm,xx,yy,zz,x2exp,y2exp,z2exp,   &
                  cpx,spx,mean_field,chem1,phase_x,pi2,x_shift,&
                  y_shift,sm2,dxx,dyy
       INTEGER :: ix,iy,iz,ixs,iys
!
       if(iwftyp == 1) then
             xnorm = ((2.e0*alpha)/pi)**0.75
             do ix = 1,nrp1x
                 csx = cos(v_0*xn(ix))
                 ssx = sin(v_0*xn(ix))
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      wfr(ix,iy,iz) =  exp(-alpha*(xn(ix)**2   &
                                       + yn(iy)**2 + zn(iz)**2))
                      wfi(ix,iy,iz) = xnorm*ssx*wfr(ix,iy,iz)
                      wfr(ix,iy,iz) = xnorm*csx*wfr(ix,iy,iz)
                   enddo
                enddo
              enddo
       else if(iwftyp == 2) then
             open(44,file='wfn.tmp',status='old')
             read(44,*) chem1
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                       read(44,*) xx,yy,zz,wfr(ix,iy,iz)
                        wfi(ix,iy,iz) = 0.e0
                   enddo
                enddo
              enddo
       else if(iwftyp == 3) then
             open(44,file='wfn.tmp',status='old')
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                       read(44,*) xx,yy,zz,wfr(ix,iy,iz), &
                                  wfi(ix,iy,iz)
                   enddo
                enddo
              enddo
       endif
       call flush(6)
!
!
!! Scale the Initial Wavefunction by a phase exp(i*phase)
!
    if(mask == 1) then
       do ix = 1,nrp1x
          phase_x = phase/(exp(phase_coeff*xn(ix)) +1.e0)
          write(47,780) xn(ix),phase_x
          cpx = cos(phase_x)
          spx = sin(phase_x)
          do iy = 1,nrp1y
             do iz = 1,nrp1z
                wfi(ix,iy,iz) = spx*wfr(ix,iy,iz)
                wfr(ix,iy,iz) = cpx*wfr(ix,iy,iz)
             enddo
          enddo
       enddo
    elseif (mask == 2) then
       pi2 = 2.e0*pi
       do iz = 1,nrp1z
          do ix = 1,nrp1x
             do iy = 1,nrp1y
                phase_x = (phase/pi2)*datan2(yn(iy),xn(ix))
                cpx = cos(phase_x)
                spx = sin(phase_x)
                wfi(ix,iy,iz) = spx*wfr(ix,iy,iz)
                wfr(ix,iy,iz) = cpx*wfr(ix,iy,iz)
             enddo
          enddo
       enddo
    elseif (mask == 3) then
       pi2 = 2.e0*pi
       do iz = 1,nrp1z
          do ix = 1,nrp1x
             do iy = 1,nrp1y
                phase_x = 0.e0
                do ixs = -n_shift,n_shift
                   x_shift = ixs*delrx
                   do iys = -n_shift,n_shift
                      y_shift = iys*delry
                      phase_x = phase_x + (phase/pi2)*  &
                                datan2(yn(iy)+y_shift,xn(ix)+x_shift)
                   enddo
                enddo
                phase_x = phase_x/((2*n_shift+1)**2)
                cpx = cos(phase_x)
                spx = sin(phase_x)
                wfi(ix,iy,iz) = spx*wfr(ix,iy,iz)
                wfr(ix,iy,iz) = cpx*wfr(ix,iy,iz)
             enddo
          enddo
       enddo
    elseif (mask == 4) then
       sm2 = 2.e+0/(pi*phase_coeff)
          do ix = 1,nrp1x
             do iy = 1,nrp1y
                phase_x = 0.e0
                do ixs = 1,nrp1x
                   do iys = 1,nrp1y
                      dxx = xn(ixs) - xn(ix)
                      dyy = yn(iys) - yn(iy)
                      phase_x = phase_x + datan2(dyy,dxx)*  &
                                exp(-(xn(ixs)**2 + yn(iys)**2)/phase_coeff)*&
                                wtnx(ixs)*wtny(iys)
                   enddo
                enddo
                phxy(ix,iy) = phase_x*(sm2)
             enddo
          enddo
       do iz = 1,nrp1z
          do ix = 1,nrp1x
             do iy = 1,nrp1y
                phase_x = phxy(ix,iy)
                cpx = cos(phase_x)
                spx = sin(phase_x)
                wfi(ix,iy,iz) = spx*wfr(ix,iy,iz)
                wfr(ix,iy,iz) = cpx*wfr(ix,iy,iz)
             enddo
          enddo
       enddo
    endif
!
! check the normalization of the initial function
!      and re-normalize if necessary
!
!
       call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp, &
                   x2exp,y2exp,z2exp,mean_field)
       write(6,113) pbb
       write(6,779) xexp,yexp,zexp,x2exp,y2exp,z2exp,mean_field
       scl = 1.e0/sqrt(pbb)
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      wfr(ix,iy,iz) = wfr(ix,iy,iz)*scl 
                      wfi(ix,iy,iz) = wfi(ix,iy,iz)*scl
                   enddo
                enddo
              enddo
       call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp, &
                   x2exp,y2exp,z2exp,mean_field)
       write(6,113) pbb
       write(6,779) xexp,yexp,zexp,x2exp,y2exp,z2exp,mean_field
!
!
!
 113    format(///' norm at initial time='/3e15.6//)
 779   format(' <x>=',e15.6,5x,'<y>=',e15.6,5x,'<z>=',e15.6/ &
              ' <x2>=',e15.6,5x,'<y2>=',e15.6,5x,'<z2>=',e15.6/ &
              ' Mean_field =',e15.6//) 
 780   format(2e15.6)
!
       return
       end
       subroutine pot_trap(xx,yy,zz,gamma,alpha,beta,vvv)
!
       USE gen_data
!
       IMPLICIT NONE
!
       REAL*8, INTENT(IN) :: xx,yy,zz,alpha,beta,gamma
       REAL*8, INTENT(OUT) :: vvv
!
!
           vvv = 0.5e0*((gamma*xx)**2 + (alpha*yy)**2 +(beta*zz)**2)
!
!
      return
      end
      subroutine expk(nrp1x,nrp1y,nrp1z,icord,c2,s2,c4,s4)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nrp1x,nrp1y,nrp1z,icord
       REAL*8, INTENT(IN) :: c2,s2,c4,s4
!
       INTEGER :: ixo0,ixm,ixc,ixo2,ixe0,ixe2,ixem,ixemm,ixom, &
                  iyo0,iym,iyc,iyo2,iye0,iye2,iyem,iyemm,iyom, &
                  izo0,izm,izc,izo2,ize0,ize2,izem,izemm,izom
!
!
! split propagator for the kinetic energy
!
!
       if(icord .eq. 1) then
          ixo0 = 1
          ixm = nrp1x
          ixc = 2
          ixo2 = 2
          ixe0 = 2
          ixe2 = 3
          iyo0 = 1
          iym = nrp1y
          iyc = 1
          iyo2 = 1
          iye0 = 1
          iye2 = 1
          izo0 = 1
          izm = nrp1z
          izc = 1
          izo2 = 1
          ize0 = 1
          ize2 = 1
          ixom = nrp1x
          iyom=1
          izom = 1
          ixem = 1
          ixemm = 1
          iyem = 1
          iyemm = nrp1y
          izem = 1
          izemm = nrp1z
       else if(icord .eq. 2) then
          ixo0 = 1
          ixm = nrp1x
          ixc = 1
          ixo2 = 1
          ixe0 = 1
          ixe2 = 1
          iyo0 = 1
          iym = nrp1y
          iyc = 2
          iyo2 = 2
          iye0 = 2
          iye2 = 3
          izo0 = 1
          izm = nrp1z
          izc = 1
          izo2 = 1
          ize0 = 1
          ize2 = 1
          ixom = 1
          iyom = nrp1y
          izom = 1
          ixem = 1
          ixemm = nrp1x
          iyem = 1
          iyemm = 1
          izem = 1
          izemm = nrp1z
       else if(icord .eq. 3) then
          ixo0 = 1
          ixm = nrp1x
          ixc = 1
          ixo2 = 1
          ixe0 = 1
          ixe2 = 1
          iyo0 = 1
          iym = nrp1y
          iyc = 1
          iyo2 = 1
          iye0 = 1
          iye2 = 1
          izo0 = 1
          izm = nrp1z
          izc = 2
          izo2 = 2
          ize0 = 2
          ize2 = 3
          ixom = 1
          iyom = 1
          izom = nrp1z
          ixem = 1
          ixemm = nrp1x
          iyem = 1
          iyemm = nrp1y
          izem = 1
          izemm = 1
       endif
!
! operations performed on pairs of radial points
!
!
!* operate using k(odd) on delt/4
!
         tm1(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =   &
                     c4*wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   - s4*wfi(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm2(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =   &
                     c4*wfi(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   + s4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm3(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                     -s4*wfi(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   + c4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm4(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                     +s4*wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   + c4*wfi(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =    &
                        tm1(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) 
         wfi(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =    &
                        tm2(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) 
         wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                        tm3(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfi(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                        tm4(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixom:ixm:ixc,iyom:iym:iyc,izom:izm:izc) = zero
         wfi(ixom:ixm:ixc,iyom:iym:iyc,izom:izm:izc) = zero
!
!
!
!
!* operate with k(even) for delt/2
!
!
         wfr(ixem:ixemm:ixc,iyem:iyemm:iyc,izem:izemm:izc) = zero
         wfi(ixem:ixemm:ixc,iyem:iyemm:iyc,izem:izemm:izc) = zero
         tm1(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) =    &
                     c2*wfr(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc)& 
                   - s2*wfi(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 
         tm2(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) =    &
                     c2*wfi(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc)& 
                   + s2*wfr(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 
         tm3(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) =    &
                     -s2*wfi(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc)& 
                   + c2*wfr(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 
         tm4(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) =    &
                     s2*wfr(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc)& 
                   + c2*wfi(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 
         wfr(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) =   &
                        tm1(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) 
         wfi(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) =   &
                        tm2(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) 
         wfr(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) =  &
                        tm3(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 
         wfi(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) =  &
                        tm4(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 

!
!
!
!* operate with k(odd) on delt/4
!
         tm1(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =   &
                     c4*wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   - s4*wfi(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm2(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =   &
                     c4*wfi(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   + s4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm3(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                     -s4*wfi(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   + c4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm4(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                     +s4*wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   + c4*wfi(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =    &
                        tm1(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) 
         wfi(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =    &
                        tm2(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) 
         wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                        tm3(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfi(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                        tm4(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixom:ixm:ixc,iyom:iym:iyc,izom:izm:izc) = zero
         wfi(ixom:ixm:ixc,iyom:iym:iyc,izom:izm:izc) = zero
!
!
      return
      end
      subroutine expu(N_tot)
!
      USE gen_data
      USE temp_mats
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: N_tot
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
          tm5(1:N_tot) = (wfnlgr(1:N_tot)*upr(1:N_tot)) &
                        + (wfnlgi(1:N_tot)*upi(1:N_tot)) 
          tm6(1:N_tot) = (wfnlgi(1:N_tot)*upr(1:N_tot)) &
                        - (wfnlgr(1:N_tot)*upi(1:N_tot)) 
          wfnlgr(1:N_tot) = tm5(1:N_tot)
          wfnlgi(1:N_tot) = tm6(1:N_tot)
!
       return
       end
       subroutine timfn(t,delt,N_tot,nrp1x,nrp1y,nrp1z, &
                  term_nl,chem_g)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE 
!
!
       REAL*8, INTENT(IN) :: t,term_nl,delt
       INTEGER, INTENT(IN) :: nrp1x,nrp1y,nrp1z,N_tot
       REAL*8, INTENT(INOUT) :: chem_g
!
       REAL*8 :: pbb,xexp,yexp,zexp,x2exp,y2exp,z2exp, & 
                 ex1,ff,v_opt,v_time,mean_field,rho2,amt
       INTEGER :: i,ix,iy,iz,i_count
!
!
       call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp, &
                   x2exp,y2exp,z2exp,mean_field)
       chem_g = chem_g/pbb
!
! add temporal dependence
!
       if(t <= t_int) then
          f_time = 1.e0
       elseif(t>t_int .and. t<=t_fin) then
          f_time = (t_fin - t)/(t_fin - t_int)
       elseif(t > t_fin) then
          f_time = 0.e0
       endif
!
       i_count = 0
       amt = f_time*amp_opt
       do ix = 1,nrp1x
          do iy = 1,nrp1y
              do iz = 1,nrp1z
                 rho2 = xn(ix)**2 + yn(iy)**2
                 v_opt = f_time*amp_opt*exp(-width_opt*rho2)
                 i_count = i_count + 1
                 vpt_nl(i_count) = term_nl*(wfr(ix,iy,iz)**2 &
                                   +  wfi(ix,iy,iz)**2) +    &
                                    v_opt
              enddo
          enddo
       enddo
       vijt(1:N_tot) = vpt(1:N_tot) + (vpt_nl(1:N_tot)*delt)
       upr(1:N_tot) = cos(vijt(1:N_tot))
       upi(1:N_tot) = sin(vijt(1:N_tot))
!
 100  format(2e15.6)
!
      return
      end
      subroutine prbfn(it,ntp,nrp1x,nrp1y,nrp1z,iwftyp,t,term_nl)  
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
!
       INTEGER, INTENT(IN) :: it,nrp1x,nrp1y,nrp1z,iwftyp,ntp
       REAL*8, INTENT(IN) :: t
!
       REAL*8 :: pbb,xexp,yexp,zexp,x2exp,y2exp,z2exp, &
                 pbbx,mean_field,term_nl,fr_x,fr_y,fr_z,&
                 fi_x,fi_y,fi_z,Lr_x,Li_x,L_x,Lr_y,Li_y,&
                 L_y,Lr_z,Li_z,L_z,L_tot,delx2,dely2,delz2,&
                 Lt_x,Lt_y,Lt_z,amt
       INTEGER :: ib,ix,iy,iz,i_count,nx2,ny2,nz2,i,j,k
!
!
           write(6,210) t,it,ntp
!
!
! determine the norm of the wavefunction
!
          call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp, &
                      x2exp,y2exp,z2exp,mean_field)
!
!
! Mean Field energy
!
          mean_field = mean_field*term_nl
          write(6,711) pbb
          write(33,211) t,pbb
          write(6,779) xexp,yexp,zexp,x2exp,y2exp,z2exp,mean_field
          write(34,780) t,xexp,yexp,zexp,x2exp,y2exp,z2exp,mean_field
!
!
!
! Print to unit 42 the x-direction of the probabiltiy
!
            ny2 = nrp1y/2 + 1
            nz2 = nrp1z/2 + 1
            write(42,100)
            do ix = 1,nrp1x
               pbbx = wfr(ix,ny2,nz2)**2 + wfi(ix,ny2,nz2)**2
               write(42,101) xn(ix),pbbx
            enddo
!
! Print to unit 43 the y-direction of the probabiltiy
!
            nx2 = nrp1x/2 + 1
            nz2 = nrp1z/2 + 1
            write(43,100)
            do iy = 1,nrp1y
               pbbx = wfr(nx2,iy,nz2)**2 + wfi(nx2,iy,nz2)**2
               write(43,101) yn(iy),pbbx
            enddo
!
!
! Print to unit 44 the z-direction of the probabiltiy
!
            nx2 = nrp1x/2 + 1
            ny2 = nrp1y/2 + 1
            write(45,100)
            do iz = 1,nrp1z
               pbbx = wfr(nx2,ny2,iz)**2 + wfi(nx2,ny2,iz)**2
               write(45,101) zn(iz),pbbx
            enddo
!
!
! Integrate in z direction and write resulting function
!
            do ix = 1,nrp1x
                do iy = 1,nrp1y
                    xexp = 0.e0
                    do iz = 1,nrp1z
                       pbbx = wfr(ix,iy,iz)**2 + wfi(ix,iy,iz)**2
                       xexp = xexp + (pbbx*wtnz(iz))
                    enddo
                    write(46,101) xn(ix),yn(iy),xexp
                 enddo
             enddo
!
!! Determine the expectation value of the angular momnetum
!
      delx2 = 1.e0/(2.e0*delrx)
      dely2 = 1.e0/(2.e0*delry)
      delz2 = 1.e0/(2.e0*delrz)
      L_x = 0.e0
      L_y = 0.e0
      L_z = 0.e0
!
!   <Lz>
!
      do k = 1,nrp1z
         do i = 1,nrp1x
            do j = 1,nrp1y
               if(i == 1) then
                  fr_x = (wfr(i+1,j,k) - wfr(i,j,k))/delrx
                  fi_x = (wfi(i+1,j,k) - wfi(i,j,k))/delrx
               elseif(i >1 .and. i < nrp1x) then
                  fr_x = (wfr(i+1,j,k) - wfr(i-1,j,k))*delx2
                  fi_x = (wfi(i+1,j,k) - wfi(i-1,j,k))*delx2
               elseif(i == nrp1x) then
                  fr_x = (wfr(i,j,k) - wfr(i-1,j,k))/delrx
                  fi_x = (wfi(i,j,k) - wfi(i-1,j,k))/delrx
               endif
               if(j == 1) then
                  fr_y = (wfr(i,j+1,k) - wfr(i,j,k))/delry
                  fi_y = (wfi(i,j+1,k) - wfi(i,j,k))/delry
               elseif(j >1 .and. j < nrp1y) then
                  fr_y = (wfr(i,j+1,k) - wfr(i,j-1,k))*dely2
                  fi_y = (wfi(i,j+1,k) - wfi(i,j-1,k))*dely2
               elseif(j == nrp1y) then
                  fr_y = (wfr(i,j,k) - wfr(i,j-1,k))/delry
                  fi_y = (wfi(i,j,k) - wfi(i,j-1,k))/delry
                endif
                Lr_z = xn(i)*fr_y - (yn(j)*fr_x)
                Li_z = xn(i)*fi_y - (yn(j)*fi_x)
                Lt_z =  ((wfr(i,j,k)*Lr_z) + (wfi(i,j,k)*Li_z)) +&
                            ((wfr(i,j,k)*Li_z) - (wfi(i,j,k)*Lr_z))
                L_z = L_z + (Lt_z*wtnx(i)*wtny(j)*wtnz(k))
            enddo
         enddo
      enddo
!
!   <Ly>
!
      do j = 1,nrp1y
         do i = 1,nrp1x
            do k = 1,nrp1z
            if(i == 1) then
               fr_x = (wfr(i+1,j,k) - wfr(i,j,k))/delrx
               fi_x = (wfi(i+1,j,k) - wfi(i,j,k))/delrx
            elseif(i >1 .and. i < nrp1x) then
               fr_x = (wfr(i+1,j,k) - wfr(i-1,j,k))*delx2
               fi_x = (wfi(i+1,j,k) - wfi(i-1,j,k))*delx2
            elseif(i == nrp1x) then
               fr_x = (wfr(i,j,k) - wfr(i-1,j,k))/delrx
               fi_x = (wfi(i,j,k) - wfi(i-1,j,k))/delrx
            endif
               if(k == 1) then
                  fr_z = (wfr(i,j,k+1) - wfr(i,j,k))/delrz
                  fi_z = (wfi(i,j,k+1) - wfi(i,j,k))/delrz
               elseif(k >1 .and. k < nrp1z) then
                  fr_z = (wfr(i,j,k+1) - wfr(i,j,k-1))*delz2
                  fi_z = (wfi(i,j,k+1) - wfi(i,j,k-1))*delz2
               elseif(k == nrp1z) then
                  fr_z = (wfr(i,j,k) - wfr(i,j,k-1))/delrz
                  fi_z = (wfi(i,j,k) - wfi(i,j,k-1))/delrz
                endif
                Lr_y = zn(k)*fr_x - (xn(i)*fr_z)
                Li_y = zn(k)*fi_x - (xn(i)*fi_z)
                Lt_y = ((wfr(i,j,k)*Lr_y) + (wfi(i,j,k)*Li_y)) +&
                            ((wfr(i,j,k)*Li_y) - (wfi(i,j,k)*Lr_y))
                L_y = L_y + (Lt_y*wtnx(i)*wtny(j)*wtnz(k))
            enddo
         enddo
      enddo
!
!   <Lx>
!
      do i = 1,nrp1x
         do j = 1,nrp1y
            do k = 1,nrp1z
            if(j == 1) then
               fr_y = (wfr(i,j+1,k) - wfr(i,j,k))/delry
               fi_y = (wfi(i,j+1,k) - wfi(i,j,k))/delry
            elseif(j >1 .and. j < nrp1y) then
               fr_y = (wfr(i,j+1,k) - wfr(i,j-1,k))*dely2
               fi_y = (wfi(i,j+1,k) - wfi(i,j-1,k))*dely2
            elseif(j == nrp1y) then
               fr_y = (wfr(i,j,k) - wfr(i,j-1,k))/delry
               fi_y = (wfi(i,j,k) - wfi(i,j-1,k))/delry
            endif
               if(k == 1) then
                  fr_z = (wfr(i,j,k+1) - wfr(i,j,k))/delrz
                  fi_z = (wfi(i,j,k+1) - wfi(i,j,k))/delrz
               elseif(k >1 .and. k < nrp1z) then
                  fr_z = (wfr(i,j,k+1) - wfr(i,j,k-1))*delz2
                  fi_z = (wfi(i,j,k+1) - wfi(i,j,k-1))*delz2
               elseif(k == nrp1z) then
                  fr_z = (wfr(i,j,k) - wfr(i,j,k-1))/delrz
                  fi_z = (wfi(i,j,k) - wfi(i,j,k-1))/delrz
                endif
                Lr_x = yn(j)*fr_z - (zn(k)*fr_y)
                Li_x = yn(j)*fi_z - (zn(k)*fi_y)
                Lt_x =  ((wfr(i,j,k)*Lr_x) + (wfi(i,j,k)*Li_x)) +&
                            ((wfr(i,j,k)*Li_x) - (wfi(i,j,k)*Lr_x))
                L_x = L_x + (Lt_x*wtnx(i)*wtny(j)*wtnz(k))
            enddo
         enddo
      enddo
!
! total L
!
      L_tot = sqrt(L_x**2 + L_y**2 + L_z**2)
      write(6,781) L_tot,L_x,L_y,L_z,f_time
      call flush(6)
!
! Write the wavefunction
!  and  the velocity field in x and y at z=0
!
            nz2 = nrp1z/2 + 1
            do ix = 1,nrp1x
                do iy = 1,nrp1y
                    write(49) xn(ix),yn(iy),wfr(ix,iy,nz2),wfi(ix,iy,nz2)
                enddo
            enddo
            call flush(49)
!
!
 100   format('start')
 101   format(6e15.6)
 210   format(/' time=',f15.5,5x,' it=',i10,5x, &
              ' ntime =',i10)
 211   format(6e15.6)
 711   format('   norm=',5e15.6)
 714   format('   total =',e20.12)
 779   format(' <x>=',e15.6,5x,'<y>=',e15.6,5x,'<z>=',e15.6/ &
              ' <x2>=',e15.6,5x,'<y2>=',e15.6,5x,'<z2>=',e15.6/ &
              ' Mean_field =',e15.6/) 
 780   format(9e9.3)
 781   format(' Total Angular Momentum =',e15.6/3e15.6,5x,2f15.6//)
!
       return
       end
       subroutine wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp, &
                         x2exp,y2exp,z2exp,mean_field)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nrp1x,nrp1y,nrp1z
       REAL*8,  INTENT(OUT) :: pbb,xexp,yexp,zexp,x2exp,y2exp,z2exp
!
       REAL*8 :: pbb1,pbb2,pbbx,mean_field
       INTEGER :: ix,iy,iz,ny2,nz2
!
!
! determine the norm of the wavefunction
! rexp = expectation value of <x>
!
             pbb = zero
             mean_field = 0.e0
             xexp = zero
             yexp = zero
             zexp = zero
             x2exp = zero
             y2exp = zero
             z2exp = zero
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      pbb2 = (wfr(ix,iy,iz)**2 + wfi(ix,iy,iz)**2)
                      pbb1 = pbb2*wtnx(ix)*wtny(iy)*wtnz(iz)
                      pbb = pbb + pbb1
                      xexp = xexp + xn(ix)*pbb1
                      yexp = yexp + yn(iy)*pbb1
                      zexp = zexp + zn(iz)*pbb1
                      x2exp = x2exp + xn(ix)*xn(ix)*pbb1
                      y2exp = y2exp + yn(iy)*yn(iy)*pbb1
                      z2exp = z2exp + zn(iz)*zn(iz)*pbb1
                      mean_field = mean_field + pbb1*pbb2
                   enddo
                enddo
             enddo
!
!
       return
       end
       subroutine conv_3d1d(nx,ny,nz,ntot,icon)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nx,ny,nz,icon,ntot
       INTEGER :: i_count,ix,iy,iz
!
! Convert from 1-d form to and from 3-d form of wavefunction
!
! icon = 1  - FROM wfr TO wfnlg
!        2  - FROM wfnlg TO wfr
!
       if(icon .eq. 1) then
           i_count = 0
           do ix = 1,nx
              do iy = 1,ny
                 do iz = 1,nz
                    i_count = i_count + 1
                    wfnlgr(i_count) = wfr(ix,iy,iz)
                    wfnlgi(i_count) = wfi(ix,iy,iz)
                 enddo
              enddo
           enddo
       else if(icon .eq. 2) then
           i_count = 0
           do ix = 1,nx
              do iy = 1,ny
                 do iz = 1,nz
                    i_count = i_count + 1
                    wfr(ix,iy,iz) = wfnlgr(i_count)
                    wfi(ix,iy,iz) = wfnlgi(i_count)
                 enddo
              enddo
           enddo
       endif
!
       return
       end

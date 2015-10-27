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
!
!
   REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: xn,yn,zn,wtnx,wtny,wtnz
   REAL*8, SAVE :: delrx,delry,delrz
   REAL*8, SAVE :: omega_opt,amp_opt,t0_opt,width_opt
!
END MODULE gen_data
!
MODULE temp_mats
!
! Sets allocatable arrys
!

       REAL*8, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: wfr,wfr0
       REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: tm1,tm2
       REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: vpt_nl,vpt,upr,wfnlg, &
                                                  vijt,vptr,tm3
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
       REAL*8, PARAMETER :: tol=1.e-08
!
!
       REAL*8 :: pi,twopi,xstart,xstop,ystart,ystop,zstart,zstop,  &
                 tstart,tstop,delt,atoms_N,a_d0,gauss,gamma,alpha, &
                 beta,chem_g,term_nl,capdl2x,capdl2y,capdl2z,dddx, &
                 dddy,dddz,xxx,tvvdx,tvvdy,tvvdz,tvvd,tt2,tt4,ust, &
                 uuu,cc2x,fc4x,fs4x,cc3x,fc2x,fs2x,cc2y,fc4y,fs4y, &
                 cc3y,fc2y,fs2y,cc2z,fc4z,fs4z,cc3z,fc2z,fs2z,t,   &
                 ttt,tav,etot,phase,iorth
       REAL*4 :: ttx1,ttx2,tarry(2),etime
       INTEGER :: iflag,nrad,nrpx,nrpy,nrpz,nrp1x,nrp1y,nrp1z,ntop, &
                  ntp,ntime,iwftyp,N_tot,i_count,ix,iy,iz,it,icord, &
                  i_con,ir
!
       open (5,file='rs2nl.dat',status='old')
       open (6,file='rs2nl.out',status='unknown')
       open(22,file='wfn1.fle',status='unknown')
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
       open(60,file='wav.fle',status='unknown')
!
       pi = acos(-one)
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
!*** Imaginary Time
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
! phase = z+ phase*pi
! iorth = 1 - force othogonality
! 
      read(5,*) atoms_N,a_d0,gauss,gamma,alpha,beta,chem_g,phase,iorth
!
      term_nl = 4.e0*pi*atoms_N*a_d0
      phase = phase*pi
!
!
!*** Parameters for a radiation interaction term
!
!    omega_opt = frequency in units of H_bar*omega_x
!    amp_opt   = amplitude of EM wave
!    t0_opt    = shift in the ramp function
!    width_opt = width in ramp function
!
     read(5,*) omega_opt,amp_opt,t0_opt,width_opt
!
!**
!** INPUT: additional parameters
!**
!
! iwftyp = type of initial wavefunction
!          1 - 3d Gaussian with width "gauss"
!          2 - input from unit 44 wfn.tmp
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
       write(6,103) atoms_N,a_d0,alpha,beta,term_nl,chem_g
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
                 wfnlg(1:N_tot+1),vpt_nl(1:N_tot+1),vijt(1:N_tot+1), &
                 tm3(1:N_tot+1))
        ALLOCATE(wfr(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 tm1(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 tm2(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1),    &
                 wfr0(1:nrp1x+1,1:nrp1y+1,1:nrp1z+1))
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
!        write(6,*) 'vpt',vpt(1:100)
!
       upr(1:N_tot) = exp(-vpt(1:N_tot))
!       write(6,*) 'upr',upr(1:100)
!
!
!
! set the initial conditions at the first time point t0
!
!
       call funint(nrp1x,nrp1y,nrp1z,iwftyp,gauss,chem_g,phase,iorth)
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
       fc4x = cosh(cc2x)
       fs4x = sinh(cc2x)
       cc3x = tt2*capdl2x
       fc2x = cosh(cc3x)
       fs2x = sinh(cc3x)
       cc2y = tt4*capdl2y
       fc4y = cosh(cc2y)
       fs4y = sinh(cc2y)
       cc3y = tt2*capdl2y
       fc2y = cosh(cc3y)
       fs2y = sinh(cc3y)
       cc2z = tt4*capdl2z
       fc4z = cosh(cc2z)
       fs4z = sinh(cc2z)
       cc3z = tt2*capdl2z
       fc2z = cosh(cc3z)
       fs2z = sinh(cc3z)
!
!
!
!
       write(32,1331)
       write(33,1331)
       write(34,1331)
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
          if(mod(it,ntime) .eq. 0) then
              call prbfn(it,ntp,nrp1x,nrp1y,nrp1z,iwftyp,t,  &
                         term_nl,atoms_N,chem_g,etot)
          endif
!
! Force orthogonality
!
          if (iorth .eq. 1) then
             call Orthog(nrp1x,nrp1y,nrp1z)
          endif
!
!
       enddo
!
! Write to unit 23 the wavefunction at the final time step
!
       write(23,*) chem_g
       do ix = 1,nrp1x
          do iy = 1,nrp1y
             do iz = 1,nrp1z
                 write(23,*) xn(ix),yn(iy),zn(iz),wfr(ix,iy,iz)
             enddo
          enddo
       enddo
!
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
                ' Chemical Potential(Initial Guess) =',1pe15.6//)  
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
 1666  format('c',6e15.6)
!
       stop
       end
       subroutine funint(nrp1x,nrp1y,nrp1z,iwftyp,alpha,chem_g, &
                         phase,iorth)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nrp1x,nrp1y,nrp1z,iwftyp,iorth
       REAL*8, INTENT(IN) :: alpha,phase
       REAL*8, INTENT(OUT) :: chem_g
       REAL*8 ::  pbt,pbb,rexp,r2exp,ebnd,xexp,yexp,zexp,scl,  &
                  xx,yy,zz,chem_g0,cpx
       INTEGER :: ix,iy,iz
!
!
       if(iwftyp .eq. 1) then
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      wfr(ix,iy,iz) =  exp(-alpha*(xn(ix)**2   &
                                       + yn(iy)**2 + zn(iz)**2))
                   enddo
                enddo
              enddo
       else if(iwftyp .eq. 2) then
             open(44,file='wfn.tmp',status='old')
             read(44,*) chem_g
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                       read(44,*) xx,yy,zz,wfr(ix,iy,iz)
                   enddo
                enddo
              enddo
       endif
!
! Read in the orbital to which orthogonality forced
!
       if(iorth == 1) then
             open(50,file='wfn.orth',status='old')
             read(50,*) chem_g0
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                       read(50,*) xx,yy,zz,wfr0(ix,iy,iz)
                   enddo
                enddo
              enddo
        endif
!
!
!! Scale the Initial Wavefunction by a phase exp(i*phase)
!
       do ix = 1,nrp1x
          cpx = 1.e0
          if(xn(ix) .gt. 0.e0) then
             cpx = cos(phase)
          endif
          do iy = 1,nrp1y
             do iz = 1,nrp1z
                wfr(ix,iy,iz) = cpx*wfr(ix,iy,iz)
             enddo
          enddo
       enddo
!
!
! check the normalization of the initial function
!      and re-normalize if necessary
!
     if(iwftyp .eq. 1) then
       call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp)
       write(6,113) pbb
       write(6,779) xexp,yexp,zexp
       scl = 1.e0/sqrt(pbb)
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      wfr(ix,iy,iz) = wfr(ix,iy,iz)*scl 
                   enddo
                enddo
              enddo
       call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp)
       write(6,113) pbb
       write(6,779) xexp,yexp,zexp
     endif
!
!
!
 113    format(///' norm at initial time='/3e15.6//)
 779   format('   <x>=',e15.6/      &
              '   <y>=',e15.6/      &
              '   <z>=',e15.6////)      
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
                   + s4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm2(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                     s4*wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc)& 
                   + c4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =    &
                        tm1(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) 
         wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =    &
                        tm2(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixom:ixm:ixc,iyom:iym:iyc,izom:izm:izc) = zero
!
!
!
!
!* operate with k(even) for delt/2
!
!
         wfr(ixem:ixemm:ixc,iyem:iyemm:iyc,izem:izemm:izc) = zero
         tm1(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) =    &
                     c2*wfr(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc)& 
                   + s2*wfr(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 
         tm2(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) =    &
                     s2*wfr(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc)& 
                   + c2*wfr(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 
         wfr(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) =   &
                        tm1(ixe0:ixm:ixc,iye0:iym:iyc,ize0:izm:izc) 
         wfr(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) =  &
                        tm2(ixe2:ixm:ixc,iye2:iym:iyc,ize2:izm:izc) 

!
!
!
!* operate with k(odd) on delt/4
!
         tm1(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =   &
                     c4*wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) &
                   + s4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         tm2(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =  &
                     s4*wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) &
                   + c4*wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) =   &
                        tm1(ixo0:ixm:ixc,iyo0:iym:iyc,izo0:izm:izc) 
         wfr(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) =   &
                        tm2(ixo2:ixm:ixc,iyo2:iym:iyc,izo2:izm:izc) 
         wfr(ixom:ixm:ixc,iyom:iym:iyc,izom:izm:izc) = zero
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
          tm3(1:N_tot) = (wfnlg(1:N_tot)*upr(1:N_tot)) 
          wfnlg(1:N_tot) = tm3(1:N_tot)
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
       REAL*8 :: pbb,xexp,yexp,zexp,ex1,ff,v_opt,rho2
       INTEGER :: i,ix,iy,iz,i_count
!
!
       call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp)
       chem_g = chem_g/pbb
!
! add temporal dependence
!
       i_count = 0
       do ix = 1,nrp1x
          do iy = 1,nrp1y
              do iz = 1,nrp1z
                 rho2 = xn(ix)**2 + yn(iy)**2
                 v_opt = amp_opt*exp(-width_opt*rho2)
                 i_count = i_count + 1
                 vpt_nl(i_count) = (term_nl*((wfr(ix,iy,iz))**2) +&
                                    v_opt)
              enddo
          enddo
       enddo
       vijt(1:N_tot) = vpt(1:N_tot) + ((vpt_nl(1:N_tot) - chem_g)*delt)
       upr(1:N_tot) = exp(-vijt(1:N_tot))
!
      return
      end
      subroutine prbfn(it,ntp,nrp1x,nrp1y,nrp1z,iwftyp,t,  &
                       term_nl,atoms_N,chem_g,etot)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
!
       INTEGER, INTENT(IN) :: it,nrp1x,nrp1y,nrp1z,iwftyp,ntp
       REAL*8, INTENT(IN) :: t,term_nl,atoms_N,chem_g
       REAL*8, INTENT(OUT) :: etot
!
       REAL*8 :: pbb,xexp,yexp,zexp,wf2,dr,atoms_eff,wtt
       REAL*8 :: eng_ham,eng_hamv,eng_hamt,eng_hamtx,     &
                 eng_hamty,eng_hamtz,pbbx
       INTEGER :: ib,ix,iy,iz,i_count,nx2,ny2,nz2
!
!
           write(6,210) t,it,ntp
!
!
! determine the norm of the wavefunction
!
          call wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp)
!
          write(6,711) pbb
          write(33,211) t,pbb
          write(6,779) xexp,yexp,zexp
!
!
! Print to unit 42 the x-direction of the probabiltiy
!
            ny2 = nrp1y/2 + 1
            nz2 = nrp1z/2 + 1
            write(42,100)
            do ix = 1,nrp1x
               pbbx = wfr(ix,ny2,nz2)**2 
               write(42,101) xn(ix),pbbx
            enddo
!
! Print to unit 43 the y-direction of the probabiltiy
!
            nx2 = nrp1x/2 + 1
            nz2 = nrp1z/2 + 1
            write(43,100)
            do iy = 1,nrp1y
               pbbx = wfr(nx2,iy,nz2)**2 
               write(43,101) yn(iy),pbbx
            enddo
!
!
! Print to unit 45 the z-direction of the probabiltiy
!
            nx2 = nrp1x/2 + 1
            ny2 = nrp1y/2 + 1
            write(45,100)
            do iz = 1,nrp1z
               pbbx = wfr(nx2,ny2,iz)**2 
               write(45,101) zn(iz),pbbx
            enddo
!
!
! Print to unit 46 the contour z=0  of the probabiltiy
!
            nz2 = nrp1z/2 + 1
            write(46,100)
            do ix = 1,nrp1x
               do iy = 1,nrp1y
                  pbbx = wfr(ix,iy,nz2)**2 
                  write(46,101) xn(ix),yn(iy),pbbx
              enddo
           enddo
!
!
!
! Expectation value of d^2/dx^2
!
          dr = delrx
          eng_hamtx = 0.e0
          if(nrp1x .gt. 1) then
             do ix = 1,nrp1x 
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      wtt = wtnx(ix)*wtny(iy)*wtnz(iz)
                      wf2 = 0.e0
                      if(ix .eq. 1) then
                         wf2=(wfr(3,iy,iz) - 2.0e0*wfr(2,iy,iz)  &
                                        + wfr(1,iy,iz))/(dr**2) 
                      elseif(ix .gt. 1 .and. ix .lt. nrp1x-1) then
                          wf2=(wfr(ix+1,iy,iz) - 2.e0*wfr(ix,iy,iz) &
                                         + wfr(ix-1,iy,iz))/(dr**2)
                      elseif(ix .eq. nrp1x) then
                          wf2 = (wfr(nrp1x,iy,iz)  &
                                 - 2.e0*wfr(nrp1x-1,iy,iz) &
                                 + wfr(nrp1x-2,iy,iz))/(dr**2)
                      endif
                      eng_hamtx = eng_hamtx + (wfr(ix,iy,iz)*wf2*wtt)
                   enddo
                enddo
             enddo
          endif
!
! Expectation value of d^2/dy^2
!
          dr = delry
          eng_hamty = 0.e0
          if(nrp1y .gt. 1) then
             do ix = 1,nrp1x 
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      wtt = wtnx(ix)*wtny(iy)*wtnz(iz)
                      wf2 = 0.e0
                      if(iy .eq. 1) then
                         wf2=(wfr(ix,3,iz) - 2.0e0*wfr(ix,2,iz)  &
                                        + wfr(ix,1,iz))/(dr**2) 
                      elseif(iy .gt. 1 .and. iy .lt. nrp1y-1) then
                          wf2=(wfr(ix,iy+1,iz) - 2.e0*wfr(ix,iy,iz) &
                                         + wfr(ix,iy-1,iz))/(dr**2)
                      elseif(iy .eq. nrp1y) then
                          wf2 = (wfr(ix,nrp1y,iz)  &
                                 - 2.e0*wfr(ix,nrp1y-1,iz) &
                                 + wfr(ix,nrp1y-2,iz))/(dr**2)
                      endif
                      eng_hamty = eng_hamty + (wfr(ix,iy,iz)*wf2*wtt)
                   enddo
                enddo
             enddo
          endif
!
! Expectation value of d^2/dz^2
!
          dr = delrz
          eng_hamtz = 0.e0
          if(nrp1z .gt. 1) then
             do ix = 1,nrp1x 
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      wtt = wtnx(ix)*wtny(iy)*wtnz(iz)
                      wf2 = 0.e0
                      if(iz .eq. 1) then
                         wf2=(wfr(ix,iy,3) - 2.0e0*wfr(ix,iy,2)  &
                                        + wfr(ix,iy,1))/(dr**2) 
                      elseif(iz .gt. 1 .and. iz .lt. nrp1z-1) then
                          wf2=(wfr(ix,iy,iz+1) - 2.e0*wfr(ix,iy,iz) &
                                         + wfr(ix,iy,iz-1))/(dr**2)
                      elseif(iz .eq. nrp1z) then
                          wf2 = (wfr(ix,iy,nrp1z)  &
                                 - 2.e0*wfr(ix,iy,nrp1z-1) &
                                 + wfr(ix,iy,nrp1z-2))/(dr**2)
                      endif
                      eng_hamtz = eng_hamtz + (wfr(ix,iy,iz)*wf2*wtt)
                   enddo
                enddo
             enddo
          endif
!
          eng_hamt = -0.5e0*(eng_hamtx + eng_hamty + eng_hamtz)
!
! Expectation value of the potential
!
             eng_hamv = 0.e0 
             i_count = 0
             do ix = 1,nrp1x 
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      i_count = i_count + 1
                      wtt = wtnx(ix)*wtny(iy)*wtnz(iz)
                      wf2 = wfr(ix,iy,iz)**2
                      eng_hamv = eng_hamv +                &
                         wf2*(vptr(i_count) + (term_nl*wf2))*wtt
                   enddo
                enddo
             enddo
!
             eng_ham = (eng_hamt + eng_hamv)/pbb

          atoms_eff = pbb*atoms_N
          etot = eng_ham
          write(6,723) atoms_eff,eng_ham,eng_hamt,eng_hamv,chem_g
          write(32,722) it,t,pbb,eng_ham,atoms_eff,chem_g
!          write(6,*) eng_hamtx,eng_hamty,eng_hamtz
!
!
!
!
 100   format('start')
 101   format(3e15.6)
 210   format(/' time=',f15.5,5x,' it=',i10,5x, &
              ' ntime =',i10)
 211   format(6e15.6)
 711   format('   norm=',5e15.6)
 714   format('   total =',e20.12)
 779   format('   <r>=',5e15.6//)
 722   format(i10,5e15.6)
 723   format(//' Effective Number of Trap atoms =',1pe15.5//   &
                   '   energy (ham-exp) =',1pe15.6//1p2e15.6//  &
                   '   chemical potential =',1pe15.6/)
 724   format(5e15.6)
!
       return
       end
       subroutine wfnorm(nrp1x,nrp1y,nrp1z,pbb,xexp,yexp,zexp)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nrp1x,nrp1y,nrp1z
       REAL*8,  INTENT(OUT) :: pbb,xexp,yexp,zexp
!
       REAL*8 :: pbb1
       INTEGER :: ix,iy,iz
!
!
! determine the norm of the wavefunction
! rexp = expectation value of <x>
!
             pbb = zero
             xexp = zero
             yexp = zero
             zexp = zero
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      pbb1 = (wfr(ix,iy,iz)**2)*wtnx(ix)*   &
                              wtny(iy)*wtnz(iz)
                      pbb = pbb + pbb1
                      xexp = xexp + xn(ix)*pbb1
                      yexp = yexp + yn(iy)*pbb1
                      zexp = zexp + zn(iz)*pbb1
                   enddo
                enddo
             enddo
!
            if(abs(pbb) .gt. 0.e0) then
                 xexp = xexp/pbb
                 yexp = yexp/pbb
                 zexp = zexp/pbb
            endif
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
                    wfnlg(i_count) = wfr(ix,iy,iz)
                 enddo
              enddo
           enddo
       else if(icon .eq. 2) then
           i_count = 0
           do ix = 1,nx
              do iy = 1,ny
                 do iz = 1,nz
                    i_count = i_count + 1
                    wfr(ix,iy,iz) = wfnlg(i_count)
                 enddo
              enddo
           enddo
       endif
!
       return
       end
       subroutine Orthog(nrp1x,nrp1y,nrp1z)
!
       USE gen_data
       USE temp_mats
!
       IMPLICIT NONE
!
       INTEGER, INTENT(IN) :: nrp1x,nrp1y,nrp1z
!
       REAL*8 :: ovr,ovr1,xnn,ovr_old
       INTEGER :: ix,iy,iz
!
!
! determine the overlap of the wavefunctions
!
             ovr = zero
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      ovr1 = wfr(ix,iy,iz)*wtnx(ix)*wtny(iy)* &
                             wtnz(iz)*wfr0(ix,iy,iz)
                      ovr = ovr + ovr1
                   enddo
                enddo
             enddo
             ovr_old = ovr
!
! Construct the orthonormal wavefunction
!
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                       wfr(ix,iy,iz) = wfr(ix,iy,iz) -  &
                             ovr*wfr0(ix,iy,iz)
                   enddo
                enddo
             enddo
!
! Check orthogonality
!
             ovr = zero
             do ix = 1,nrp1x
                do iy = 1,nrp1y
                   do iz = 1,nrp1z
                      ovr1 = wfr(ix,iy,iz)*wtnx(ix)*wtny(iy)* &
                             wtnz(iz)*wfr0(ix,iy,iz)
                      ovr = ovr + ovr1
                   enddo
                enddo
             enddo
!
      write(6,100) ovr_old,ovr
 100  format(/' overlap =',1pe15.6/' orthogonality =',1pe15.6)
!
       return
       end



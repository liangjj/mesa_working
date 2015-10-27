*deck @(#)rhomax.f	5.2 2/5/95
      function rhomax(type,atnum)
c***begin prologue     rhomax.f
c***date written       950120  
c***revision date      2/5/95      
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             @(#)rhomax.f	5.2   2/5/95
c***purpose            
c***description
c     
c                      return an estimate of the maximum in the
c                      atomic density assuming either slater's
c                      rules, values from clementi, 
c                      or bragg-slater radii.
c
c                      the data was taken from G92/DFT.
c
c***references
c
c***routines called
c
c***end prologue       rhomax.f
      implicit none
      real*8 rhomax
c     --- input variables -----
      integer atnum
      character*(*) type
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      logical called
      integer itype
      real*8 rhomx1(3,18), rhomx2(3,18), rhomx3(2,18), rhomx4(2,17)
      real*8 rhomx5(2,15), rhomx6(8)
      real*8 toang
c
      data called/.false./
      save called,rhomx1, rhomx2, rhomx3, rhomx4, rhomx5, rhomx6
c
c     --- distances at which the densities exhibit maxima.
c                  Slater-Bragg,         Clementi,   Slater's Rules,
      data rhomx1 / 0.529177d0,          0.529177d0,   0.529177d0,      H
     $              0.31d0,              0.31d0,       0.31d0,          He*
     $              1.45d0,              1.67d0,       1.63d0,          Li
     $              1.05d0,              1.12d0,       1.08d0,          Be
     $              0.85d0,              0.87d0,       0.81d0,          B
     $              0.70d0,              0.67d0,       0.65d0,          C
     $              0.65d0,              0.56d0,       0.54d0,          N
     $              0.60d0,              0.48d0,       0.47d0,          O
     $              0.50d0,              0.42d0,       0.40d0,          F
     $              0.38d0,              0.38d0,       0.36d0,          Ne*
     $              1.80d0,              1.90d0,       2.16d0,          Na
     $              1.50d0,              1.45d0,       1.67d0,          Mg
     $              1.25d0,              1.18d0,       1.36d0,          Al
     $              1.10d0,              1.11d0,       1.15d0,          Si
     $              1.00d0,              0.98d0,       0.99d0,          P
     $              1.00d0,              0.88d0,       0.87d0,          S
     $              1.00d0,              0.79d0,       0.78d0,          Cl
     $              0.71d0,              0.71d0,       0.70d0/          Ar*
c
      data rhomx2 / 2.20d0,              2.43d0,       3.32d0,          K
     $              1.80d0,              1.94d0,       2.56d0,          Ca
     $              1.60d0,              1.84d0,       2.43d0,          Sc
     $              1.40d0,              1.76d0,       2.32d0,          Ti
     $              1.35d0,              1.71d0,       2.22d0,          V
     $              1.40d0,              1.66d0,       2.12d0,          Cr
     $              1.40d0,              1.61d0,       2.02d0,          Mn
     $              1.40d0,              1.56d0,       1.95d0,          Fe
     $              1.35d0,              1.52d0,       1.87d0,          Co
     $              1.35d0,              1.49d0,       1.80d0,          Ni
     $              1.35d0,              1.45d0,       1.73d0,          Cu
     $              1.35d0,              1.42d0,       1.67d0,          Zn
     $              1.30d0,              1.36d0,       1.46d0,          Ga
     $              1.25d0,              1.25d0,       1.29d0,          Ge
     $              1.15d0,              1.14d0,       1.16d0,          As
     $              1.15d0,              1.03d0,       1.05d0,          Se
     $              1.15d0,              0.94d0,       0.96d0,          Br
     $              0.88d0,              0.88d0,       0.88d0/          Kr*
c
      data rhomx3 / 2.35d0,              2.65d0,                        Rb
     $              2.00d0,              2.19d0,                        Sr
     $              1.80d0,              2.12d0,                        Y
     $              1.55d0,              2.06d0,                        Zr
     $              1.45d0,              1.98d0,                        Nb
     $              1.45d0,              1.90d0,                        Mo
     $              1.35d0,              1.83d0,                        Tc
     $              1.30d0,              1.78d0,                        Ru
     $              1.35d0,              1.73d0,                        Rh
     $              1.40d0,              1.69d0,                        Pd
     $              1.60d0,              1.65d0,                        Ag
     $              1.55d0,              1.61d0,                        Cd
     $              1.55d0,              1.56d0,                        In
     $              1.45d0,              1.45d0,                        Sn
     $              1.45d0,              1.33d0,                        Sb
     $              1.40d0,              1.23d0,                        Te
     $              1.40d0,              1.15d0,                        I
     $              1.08d0,              1.08d0/                        Xe*
c
      data rhomx4 / 2.60d0,              2.98d0,                        Cs
     $              2.15d0,              2.53d0,                        Ba
     $              1.95d0,              6.22d0,                        La
     $              1.85d0,              5.05d0,                        Ce
     $              1.85d0,              2.47d0,                        Pr
     $              1.85d0,              2.06d0,                        Nd
     $              1.85d0,              2.05d0,                        Pm
     $              1.85d0,              2.38d0,                        Sm
     $              1.85d0,              2.31d0,                        Eu
     $              1.80d0,              2.33d0,                        Gd
     $              1.75d0,              2.25d0,                        Tb
     $              1.75d0,              2.28d0,                        Dy
     $              1.75d0,              2.26d0,                        Ho
     $              1.75d0,              2.26d0,                        Er
     $              1.75d0,              2.22d0,                        Tm
     $              1.75d0,              2.22d0,                        Yb
     $              1.75d0,              2.17d0/                        Lu
c
      data rhomx5 / 1.55d0,              2.08d0,                        Hf
     $              1.45d0,              2.00d0,                        Ta
     $              1.35d0,              1.93d0,                        W
     $              1.35d0,              1.88d0,                        Re
     $              1.30d0,              1.85d0,                        Os
     $              1.35d0,              1.80d0,                        Ir
     $              1.35d0,              1.77d0,                        Pt
     $              1.35d0,              1.74d0,                        Au
     $              1.50d0,              1.71d0,                        Hg
     $              1.90d0,              1.56d0,                        Tl
     $              1.80d0,              1.54d0,                        Pb
     $              1.60d0,              1.43d0,                        Bi
     $              1.90d0,              1.35d0,                        Po
     $              1.27d0,              1.27d0,                        At*
     $              1.20d0,              1.20d0/                        Rn*
c
c                   Slater-Bragg
c                   ?d0,                                                Fr
      data rhomx6 / 2.15d0,                                             Ra
     $              1.95d0,                                             Ac
     $              1.80d0,                                             Th
     $              1.75d0,                                             Pa
     $              1.75d0,                                             U
     $              1.75d0,                                             Np
     $              1.75d0,                                             Pu
     $              1.75d0/                                             Am
c
c     $              ?d0,                                               Cm
c     $              ?d0,                                               Bk
c     $              ?d0,                                               Cf
c     $              ?d0,                                               Es
c     $              ?d0,                                               Fm
c     $              ?d0,                                               Md
c     $              ?d0,                                               No
c     $              ?d0,                                               Lr
c
c
      if(.not.called) then
         call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
         called=.true.
      endif
c
c
      if(type.eq.'bragg') then
         itype=1
      else if(type.eq.'clementi') then
         itype=2
         if(atnum.gt.86) then
            call lnkerr('clementi radii not available for Z>86')
         endif
      else if(type.eq.'slater') then
         itype=3
         if(atnum.gt.36) then
            call lnkerr('slater radii not available for Z>36')
         endif
      endif
c
c     --- assign the radii
      if(atnum.le.18) then
         rhomax = rhomx1(itype,atnum)/toang
      else if(atnum.le.36) then
         rhomax = rhomx2(itype,atnum-18)/toang
      else if(atnum.le.54) then
         rhomax = rhomx3(itype,atnum-36)/toang
      else if(atnum.le.71) then
         rhomax = rhomx4(itype,atnum-54)/toang
      else if(atnum.le.86) then
         rhomax = rhomx5(itype,atnum-71)/toang
      else if(atnum.eq.87) then
         call lnkerr('radii missing for Z=87 in rhomax')
      else if(atnum.le.95) then
         rhomax = rhomx6(atnum-87)/toang
      else
         call lnkerr('radii missing for Z>95 in rhomax')
      endif
c
c
      return
      end

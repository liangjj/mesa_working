$route

$end
          
$nonstd
11//04;
20//01;
$end

$title
time-propagation code
$end

$dvrprop_basis
 number-of-space-variables=2 coordinate-system=spherical
 space-variable-1=r_1 space-variable-2=r_2 use-atomic-units 
 kinetic-energy-type=packed
 xplot get-eigenpairs xprint=all-details keep-diagonals plot_step=20
$end

$h0(r_1)
 automate number-of-major-blocks=4
 number-of-fixed-points=2 left-fixed-point right-fixed-point 
 drop-first-function drop-last-function xdo-not-diagonalize
 region-boundaries=(1.d-15,300.d0) reuse-space-data
 print=(eigenvalues) xsector-print=sector-details 
$end

 $block-1
  default-order=10 number-of-subregions=1
  left-boundary=1.d-15 right-boundary=5.0d0
 $end

 $block-2
  default-order=5 number-of-subregions=3
  left-boundary=5.d0 right-boundary=10.d0
 $end

 $block-3
  default-order=4 number-of-subregions=25
  left-boundary=10.d0 right-boundary=50.d0
 $end

 $block-4
  default-order=3 number-of-subregions=25
  left-boundary=50.d0 right-boundary=300.d0
 $end

$v_reg_1(r_1)
 use-atomic-units potential=coulomb
$end

$v_reg_2(r_1)
 use-atomic-units potential=coulomb
$end

$v_reg_3(r_1)
 use-atomic-units potential=coulomb
$end

$h0(r_2)
 automate number-of-major-blocks=1
 number-of-fixed-points=2 left-fixed-point right-fixed-point 
 drop-first-function drop-last-function xdo-not-diagonalize
 region-boundaries=(1.d-15,200.d0) reuse-space-data
 print=(xall) xsector-print=sector-details 
$end

 $block-1
  default-order=10 number-of-subregions=1
 left-boundary=1.d-15 right-boundary=5.0d0
 $end

 $block-2
  default-order=5 number-of-subregions=3
  left-boundary=5.d0 right-boundary=10.d0
 $end

 $block-3
  default-order=3 number-of-subregions=25
  left-boundary=10.d0 right-boundary=50.d0
 $end

 $block-4
  default-order=3 number-of-subregions=25
  left-boundary=50.d0 right-boundary=200.d0
 $end

$v_reg_1(r_2)
 use-atomic-units potential=none
$end
 
$time
 automate number-of-time-regions=100 first-time=0.0 
 time-interval=.001d0
 xprint=(main=(xpointers,xpotential,xinitial-state,solution,
              xnon-linear-potential,xh-on-initial-state))
$end

$v0(t1)
 potential=none
$end

$initial-state
 driver=gaussian-pulse sigma=(1.d0,1.d0) alpha=(2.d0,2.d0) 
 beta=(0.d0,0.d0) x_0=(0.d0,0.d0) xprint=on
$end

$v_couple
  v-space-time=none
$end


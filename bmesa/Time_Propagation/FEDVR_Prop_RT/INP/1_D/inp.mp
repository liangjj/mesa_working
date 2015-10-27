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
 number-of-space-variables=1 coordinate-system=cartesian
 space-variable-1=x xunits=atomic-units kinetic-energy-type=packed
 plot get-eigenpairs diagonal-modification=add_potential_energy
 plot_step=100 print=all-details
 propagation-order=4 xnon-linear-equations
$end

$h0(x)
 automate number-of-major-blocks=3
 number-of-fixed-points=2 left-fixed-point right-fixed-point 
 drop-left-function drop-right-function xdo-not-diagonalize
 order-of-finite-difference-formula=3 number-of-steps=50
 region-boundaries=(0.d0,10.d0)
 print=all sector-print=sector-details xreuse-space-data
$end

 $block-1
  default-order=10 number-of-subregions=1
  left-boundary=0.d0 right-boundary=1.d0
 $end

$v_reg_1(x)
 use-atomic-units potential=well well-depth=-2.5d0
 well-size=1.d0
 $end

  $block-2
  default-order=150 number-of-subregions=1
  left-boundary=1.d0 right-boundary=2.d0
 $end

$v_reg_2(x)
 use-atomic-units potential=well well-depth=0.0d0
 $end
   $block-3
  default-order=180 number-of-subregions=1
  left-boundary=2.d0 right-boundary=10.d0
 $end

$v_reg_3(x)
 use-atomic-units potential=well well-depth=0.0d0
$end

$time
 number-of-time-regions=10 initial-time=0.d0 time-interval=.001d0
 xprint=(main=(xpointers,xpotential,xinitial-state,xnon-linear-potential,
              xh-on-initial-state))
$end

$v0(t1)
 potential=none
$end

$initial-state
 driver=perturbed-state-vector initial-state=0  
$end

$v_couple
  v-space-time=none
$end

$nlse
 oscillator-frequencies-in-hertz=(27.15d0) number-of-atoms=1.d+06
 mass=3.8192d-26  scattering-length=2.79d-09
$end

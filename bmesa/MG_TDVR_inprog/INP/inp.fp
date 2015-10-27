$route
 number-of-space-dimensions=1 print=(sector=all,tdvr=all) 
 use-atomic-units coordinate-system= cartesian dimension-1=x 
 number-of-time-regions=1 xiterative-linear-system-solve m8001=punch
$end
          
$nonstd
80//01;
20//01;
$end

$title
time-propagation code
$end

$h0(x)
 number-of-regions=1 region-boundaries=(-.5d0,.5d0)
 polynomial-order-per-region=(10) 
$end

$v0(x)
 use-atomic-units potential=none
$end 
 
$h0(t1)
 region-boundaries=(0.d0,1.d0) polynomial-order=10
 print=(sector=xall,tdvr=xall)
$end

$v0(t1)
 potential=none
$end

$gmres
 overlap-tolerance=1.e-10 convergence=1.e-08
 preconditioner=block maximum-size-of-preconditioning-block=50
 xprint=(iterative-solve=all) maximum-number-of-iterations=3
 maximum-number-of-vectors=2 
$end

$trials
 type-of-trial-vectors=unit number-of-trial-vectors=1
$end

$initial-state
 driver=gaussian-pulse sigma=1.0 beta=4.0
$end

$vtnlse(1,1,t1)

$end

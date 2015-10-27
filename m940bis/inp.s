$route
2s+1=1
hf
scf=(pulay,convergence=5)
maxsiz=2000000
guess=(chk) drop
ksym=(symmetry=1,psalc)
print=(guess=all,scf=vector,nodis,noang,kohn=p-space)
geom=(coord,inau,nocrowd)
drt=(sets=2,key-1=target,key-2=p-and-q)
target=(ngroups=5,nrefs=5) p-and-q=(ngroups=5,nrefs=6)
ci=(target,form,diagonalize,target-roots=6,cycles=60) internal-orbitals=3
sym=norotate
int=(drt=(key=(target)),reuse,zint) 
kohn=(nsmall=3,ncsfs=10,nl2=3,freeze,hqq,
        energy=(.4410,.55125,.66152,.735))
$end
 
$nonstd
1//1,2;
3//1,2,12,30;
4//1;
8//6,11,19,22;
9//40,21,29;
20//1;
$end
 
$title
H2 Two State Calculation
$end
 
$geom
h1(basis=small) 0.0 0.0  .7
h2(basis=small) 0.0 0.0 -.7
$end
 
$small h
/ old H2
  type=s
 48.4479    1.00
  type=s
 7.28346    1.00
  type=s
 1.65139    1.00
  type=s
 .462447    1.00
  type=s
 .145885    1.00
  type=s
  1.00      1.00
  type=s
  .25       1.00
  type=s
  .0725     1.00
  type=s
  .03625    1.00
  type=s
  .018      1.00
  type=p
  1.5000    1.00
  type=p
  .50       1.00
  type=p
  .25       1.00
  type=p
  .125      1.00
  type=p
  .0625     1.00
  type=p
  .0312     1.00
  type=p
  .011      1.00
$end
 
$drop
 atom=1 type=p keep=(z)
 atom=2 type=p keep=(z)
 atom=3 type=p keep=(z)
$end

$target
1typ1;1 1typ2;1 1typ3;1 1typ4;1 30typ5;1
na=1 nb=1 ns=1
$end
 
$target-groups
numel=(2,0,0,1,0)
numel=(1,1,0,1,0) spin=(1,2,2,1,1)
numel=(1,0,1,1,0)
numel=(0,0,2,1,0)
numel=(0,2,0,1,0)
$end

$p-and-q
1typ1;1 1typ2;1 1typ3;1 3typ4;1 28typ5;1
na=1 nb=1 ns=1
$end
 
$p-and-q-groups
numel=(2,0,0,1,0)
numel=(1,1,0,1,0) spin=(1,2,2,1,1)
numel=(1,0,1,1,0)
numel=(0,0,2,1,0)
numel=(0,2,0,1,0)
numel=(1,2,0,0,0) qspace
$end
 
$kohndt
 point-buffer=5000 no-contracted=47 no-open-channels=2
 maximum-l-value=11 maximum-m-value=0
 maximum-r-value=100.  points-per-interval-for-spline-coefficients=10
 maximum-lm-chan=6 maximum-k**2=2.0
$end







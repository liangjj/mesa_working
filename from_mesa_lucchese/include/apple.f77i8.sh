export FC=f77
export LD=f77
export FC90=f90
export BUFLEN=1024
export MAXCOR=200000000
export MAXNBFKOHN=300
export MAXPRIMKOHN=500
export MAXCHAN=20
export MAXLTOP=10
export MAXLMTOP=40
export MAXSMALL=200
export MAXBIG=800       # maximum number of parital waves + scattering functions summed over channels

export LDFLAGS="-lU77 -i8 -m64 -X-map -Xmap.out"

export FFLAGS90="-i8 -m64"

export FFLAGSC="-c -g -s -f -N15 -i8 -m64"
export FFLAGSO="-c -g -O -s -f -N15 -i8 -m64"
export FFLAGSO2="-c -g -O2 -s -f -N15 -i8 -m64"

export FFLAGSOZ="-c -g -O -s -f -N15 -i8 -m64"
export FFLAGSOZB="-c -g -O -s -f -N15 -i8 -m64 -C"
export FFLAGSO2Z="-c -g -O2 -s -f -N15 -i8 -m64"
export FFLAGSO5Z="-c -g -O3 -s -f -N15 -i8 -m64"

export FFLAGSZ="-c -g -s -f -N15 -i8 -m64"
export FFLAGSZA="-c -g -s -f -N15 -i8 -m64"
export FFLAGSZB="-c -g -s -f -N15 -i8 -m64"

export BLASUSE=" -L/Users/lucchese/Applications/lapacki8/lapack-3.4.2 -lrefblas"
export BLASUSEB=$BLASUSE

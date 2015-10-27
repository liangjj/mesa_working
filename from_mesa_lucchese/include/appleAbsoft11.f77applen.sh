export FC=f77
export LD=f77
export FC90=f90
export BUFLEN=2048
export MAXCOR=120000000
export MAXNBFKOHN=300
export MAXPRIMKOHN=500
export MAXCHAN=20
export MAXLTOP=10
export MAXLMTOP=40
export MAXSMALL=200
export MAXBIG=800       # maximum number of parital waves + scattering functions summed over channels

export LDFLAGS="-lU77"

export FFLAGS90=

export FFLAGSC="-c -g -s -f -N15 "
export FFLAGSO="-c -g -O -s -f -N15 "
export FFLAGSO2="-c -g -O2 -s -f -N15"

export FFLAGSOZ="-c -g -O -s -f -N15"
export FFLAGSO2Z="-c -g -O2 -s -f -N15"
export FFLAGSO5Z="-c -g -O3 -s -f -N15"

export FFLAGSZ="-c -g -s -f -N15"
export FFLAGSZA="-c -g -s -f -N15"
export FFLAGSZB="-c -g -s -f -N15"

export BLASUSE=" -X -framework -X Accelerate"
export BLASUSEB=" -X -framework -X Accelerate"

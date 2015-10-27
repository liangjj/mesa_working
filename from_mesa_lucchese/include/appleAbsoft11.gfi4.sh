#
# works with gfortran ver 4.6.2 and the blas library compiled with the
# same -fdefault-integer-8 flag
#

export FC=gfortran
export LD=gfortran
export FC90=gfortran
export BUFLEN=2048
export MAXCOR=120000000
export MAXNBFKOHN=300
export MAXPRIMKOHN=500
export MAXCHAN=20
export MAXLTOP=10
export MAXLMTOP=40
export MAXSMALL=200
export MAXBIG=800       # maximum number of parital waves + scattering functions summed over channels

export LDFLAGS=

export FFLAGS90=

export FFLAGSC="-c -g -finit-local-zero"
export FFLAGSO="$FFLAGSC -O"
export FFLAGSO2="$FFLAGSC -O2"

export FFLAGSOZ="$FFLAGSC -O"
export FFLAGSO2Z="$FFLAGSC -O2"
export FFLAGSO5Z="$FFLAGSC -O3"

export FFLAGSZ="$FFLAGSC"
export FFLAGSZA="$FFLAGSC -fno-automatic"
export FFLAGSZB="$FFLAGSC -fno-backslash"

export BLASUSE="-L/Users/lucchese/Applications/lapackgfi4/lapack-3.4.2 -lrefblas"
export BLASUSEB=$BLASUSE

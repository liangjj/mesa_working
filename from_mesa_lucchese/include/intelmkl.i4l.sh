export FC=ifort
export LD=ifort
export FC90=ifort
export BUFLEN=2048
export MAXCOR=500000000
export MAXNBFKOHN=300
export MAXPRIMKOHN=500
export MAXCHAN=20
export MAXLTOP=10
export MAXLMTOP=40
export MAXSMALL=200
export MAXBIG=800       # maximum number of parital waves + scattering functions summed over channels

export LDFLAGS=

export FFLAGS90="-g"

#
# -fbounds-check for array bound checking
#

export FFLAGSC="-c -g -72 -save -zero -assume byterecl"
export FFLAGSO="$FFLAGSC -O"
export FFLAGSO2="$FFLAGSC -O2"

export FFLAGSOZ="$FFLAGSC -O"
export FFLAGSO2Z="$FFLAGSC -O2"
export FFLAGSO5Z="$FFLAGSC -O3"

export FFLAGSZ=$FFLAGSC
export FFLAGSZA=$FFLAGSC
export FFLAGSZB=$FFLAGSC

export BLASUSE="-L/home/r-lucchese/Applications/lapacki4/lapack-3.1.1 -lblas_LINUX"
export BLASUSEB=$BLASUSE
# to add link map -Xlinker -M

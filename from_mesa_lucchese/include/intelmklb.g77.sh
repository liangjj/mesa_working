export FC=g77
export LD=g77
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

export FFLAGS90="-g"

#
# -fbounds-check for array bound checking
#

export FFLAGSC="-c -g -Wno-globals"
export FFLAGSO="-c -g -Wno-globals -O"
export FFLAGSO2="-c -g -Wno-globals -O2"

export FFLAGSOZ="-c -g -Wno-globals -O -finit-local-zero"
export FFLAGSO2Z="-c -g -Wno-globals -O2 -finit-local-zero"
export FFLAGSO5Z="-c -g -Wno-globals -O5 -finit-local-zero"

export FFLAGSZ="-c -g -Wno-globals -finit-local-zero"
export FFLAGSZA="-c -g -Wno-globals -finit-local-zero -fno-automatic"
export FFLAGSZB="-c -g -Wno-globals -finit-local-zero -fno-backslash"

export BLASUSE="-L/scratch/rrl581a/lapack-3.1.1 -llapack -lblas"
export BLASUSEB=$BLASUSE
# to add link map -Xlinker -M

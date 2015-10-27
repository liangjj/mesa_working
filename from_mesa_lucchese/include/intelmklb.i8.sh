export FC=ifort
export LD=ifort
export FC90=ifort
export BUFLEN=1024

export MAXCOR=625000000 # size of balank common in mesa, in interger words (i4 or i8)
export MAXNBFKOHN=500   # maximum number of basis functions
export MAXPRIMKOHN=700  # maximum number of gaussian primitives
export MAXCHAN=30       # maximum number of coupled channels
export MAXLTOP=20       # maximum asymptotic L
export MAXLMTOP=60      # maximum number of partial waves (with different l and m) in any channel
export MAXSMALL=400     # maximum total number of partial waves summed over all channels, 
                        # i.e. size of the asymptotic scattering matrix
export MAXBIG=1000      # maximum number of parital waves + scattering functions summed over channels

export LDFLAGS="-g -72 -i8 -heap-arrays 10 -shared-intel -mcmodel=medium -save -zero -assume byterecl"

export FFLAGS90="-g -i8"

#
# -fbounds-check for array bound checking
#

export FFLAGSC="-c -g -72 -i8 -heap-arrays 10 -shared-intel -mcmodel=medium -save -zero -assume byterecl"
export FFLAGSO="$FFLAGSC -O"
export FFLAGSO2="$FFLAGSC -O2"

export FFLAGSOZ="$FFLAGSC -O"
export FFLAGSO2Z="$FFLAGSC -O2"
export FFLAGSO5Z="$FFLAGSC -O3"

export FFLAGSZ=$FFLAGSC
export FFLAGSZA=$FFLAGSC
export FFLAGSZB=$FFLAGSC

export BLASUSE="-L/scratch/rrl581a/lapacki8/lapack-3.4.2 -lrefblas"
export BLASUSEB=$BLASUSE
# to add link map -Xlinker -M

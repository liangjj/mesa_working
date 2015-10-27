#
# works with gfortran ver 4.6.2 and the blas library compiled with the
# same -fdefault-integer-8 flag
#

export FC=gfortran
export LD=gfortran
export FC90=gfortran
export BUFLEN=1024

export MAXCOR=625000000 # size of blank common in mesa, in integer words (i4 or i8)
export MAXNBFKOHN=300   # maximum number of basis functions
export MAXPRIMKOHN=500  # maximum number of gaussian primitives
export MAXCHAN=20       # maximum number of coupled channels
export MAXLTOP=10       # maximum asymptotic L
export MAXLMTOP=40      # maximum number of partial waves (with different l and m) in any channel
export MAXSMALL=200     # maximum total number of partial waves summed over all channels, 
                        # i.e. size of the asymptotic scattering matrix
export MAXBIG=3500      # maximum number of partial waves + scattering functions summed over channels

export LDFLAGS=

export FFLAGS90="-fdefault-integer-8"

export FFLAGSC="-c -g -fdefault-integer-8 -finit-local-zero"
export FFLAGSO="$FFLAGSC -O"
export FFLAGSO2="$FFLAGSC -O2"

export FFLAGSOZ="$FFLAGSC -O"
export FFLAGSO2Z="$FFLAGSC -O2"
export FFLAGSO5Z="$FFLAGSC -O3"

export FFLAGSZ="$FFLAGSC"
export FFLAGSZA="$FFLAGSC -fno-automatic"
export FFLAGSZB="$FFLAGSC -fno-backslash"

export BLASUSE="-L/Users/lucchese/Applications/lapackgfi8/lapack-3.4.2 -lrefblas"
export BLASUSEB=$BLASUSE

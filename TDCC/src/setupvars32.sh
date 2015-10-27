machine=$(uname -n)
#if [ "x$machine" == "xraritan.nist.gov" ]; then
if [ "x$PBS_ROOT" != "x" ]; then
    export PREFERRED_FORTRAN=pgi
    export PREFERRED_ARCH=32
    export PREFERRED_MPI=mpich
    . /usr/local/etc/mpi_fort_env_new.sh

    LAPACK95=~/soft/lib/lapack95
    LAPACK95_MODULES=${LAPACK95}/modules
    LAPACK95_LIB=${LAPACK95}/${PREFERRED_ARCH}/lapack95.a
    export FFLAGS="-I${LAPACK95_MODULES} -Mvect=sse -Mscalarsse -Mcache_align -O2"
    #export LDFLAGS="-llapack -lblas"
    export LDFLAGS="${LAPACK95_LIB} -lacml"
elif [ "x$machine" == "xmaster" ]; then # klaus bartschat's machine at drake
    export FFLAGS="-Mvect=sse -Mscalarsse -Mcache_align -O2 -tp=k8-64"
    export LDFLAGS="-llapack -lblas"
fi

#export PREFERRED_FORTRAN=intel
#export FFLAGS="-O3 -ip -no-prec-div -static -xN -unroll -ftz"


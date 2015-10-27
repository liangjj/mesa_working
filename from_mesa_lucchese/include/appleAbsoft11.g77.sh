export FC=g77
export LD=g77
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

export LDFLAGS=

export FFLAGS90=

export FFLAGSC="-c -g"
export FFLAGSO="-c -g -O"
export FFLAGSO2="-c -g -O2"

export FFLAGSOZ="-c -g -O -finit-local-zero"
export FFLAGSO2Z="-c -g -O2 -finit-local-zero"
export FFLAGSO5Z="-c -g -O5 -finit-local-zero"

export FFLAGSZ="-c -g -finit-local-zero"
export FFLAGSZA="-c -g -finit-local-zero -fno-automatic"
export FFLAGSZB="-c -g -finit-local-zero -fno-backslash"

export BLASUSE="/System/Library/Frameworks/vecLib.framework/versions/A/vecLib"
export BLASUSEB="/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate"

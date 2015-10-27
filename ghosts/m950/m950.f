*deck %W%  %G%
      program m950
c
c  **************kohn interface code*****************
c
c   5/1/90  -- changes for symmetry  bhl
c              nl2 read from route and nmotot = nl2 + nsmall
c              written to bndmat   
c
c   7/19/89 -- changes for multiple scattering energies and for
c               getting ci energies from mesa files
c
c   7/23/89 -- changes to read a print option ... bhl at llnl
c
c
c   makes input files for "quad" suite of kohn scattering codes
c
c  from prototype interface by byron lengsfield.  modified and
c  extended by cwm and tnr 7/5/89
c
c   complete set of output files specified as of 7/12/89
c
c    geobas = geometry and basis set
c    denmat = transition density matrices
c    bndmat  = transformation matrix
c              overlap matrix
c              direct hamiltonians for channel pairs
c              hpp-e
c              hopt(e)
c
      common // ia(1)
      dimension z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      integer wpadti
c
c
      call drum
      call getscm(need,z,maxcor,'m950',0)
      intoff=wpadti(ioff)
      call pm950(z(ioff),ia(intoff),maxcor)
      call chainx(0)
c
c
      stop
      end

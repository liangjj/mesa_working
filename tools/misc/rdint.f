      program rdint
      implicit integer(a-z)
c
c     read one and two electron integrals from a formatted file. 
      parameter (maxcor=1000000)
      character*80 title
      real*8 z
      common//z(maxcor)
c
c
      intout=50
      open(unit=intout,file='intout',access='sequential',
     $     form='unformatted',err=99,status='old')
c
c     read the title
      read(intout) title
      write(6,*) title
c
c     read the number of orbitals.
      read(intout) nbf
      write(6,*) 'number of orbitals:',nbf
      nnp=nbf*(nbf+1)/2
c
c     read one-electron integrals.
c     these come in as the lower triangle, i>=j.
      read(intout) (z(i),i=1,nnp)
      ij=0
      do 10 i=1,nbf
         do 9 j=1,i
            ij=ij+1
            write(6,*) i,j,z(ij)
    9    continue
   10 continue
c
c     read the two-electron integrals. these are stored as a triangle of triangles.
c     for each index ij, i>=j, there is a record which contains all the triangles
c     kl, k>=l.  Note this form "doubles" the length of the integral list.
c
      do 20 kl=1,nnp
         read(intout) (z(i),i=1,nnp)
         write(6,*) 'two-electron block:  kl=',kl
         ij=0
         do 19 i=1,nbf
            do 18 j=1,i
               ij=ij+1
               write(6,*) i,j,z(ij)
   18       continue
   19    continue
   20 continue
c
c
      close(unit=intout)
   99 continue
      stop
      end

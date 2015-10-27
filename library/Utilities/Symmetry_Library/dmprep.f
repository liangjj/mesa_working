*deck @(#)dmprep.f	5.1  11/6/94
      subroutine dmprep(natoms,iprmut)
      implicit integer(a-z)
c
c     this routine dumps the contents of /repcom/.
c
      real*8 symops,chrtbl
      character*4 repnam
      common/repcom/nsymop,nreps,lblrep(32),chrtbl(10,16),symops(9,10)
      common/reploc/nsave,ncnt,isave(10)
      common/repnam/repnam(32)
      common/io/in,iout
      dimension iprmut(natoms,1)
c
 1000 format('  nsymop =',i2,' nreps =',i2/'  orbital labels:',32a4)
 1010 format('  character table:')
 1020 format(1x,10f7.1)
 1030 format('  nr. of saved operations =',i2,' saved = ',10i3)
 1040 format(3(10x,3d15.4,/))
 1050 format('  nuclear permutations:',10i3)
c
c     call rtrace(6hdmprep,1)
      write(iout,1000) nsymop, nreps, (repnam(i),i=1,nreps)
      write(iout,1010)
      do 100 i = 1, nreps
 100     write(iout,1020) (chrtbl(j,i),j=1,nsymop)
         write(iout,1030) nsave,(isave(i),i=1,nsave)
         do 200 i=1,nsymop
            write(iout,1040)(symops(j,i),j=1,9)
 200        write(iout,1050)(iprmut(j,i),j=1,10)
            return
            end

*deck lstwrt.f
c***begin prologue     lstwrt
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           last write hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            last dump of non-zero elements of h
c***                   
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       lstwrt
      subroutine lstwrt(ibuf,rbuf,header,lenbuf,n)
      implicit integer (a-z)
      real*8 rbuf
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      character*(*) header
      common/io/inp, iout 
      call iosys('write integer '//header//
     #           ' to hamiltonian without rewinding',
     #             2*n,ibuf,0,' ')
      call iosys('write integer '//header//
     #           ' to hamiltonian without rewinding',
     #             wptoin(n),rbuf,0,' ')
      return
      end       



















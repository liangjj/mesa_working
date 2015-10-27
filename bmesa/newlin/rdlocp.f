c $Header$
*deck rdlocp.f
c***begin prologue     rdlocp
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdlocp, link 6201, local, potential
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            read in local channel potentials
c***                   in co-ordinate space.       
c***references         
c***routines called
c***end prologue       rdlocp
      subroutine rdlocp (vloc,nchan,ns,nr,nth,nph,title)
      implicit integer(a-z)
      dimension vloc(*), nr(ns), nth(ns), nph(ns)
      character*8 vloc
      character*(*) title
c     the potential is written and read in packed form.
c     it is assumed this is done for each radial shell
c     thus the number of words in the file is the sum for all
c     shells of the product of nr*nth*nph*nchan(nchan+1)/2
      call iosys ('rewind '//title//' on lamdat read-and-write',0,0,0,
     1            ' ')  
      call iosys ('get length of '//title//' on lamdat',size,0,0,' ')
      ntri=nchan*(nchan+1)/2
      len=0
      do 10 is=1,ns
         words=ntri*nr(is)*nth(is)*nph(is)
         call iosys ('read real '//title//' from lamdat without '//
     1               'rewinding',words,vloc(len+1),0,' ')
c        multiply the potential by 2.d0 which is necessary for the solution
c        of the differential equation later on in the calculation.
         call sscal(words,2.d0,vloc(len+1),1)
         len=len+words
   10 continue
      if (len.gt.size) then
          call lnkerr ('error in local potential file size')
      endif
      return
      end

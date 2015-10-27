*deck ntrial.f
c***begin prologue     ntrial
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           vectors, rearrange
c***author             schneider, barry (nsf)
c***purpose            make a new set of trial vectors for
c***                   the next set of roots.
c***description
c                      vec     = vectors holding the trials
c                      n       = vector length
c                      size    = size of small matrix diagonalized in current
c                                pass
c                      con     = number of converged vectors on the 
c                                current pass
c                      nin     = number of trial vectors used
c                      trial   = number of trial vectors available
c                      maxvec  = maximum number of allowable vectors
c
c***source             itlib
c***                   
c***references         
c
c***routines called    
c***end prologue       ntrial
      subroutine ntrial(vec,n,size,con,nin,trial,maxvec)
      implicit integer (a-z)
      real*8 vec
      character*4 itoc
      dimension vec(n,*)
      common/io/inp, iout
c
c     there are,
c
      navail = size - con
c
c     vectors from the diagonalization of the small matrix which are
c     rayleigh-ritz approximations to the next set of roots.
c
      n2cpy=min(navail,trial,maxvec)
      if(n2cpy.gt.0) then
         call copy(vec(1,con+1),vec,n*n2cpy)
      endif
c
c     if there is still room left, add the unused trial vectors, provided
c     there are any left.
c
      nout=n2cpy
      toadd=trial- nin
      if(toadd.gt.0) then
         nout=min(n2cpy+toadd,trial,maxvec)
         add=nout-n2cpy
         nout=n2cpy+add
         if(add.gt.n2cpy) then
            point=nin+1
            do 10 i=n2cpy+1,nout
               call iosys('read real "trial:'//itoc(point)
     1                                       //'" from rwf',
     2                     n,vec(1,i),0,' ')
               point=point+1
 10         continue   
         endif
      endif
      write(iout,1) n2cpy, add, nout
      if(nout.ne.0) then
         do 20 i=1,nout
            call iosys('write real "trial:'//itoc(i)
     1                                     //'" to rwf',
     2                  n,vec(1,i),0,' ')
 20      continue
      else
         write(iout,2)
         call lnkerr('quit:no trials')
      endif
c
c     update values of trial and nin
c
      trial=nout
      nin=nout
      return
 1    format(/,1x,'creating a new set of trial vectors',
     1       /,1x,'number of unconverged roots copied = ',i6,
     2       /,1x,'number of unused trials copied     = ',i6,
     3       /,1x,'total number of new trials         = ',i6)
 2    format(/,1x,'no trial vectors are left:quit')
      end



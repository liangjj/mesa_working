*deck @(#)schvec.f	1.3  8/8/91
      subroutine schvec(c,stri,s,t2,t3,num,nnp)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)schvec.f	1.3   8/8/91
c
c***purpose            to orthonormalize the scf vector using
c                      the condition csc=1
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit integer (a-z)
c
      real*8 c(num,num),stri(nnp),s(num,num),t2(num,num),t3(num,num)
      real*8 sdot
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      common /io/ inp,iout
c
c
      call trtosq(s,stri,num,nnp)
c
c     ----- loop over orbitals -----
c
      do 10 vec=1,num
c
c     ----- orthogonalize to previous vectors -----
c
         if(vec.eq.1) go to 6
         do 5 old=1,vec-1
            call saxpy(num,-sdot(num,t2(1,old),1,c(1,vec),1),c(1,old),
     $           1,c(1,vec),1)
    5    continue
 6       continue
c
c     ----- normalize -----
c
         call sgemv('n',num,num,one,s,num,c(1,vec),1,zero,t3,1)
         call sscal(num,1.0/sqrt(sdot(num,t3,1,c(1,vec),1)),c(1,vec),1,
     $        c(1,vec),1)
c
c     ----- put c(vec)s in t2 -----
c
c        call sgemv('n',num,num,one,s,num,c(1,vec),1,zero,t2(1,vec),1)
 10   continue
c
c     ----- check that we've done it correctly -----
c
      call ebc(t2,s,c,num,num,num)
      call ebtc(t3,c,t2,num,num,num)
      err=0
      do 20 i=1,num
         do 15 j=1,i-1
            if (abs(t3(i,j)).gt.1.0d-06) then
               write (iout,11) i,j,t3(i,j)
 11            format (' ***** error,',2i4,' are not orthogonal',e15.3)
               err=err+1
            end if
 15      continue
         if (abs(t3(i,i)-1.0d+00).gt.1.0d-06) then
            write (iout,16) i,t3(i,i)
 16         format (' ***** error,',i8,' is not normalized',e15.3)
            err=err+1
         end if
 20   continue
      if (err.gt.0) stop 2
c
c
      return
      end

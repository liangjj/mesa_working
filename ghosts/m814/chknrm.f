*deck @(#)chknrm.f	1.1  11/30/90
      subroutine chknrm(c,s,t1,t2,t3,num,nnp)
c
c***purpose: check the orthonormality of the scf vector.
c
c paul saxe                  21 august 1984                 lanl
c
      implicit integer (a-z)
c
      real*8 c(num,num),s(nnp),t1(num,num),t2(num,num),t3(num,num)
      real*8 thresh,maxerr
c
      common /io/     inp,ioutpt
c
      data thresh /1.0d-05/
c
      call trtosq(t1,s,num,nnp)
      call ebtc(t2,c,t1,num,num,num)
      call ebc(t3,t2,c,num,num,num)
c
c     ----- check that t3 is a unit matrix -----
c
      maxerr=0.0d+00
      err=0
      do 2 i=1,num
         do 1 j=1,num
            if (i.eq.j) then
               maxerr=max(maxerr,abs(t3(i,j)-1.0d+00))
               if (abs(t3(i,j)-1.0d+00).gt.thresh) then
                  write (ioutpt,90) i,t3(i,i)
   90             format (' orbital ',i3,' is not normalized',f15.9)
                  err=err+1
               end if
            else
               maxerr=max(maxerr,abs(t3(i,j)))
               if (abs(t3(i,j)).gt.thresh) then
                  write (ioutpt,91) i,j,t3(i,j)
   91             format (' orbitals ',i3,' and ',i3,
     #                    ' are not orthogonal:',f15.9)
                  err=err+1
               end if
            end if
    1    continue
    2 continue
c
      if (err.gt.0) then
         write (ioutpt,92) err
   92    format (//,' ##### chknrm: the scf vector is not normalized:',
     #           i5,//)
         stop 92
      end if
c
      write (ioutpt,94) maxerr
   94 format (' maximum deviation from orthonormality:',e14.2)
c
c     ----- schmidt orthonormalize the vectors -----
c
      maxerr=1.0d-08
      call schmdt(c,s,t1,t2,t3,num,num,nnp,maxerr)
c
c
      return
      end

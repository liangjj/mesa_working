*deck @(#)chknrm.f	5.1  11/6/94
      subroutine chknrm(c,s,t1,t2,t3,num,nnp)
c
c***begin prologue     chknrm
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)chknrm.f	5.1   11/6/94
c
c***purpose            check the orthonormality of the scf vector
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       chknrm
c
      implicit integer (a-z)
c
      real*8 c(num,num),s(nnp),t1(num,num),t2(num,num),t3(num,num)
      real*8 thresh
c
      common /io/ inp,iout
c
      data thresh /1.0d-07/
      save thresh
c
      call trtosq(t1,s,num,nnp)
      call ebtc(t2,c,t1,num,num,num)
      call ebc(t3,t2,c,num,num,num)
c
c     ----- check that t3 is a unit matrix -----
c
      err=0
      do 2 i=1,num
         do 1 j=1,num
            if (i.eq.j) then
               if (abs(t3(i,j)-1.0d+00).gt.thresh) then
                  write (iout,90) i,t3(i,i)
 90               format (' orbital ',i3,' is not normalized',f15.9)
                  err=err+1
               end if
            else
               if (abs(t3(i,j)).gt.thresh) then
                  write (iout,91) i,j,t3(i,j)
 91               format (' orbitals ',i3,' and ',i3,
     $                 'are not orthogonal:',f15.9)
                  err=err+1
               end if
            end if
    1    continue
    2 continue
c
      if (err.gt.0) then
         write (iout,92) err
 92      format (//,' ##### chknrm: the scf vector is not normalized:',
     $        i5,//)
         stop 92
      end if
c
c
      return
      end

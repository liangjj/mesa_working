*deck @(#)schmit.f	5.1  11/6/94
      subroutine schmit(b,mdim,iter,mvec,mmvec,s,eps,ierr)
      implicit real*8(a-h,o-z)
c
      real*8 b(mdim,*),s(*)
      common /io/ inp,iout
c
      epss=.001d0*eps
c
c     write(iout,*)' schmit: input vectors ',iter,mvec
c     call matout(b,mdim,iter+mvec,mdim,iter+mvec,iout)
c
      mmvec=mvec
      ix=0
c
      do 1 i=1,mvec
         ix=ix+1
         xx=sqrt(sdot(mdim,b(1,iter+ix),1,b(1,iter+ix),1))
         if(xx.lt.epss) then
            ileft=mvec-i
            mmvec=mmvec-1
            if(ileft.ne.0) then
               call scopy(ileft*mdim,b(1,iter+ix+1),1,b(1,iter+ix),1)
            endif
            ix=ix-1
         else
            call sscal(mdim,1.d0/xx,b(1,iter+ix),1)
         end if
  1   continue
c
c     ------ calculate the overlap matrix in s -----
c
      ix=iter
      mvec=mmvec
c
      if(ix.ne.0) then
         do 3 i=1,mvec
            jx=ix
            ix=ix+1
            call ebtc(s,b,b(1,ix),jx,mdim,1)
            do 2 j=1,jx
               s(j)=-s(j)
  2         continue
            call smxpy(mdim,b(1,ix),jx,mdim,s,b)
            xx=sqrt(sdot(mdim,b(1,ix),1,b(1,ix),1))
            if(xx.lt.eps) then
               mmvec=mmvec-1
               ileft=mvec-i
               if(ileft.ne.0) then
                  call scopy(ileft*mdim,b(1,iter+ix+1),1,b(1,iter+ix),1)
               endif
               ix=ix-1
            else
               call sscal(mdim,1.d0/xx,b(1,ix),1)
            end if
   3     continue
      else
         do 6 i=1,mvec
            jx=ix
            ix=ix+1
c
            if(jx.ne.0) then
               call ebtc(s,b,b(1,ix),jx,mdim,1)
               do 5 j=1,jx
                  s(j)=-s(j)
  5            continue
               call smxpy(mdim,b(1,ix),jx,mdim,s,b)
            end if
c
            xx=sqrt(sdot(mdim,b(1,ix),1,b(1,ix),1))
            if(xx.lt.eps) then
               mmvec=mmvec-1
               ileft=mvec-i
               if(ileft.ne.0) then
                  call scopy(ileft*mdim,b(1,ix),1,b(1,ix+1),1)
               end if
               ix=ix-1
            else
               call sscal(mdim,1.d0/xx,b(1,ix),1)
            end if
   6     continue
      end if
c
      ierr=0
      if(mmvec.eq.0) then
         ierr=1
      end if
c
      return
      end

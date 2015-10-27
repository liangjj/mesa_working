*deck @(#)redund.f	5.1  11/6/94
      subroutine redund(symcen,angmom,dercen,dermom,npass,flip)
c***begin prologue     redund.f
c***date written       860101   (yymmdd)
c***revision date      11/6/94 
c  20 november 1987    pws at lanl
c    adding a variable 'flip' to keep track of how we've permuted the
c    indices.
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)redund.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c  'flip'        indices
c -------       ---------
c    1            ijkl
c    2            jikl
c    3            ijlk
c    4            jilk
c    5            klij
c    6            klji
c    7            lkij
c    8            lkji
c
c***references
c
c***routines called    (none)
c
c***end prologue       redund.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
      integer symcen(4),angmom(4)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer dercen(4),dermom(4)
c     --- output variables ---
      integer npass,flip
c     --- scratch arrays ---
c     --- local variables ---
      logical iandj,iandk,iandl,jandk,jandl,kandl,ikajl
c
      iandj=symcen(1).eq.symcen(2)
      iandk=symcen(1).eq.symcen(3)
      iandl=symcen(1).eq.symcen(4)
      jandk=symcen(2).eq.symcen(3)
      jandl=symcen(2).eq.symcen(4)
      kandl=symcen(3).eq.symcen(4)
      ikajl=iandk.and.jandl
c
      if (iandl.and.iandk.and.iandj) then
c
c        ----- [ii;ii] -----
         npass=0
         dercen(1)=symcen(1)
         dercen(2)=symcen(1)
         dercen(3)=symcen(1)
         dercen(4)=symcen(1)
         dermom(1)=angmom(1)
         dermom(2)=angmom(2)
         dermom(3)=angmom(3)
         dermom(4)=angmom(4)
         flip=1
      else if (iandk) then
         if (iandl) then
c
c           ----- [ij;ii] -----
            npass=4
            dercen(1)=symcen(2)
            dercen(2)=symcen(1)
            dercen(3)=symcen(1)
            dercen(4)=symcen(1)
            dermom(1)=angmom(2)
            dermom(2)=angmom(1)
            dermom(3)=angmom(3)
            dermom(4)=angmom(4)
            flip=2
         else if (iandj) then
c
c           ----- [ii;ij] -----
            npass=4
            dercen(1)=symcen(4)
            dercen(2)=symcen(1)
            dercen(3)=symcen(1)
            dercen(4)=symcen(1)
            dermom(1)=angmom(4)
            dermom(2)=angmom(3)
            dermom(3)=angmom(2)
            dermom(4)=angmom(1)
            flip=8
         else if (jandl) then
c
c           ----- [ij;ij] -----
            npass=3
            dercen(1)=symcen(1)
            dercen(2)=symcen(2)
            dercen(3)=symcen(1)
            dercen(4)=symcen(2)
            dermom(1)=angmom(1)
            dermom(2)=angmom(2)
            dermom(3)=angmom(3)
            dermom(4)=angmom(4)
            flip=1
         else
c
c           ----- [ij;ik] -----
            npass=3
            dercen(1)=symcen(2)
            dercen(2)=symcen(1)
            dercen(3)=symcen(4)
            dercen(4)=symcen(3)
            dermom(1)=angmom(2)
            dermom(2)=angmom(1)
            dermom(3)=angmom(4)
            dermom(4)=angmom(3)
            flip=4
         end if
      else if (iandj) then
         if (kandl) then
c
c           ----- [ii;jj] -----
            npass=2
            dercen(1)=symcen(1)
            dercen(2)=symcen(1)
            dercen(3)=symcen(3)
            dercen(4)=symcen(3)
            dermom(1)=angmom(1)
            dermom(2)=angmom(2)
            dermom(3)=angmom(3)
            dermom(4)=angmom(4)
            flip=1
         else
c
c           ----- [ii;jk] -----
            npass=2
            dercen(1)=symcen(3)
            dercen(2)=symcen(4)
            dercen(3)=symcen(1)
            dercen(4)=symcen(1)
            dermom(1)=angmom(3)
            dermom(2)=angmom(4)
            dermom(3)=angmom(2)
            dermom(4)=angmom(1)
            flip=6
         end if
      else
         if (jandl) then
            if (jandk) then
c
c              ----- [ij;jj] -----
               npass=4
               dercen(1)=symcen(1)
               dercen(2)=symcen(2)
               dercen(3)=symcen(2)
               dercen(4)=symcen(2)
               dermom(1)=angmom(1)
               dermom(2)=angmom(2)
               dermom(3)=angmom(3)
               dermom(4)=angmom(4)
               flip=1
            else
c
c              ----- [ik;jk] -----
               npass=3
               dercen(1)=symcen(1)
               dercen(2)=symcen(2)
               dercen(3)=symcen(3)
               dercen(4)=symcen(4)
               dermom(1)=angmom(1)
               dermom(2)=angmom(2)
               dermom(3)=angmom(3)
               dermom(4)=angmom(4)
               flip=1
            end if
         else if (jandk) then
c
c           ----- [ij;jk] -----
            npass=3
            dercen(1)=symcen(1)
            dercen(2)=symcen(2)
            dercen(3)=symcen(4)
            dercen(4)=symcen(3)
            dermom(1)=angmom(1)
            dermom(2)=angmom(2)
            dermom(3)=angmom(4)
            dermom(4)=angmom(3)
            flip=3
         else if (kandl) then
c
c           ----- [ij;kk] -----
            npass=2
            dercen(1)=symcen(1)
            dercen(2)=symcen(2)
            dercen(3)=symcen(3)
            dercen(4)=symcen(3)
            dermom(1)=angmom(1)
            dermom(2)=angmom(2)
            dermom(3)=angmom(3)
            dermom(4)=angmom(4)
            flip=1
         else
c
c           ----- [ij;kl] -----
            npass=1
            dercen(1)=symcen(1)
            dercen(2)=symcen(2)
            dercen(3)=symcen(3)
            dercen(4)=symcen(4)
            dermom(1)=angmom(1)
            dermom(2)=angmom(2)
            dermom(3)=angmom(3)
            dermom(4)=angmom(4)
            flip=1
         end if
      end if
c
c
      return
      end

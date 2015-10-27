*deck @@(#)srone.f	5.1  11/6/94
      subroutine srone(ni,is,bftorb,h,nii,ijadd,numij,kadd,ladd,
     #                 norbs,val,lab,bin,lenbin,nbin,isort,rsort,
     #                 lnsort,ops)
c
      implicit integer (a-z)
c
      real*8 h(nii),val(lenbin),rsort(*)
      integer isort(lnsort)
      integer bftorb(ni),ijadd(numij),kadd(norbs),ladd(norbs)
      integer lab(lenbin),bin(lenbin)
      logical logkey
      character*(*) ops
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
c
c     ----- check print option -----
      if (logkey(ops,'print=m820=one-electron-integrals',.false.,' '))
     $     then
         write (iout,9001)
 9001    format (5x,'m820: one-electron integrals (mo)')
         call print(h,nii,ni,iout)
      end if
c
c     ----- disperse the one-electron integrals -----
c
      ij=0
      do 20 i=1,ni
         iorb=bftorb(i)
         if (nbin+i.gt.lenbin) then
            call sorter('with bin',isort,rsort,0,nbin,lab,
     #                   bin,val,0,0,0,0)
            nbin=0
         end if
         do 10 j=1,i
            ij=ij+1
            jorb=bftorb(j)
            ijorb=ioff(max(iorb,jorb),max(iorb,jorb))
            nbin=nbin+1
            val(nbin)=h(ij)
            if (j.ne.i) then
               lab(nbin)=ijadd(ijorb)+kadd(max(iorb,jorb))+
     #                                ladd(min(iorb,jorb))+3
            else
               lab(nbin)=ijadd(ijorb)+kadd(iorb)+ladd(iorb)+2
            end if
   10    continue
   20 continue
c
c
      return
      end

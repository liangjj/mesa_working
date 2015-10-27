*deck  @(#)hqqsort.f	5.1 11/6/94
      subroutine hqqsort(bptr,bcntr,xbin,ibin,hbuf,ibuf,
     $                   nbins,binsiz,lnbuf,lnblk,blksiz,lsize,
     $                   start,iu,binsz2,nspin,hqq,npdim,twalks,
     #                   filtyp)
      implicit integer(a-z)
      integer bptr(nbins),bcntr(nbins),ibin(binsz2,nbins)
      integer ibuf(2,lnbuf)
      real*8 xbin(binsiz,nbins),hbuf(lnbuf),enuc,hqq(npdim,npdim)
      character*(*) filtyp
      common /io/ inp,iout
c
c     25 july 1991    rlm at lanl
c     working on 32 bit version. 
c
      next=start
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,enuc,0,' ')
c
       call rzero(hqq,npdim*npdim)
      do 1 i=1,nbins
         bptr(i)=0
         bcntr(i)=0
  1   continue
c
      call iosys('read integer "number of elements" from hamiltonian',
     $            1,ntotal,0,' ')
      call iosys('read integer "ci buffer length" from rwf',
     $            1,lnbuf,0,' ')
c
      npass=ntotal/lnbuf
      left=ntotal - lnbuf*npass
      last=lnbuf
      if(left.gt.0) then
         npass=npass+1 
         last=left
      endif
      nleft=ntotal
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      call iosys('length of buffers on hamiltonian',nlen,0,0,' ')
c
      do 101 mm=1,npass
         n=lnbuf
         if(mm.eq.npass) then
            n=last
         endif
         call iosys('read integer buffers from hamiltonian'//
     $              ' without rewinding',2*n,ibuf,0,' ')
         call iosys('read integer buffers from hamiltonian'//
     $              ' without rewinding',wptoin(n),hbuf,0,' ')
c
c         write(iout,*) npass, n, lnbuf, nleft
         do 2 i=1,min(lnbuf,nleft)
            ip=max(ibuf(1,i),ibuf(2,i))
            jp=min(ibuf(1,i),ibuf(2,i))
            if(ip.le.npdim .and. jp.le.npdim) then
               hqq(ip,jp)=hqq(ip,jp)+hbuf(i)
            end if
            ib=(ip-1)/lsize+1
            jb=(jp-1)/lsize+1
            io=ip-(ib-1)*lsize
            jo=jp-(jb-1)*lsize
            numblk=ib*(ib-1)/2+jb
            index=(jo-1)*lsize+io
            addres=(numblk-1)*lnblk+index
            nbin=(addres-1)/blksiz+1
            bindex=addres-(nbin-1)*blksiz
            count=bcntr(nbin)+1
            xbin(count,nbin)=hbuf(i)
            ibin(count,nbin)=bindex
            if(count.lt.binsiz) then
               bcntr(nbin)=bcntr(nbin)+1
            else
               ibin(binsiz+1,nbin)=binsiz
               ibin(binsiz+2,nbin)=bptr(nbin)
               bcntr(nbin)=0
               bptr(nbin)=next
               iw=next
               call putham(iu,xbin(1,nbin),wptoin(binsiz),iw,filtyp)
               iw=iw+wptoin(binsiz)
               call putham(iu,ibin(1,nbin),binsz2,iw,filtyp)
               next=iw+binsz2
            end if
c
c -- square diagonal block
c
            if(ib.eq.jb) then
               if(io.ne.jo) then
                  index=(io-1)*lsize+jo
                  addres=(numblk-1)*lnblk+index
                  nbin=(addres-1)/blksiz+1
                  bindex=addres-(nbin-1)*blksiz
                  count=bcntr(nbin)+1
                  xbin(count,nbin)=hbuf(i)
                  ibin(count,nbin)=bindex
                  if(count.lt.binsiz) then
                     bcntr(nbin)=bcntr(nbin)+1
                  else
                     ibin(binsiz+1,nbin)=binsiz
                     ibin(binsiz+2,nbin)=bptr(nbin)
                     bcntr(nbin)=0
                     bptr(nbin)=next
                     iw=next
                     call putham(iu,xbin(1,nbin),wptoin(binsiz),iw,
     #                           filtyp)
                     iw=iw+wptoin(binsiz)
                     call putham(iu,ibin(1,nbin),binsz2,iw,filtyp)
                     next=iw+binsz2
                  end if
               end if
            end if
c
c -- end square
c
  2      continue
         nleft=nleft-lnbuf
 101  continue
c
      call iosys('read real diagonals from hamiltonian',
     $           nspin,hbuf,0,' ')
c
      do 102 i=1,nspin
         ip=i
         jp=i
         if(i.le.npdim) then
            hqq(i,i)=hqq(i,i)+hbuf(i)
         end if
         ib=(ip-1)/lsize+1
         jb=(jp-1)/lsize+1
         io=ip-(ib-1)*lsize
         jo=jp-(jb-1)*lsize
         numblk=ib*(ib-1)/2+jb
         index=(jo-1)*lsize+io
         addres=(numblk-1)*lnblk+index
         nbin=(addres-1)/blksiz+1
         bindex=addres-(nbin-1)*blksiz
         count=bcntr(nbin)+1
         xbin(count,nbin)=hbuf(i)+enuc
         ibin(count,nbin)=bindex
         if(count.lt.binsiz) then
            bcntr(nbin)=bcntr(nbin)+1
         else
            ibin(binsiz+1,nbin)=binsiz
            ibin(binsiz+2,nbin)=bptr(nbin)
            bcntr(nbin)=0
            bptr(nbin)=next
            iw=next
            call putham(iu,xbin(1,nbin),wptoin(binsiz),iw,filtyp)
            iw=iw+wptoin(binsiz)
            call putham(iu,ibin(1,nbin),binsz2,iw,filtyp)
            next=iw+binsz2
         end if
  102 continue
c
c -- dump the last bin to disk
c
      do 3 i=1,nbins
         if(bcntr(i).ne.0) then
            ibin(binsiz+1,i)=bcntr(i)
            ibin(binsiz+2,i)=bptr(i)
            bcntr(i)=0
            bptr(i)=next
            iw=next
            call putham(iu,xbin(1,i),wptoin(binsiz),iw,filtyp)
            iw=iw+wptoin(binsiz)
            call putham(iu,ibin(1,i),binsz2,iw,filtyp)
            next=iw+binsz2
         end if
  3   continue
c
c     put the pointers at the beginning of the file
      call putham(iu,bptr,nbins,0,filtyp)
c
      do 1001 i=1,npdim
         do 1002 j=1,i
            hqq(j,i)=hqq(i,j)
 1002    continue
         hqq(i,i)=hqq(i,i)+enuc
 1001 continue
c
c
      return
      end

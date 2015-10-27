*deck @(#)twoint.f	1.1  8/1/91
      subroutine twoint(label,int,lenbuf,values,ptr,lenbin,sort,
     #                  asort,lnsort,n2int,ijklpt,nnqsym,nso,nsym,
     #                  itape,bins,intfil,bfsym,bfnum)
c
      implicit integer (a-z)
c
      real*8 int(lenbuf),values(lenbin),sort(*)
      integer label(lenbuf),ptr(lenbin),ijklpt(nnqsym),nso(nsym)
      integer asort(lnsort),bins(lenbin)
      integer bfsym(lenbuf),bfnum(lenbuf)
      character*(*) intfil
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
      ioff2(i,j,k,l)=ioff(ioff(i,j),ioff(k,l))
c
c     ----- initialise the sorting routines -----
      write(iout,77)
   77 format(/,10x,'two-electron integral set')
c
      call sorter('start',asort,sort,lnsort,n2int,0,0,0,0,
     #             'sorted so integrals','ints',.true.)
c
      n=0
      nints=0
      read(14) kk
    5 continue
         read(14) nbuf,ilsti,label,int
         if(nbuf.le.0) go to 50
         do 45 i=1,nbuf
            call unpack(label(i),io,jo,ko,lo,mo)
            ism=bfsym(io)
            jsm=bfsym(jo)
            ksm=bfsym(ko)
            lsm=bfsym(lo)
            ior=bfnum(io)
            jor=bfnum(jo)
            kor=bfnum(ko)
            lor=bfnum(lo)
c
c            write (iout,230) i,ism,jsm,ksm,lsm,ior,jor,kor,lor,int(i)
c 230        format (i5,4i3,5x,4i3,f12.6)
c
c           ----- flip symmetries to be canonical -----
c
            if (jsm.gt.ism) then
               junk=ism
               ism=jsm
               jsm=junk
               junk=ior
               ior=jor
               jor=junk
            end if
            if (lsm.gt.ksm) then
               junk=ksm
               ksm=lsm
               lsm=junk
               junk=kor
               kor=lor
               lor=junk
            end if
            if (ksm.gt.ism.or.(ksm.eq.ism.and.lsm.gt.jsm)) then
               junk=ism
               ism=ksm
               ksm=junk
               junk=ior
               ior=kor
               kor=junk
               junk=jsm
               jsm=lsm
               lsm=junk
               junk=jor
               jor=lor
               lor=junk
            end if
c
c           ----- check for non-canonical symmetry labels -----
c
            if (jsm.gt.ism.or.lsm.gt.ksm.or.(ism.eq.ksm.and.lsm.gt.jsm))
     #                                                              then
               call lnkerr('non-canonical symmetry types')
            end if
c
c           ----- and non-canonical orbital labels -----
c
            if ((ism.eq.jsm.and.jor.gt.ior).or.
     #          (ksm.eq.lsm.and.lor.gt.kor).or.
     #          (ism.eq.ksm.and.jsm.eq.lsm.and.ior.eq.kor
     #           .and.lor.gt.jor)) then
               call lnkerr('non-canonical orbital types')
            end if
c
            if (n+2.gt.lenbin) then
c
c              ----- dump this bin to the sort routines -----
c
               call sorter('with bin',asort,sort,0,n,ptr,bins,values,
     #                      0,0,0,0)
               nints=nints+n
               n=0
            end if
c
c           ----- some integrals must be placed two places -----
c
            if (ism.eq.ksm.and.jsm.eq.lsm) then
c
c              ----- ij;kl --> ij;kl & kl;ij
c
               ijkl=ijklpt(ioff2(ism,jsm,ksm,lsm))
               values(n+1)=int(i)
               values(n+2)=int(i)
c
               if (ism.eq.jsm) then
c
c                 ----- [ii;ii] -----
c
                  ptr(n+1)=ijkl+
     #                     (ioff(kor,lor)-1)*ioff(nso(ism),nso(jsm))+
     #                     ioff(ior,jor)
                  ptr(n+2)=ijkl+
     #                     (ioff(ior,jor)-1)*ioff(nso(ksm),nso(lsm))+
     #                     ioff(kor,lor)
c
                  if (ptr(n+1).gt.n2int.or.ptr(n+2).gt.n2int) then
                     write (iout,1020) n,ptr(n+1),ptr(n+2),n2int,
     $                    ism,jsm,ksm,lsm,ior,jor,kor,lor,ijkl,
     $                    nso(ism),nso(jsm),nso(ksm),nso(lsm),
     $                    ioff(kor,lor)-1,
     $                    ioff(nso(ism),nso(jsm)),ioff(ior,jor),
     $                    ioff(nso(ksm),nso(lsm))
 1020                format ('iiii',i5,3i9,/,1x,4i3,5x,4i3,/,1x,i19,
     $                    4x,4i4,4x,4i9)
                     call lnkerr('iiii')
                  end if
               else 
c
c                 ----- [ij;ij] -----
c
                  ptr(n+1)=ijkl+ior+nso(ism)*(jor-1+nso(jsm)*(kor-1+
     #                     nso(ksm)*(lor-1)))
                  ptr(n+2)=ijkl+kor+nso(ksm)*(lor-1+nso(lsm)*(ior-1+
     #                     nso(ism)*(jor-1)))
                  if (ptr(n+1).gt.n2int.or.ptr(n+2).gt.n2int) then
                     write (iout,1021) n,ptr(n+1),ptr(n+2),n2int,
     $                    ism,jsm,ksm,lsm,ior,jor,kor,lor,ijkl,
     $                    nso(ism),nso(jsm),nso(ksm),nso(lsm)
 1021                format ('ijij',i5,3i9,/,1x,4i3,5x,4i3,/,1x,i19,
     $                    4x,4i4,4x,4i9)
                     call lnkerr('ijij')
                  end if
               end if
c
               n=n+2
            else
c
c              ----- ij;kl --> ij;kl
c
               ijkl=ijklpt(ioff2(ism,jsm,ksm,lsm))
               values(n+1)=int(i)
               if (ism.eq.jsm.and.ksm.eq.lsm) then
c
c                 ----- [ii;jj] -----
c
                  ptr(n+1)=ijkl+ioff(ior,jor)+ioff(nso(ism),nso(ism))*
     #                          (ioff(kor,lor)-1)
                  if (ptr(n+1).gt.n2int) then
                     write (iout,1022) n,ptr(n+1),ptr(n+2),n2int,
     $                    ism,jsm,ksm,lsm,ior,jor,kor,lor,ijkl,
     $                    nso(ism),nso(jsm),nso(ksm),nso(lsm),
     $                    ioff(kor,lor)-1,
     $                    ioff(nso(ism),nso(jsm)),ioff(ior,jor),
     $                    ioff(nso(ksm),nso(lsm))
 1022                format ('iijj',i5,3i9,/,1x,4i3,5x,4i3,/,1x,i19,
     $                    4x,4i4,4x,4i9)
                     call lnkerr('iijj')
                  end if
               else
c                  
c                 ----- [ij;kl] -----
c
                  ptr(n+1)=ijkl+ior+nso(ism)*(jor-1+nso(jsm)*(kor-1+
     #                     nso(ksm)*(lor-1)))
                  if (ptr(n+1).gt.n2int) then
                     write (iout,1023) n,ptr(n+1),ptr(n+2),n2int,
     $                    ism,jsm,ksm,lsm,ior,jor,kor,lor,ijkl,
     $                    nso(ism),nso(jsm),nso(ksm),nso(lsm)
 1023                format ('ijkl',i5,3i9,/,1x,4i3,5x,4i3,/,1x,i19,
     $                    4x,4i4,4x,4i9)
                     call lnkerr('ijkl')
                  end if
               end if
               n=n+1
            end if
 45      continue
 50   if(ilsti.eq.0) go to 5
c
c     ----- finish up the sort -----
c
      call sorter('with bin',asort,sort,0,n,ptr,bins,values,0,0,0,0)
      call sorter('end',asort,sort,0,0,0,0,0,0,0,0,0)
c
      nints=nints+n
      write (iout,60) nints
   60 format (/,' number of integrals processed out: ',i10)
c
c
      return
      end

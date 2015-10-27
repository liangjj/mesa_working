*deck %W%  %G%
      subroutine dm(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,
     #              levpt,levnr,rowoff,
     #              iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,ia,ib,
     #              ilevpt,ilevnr,inrows,
     #              nrows,nlevs,norbs,
     #              acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #              levdiv,nnp,c,nwks,h,g,bcoef,
     #              prtflg,nsym,numrow,drtpt,
     #              z,sval,isval,maxwks,ijpt,nnp0,ngsym,numsym,
     #              offsym,orbsym,h1,g1,ops,dunit)
c
      implicit integer (a-z)
c
c     ----- external unmodified arrays -----
c
      real*8 h1(nnp),g1(nnp,nnp)
      real*8 h(nnp0),g(ngsym)
      integer numsym(nsym),offsym(nsym),ijpt(nnp),orbsym(norbs)
      integer uparc(4,nrows),dnarc(4,nrows),upwt(4,nrows),dnwt(4,nrows)
      integer upnwks(nrows),dnnwks(nrows),a(nrows),b(nrows)
      integer levpt(nlevs),levnr(nlevs),rowoff(nrows)
      integer iuparc(4,inrows),idnarc(4,inrows),iupwt(4,inrows)
      integer idnwt(4,inrows)
      integer iupnwk(inrows),idnnwk(inrows),ia(inrows),ib(inrows)
      integer ilevpt(nlevs,nsym),ilevnr(nlevs,nsym)
      integer numrow(nsym),drtpt(nsym)
      integer sval(nrows),isval(inrows)
c
c     ----- external scratch arrays -----
c
      real*8 c(nwks)
      real*8 acoef(nlevs),bcoef(nlevs)
      real*8 z(*)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
c
c     ----- external unmodified scalars -----
c
      character*(*) ops
      character*16 dunit
      character*8 prtflg
      integer nsym
c
c     ----- local arrays -----
c
      real*8 coeffs(20,21)
c
c     ----- local scalars -----
c
      integer nn,nntri,eigv,ierr
      real*8 t,trace
c
c     ----- external functions -----
      logical logkey
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
c     ----- build the coefficient table -----
c
      call rzero(coeffs,20*21)
c
      do 700 i=3,20
         t = float(i-2)
         coeffs(i,1) = sqrt(t/(t+1.0d+00))
         coeffs(i,2) = -coeffs(i,1)
         coeffs(i,3) = coeffs(i,1)/sqrt(2.0d+00)
         coeffs(i,4) = -coeffs(i,3)
         coeffs(i,5) = sqrt((t+1.0d+00)/t)
         coeffs(i,6) = -coeffs(i,5)
         coeffs(i,7) = coeffs(i,5)/sqrt(2.0d+00)
         coeffs(i,8) = -coeffs(i,7)
         coeffs(i,9) = sqrt((t+2.0d+00)/(t*2.0d+00))
         coeffs(i,10) = -coeffs(i,9)
         coeffs(i,11) = sqrt(t/(2.0d+00*(t+2.0d+00)))
         coeffs(i,12) = -coeffs(i,11)
         coeffs(i,13) = sqrt(2.0d+00/(t*(t+1.0d+00)))
         coeffs(i,14) = sqrt(t*(t+2.0d+00))/(t+1.0d+00)
         coeffs(i,15) = -sqrt(t*(t+2.0d+00))/(t+1.0d+00)
         coeffs(i,16) = sqrt((t-1.0d+00)*(t+2.0d+00)/(t*(t+1.0d+00)))
         coeffs(i,17) = -coeffs(i,16)
         coeffs(i,18)=-sqrt(2.0d+00/(t*(t+2.0d+00)))
         coeffs(i,19) = 1.0d+00/t
         coeffs(i,20) = -1.0d+00/t
         coeffs(i,21) = -sqrt(2.0d+00)/t
  700 continue
c
c     ----- get the main and intermediate drt's -----
c
      call gdrt(nrows,nlevs,levpt,levnr,dnarc,uparc,dnnwks,
     #            upnwks,dnwt,upwt,a,b,sval,
     #            inrows,numrow,drtpt,nsym,ilevpt,ilevnr,idnarc,
     #            iuparc,idnnwk,iupnwk,idnwt,iupwt,ia,ib,isval)
c
c     ----- create the offset pointers for the rows on the dividing level
c
      n=0
      do 1 row=levpt(levdiv)+1,levpt(levdiv)+levnr(levdiv)
         rowoff(row)=n
         n=n+upnwks(row)*dnnwks(row)
    1 continue
      if (n.ne.dnnwks(1)) then
         call lnkerr('did not get nwks while creating offsets')
      end if
c
      n=0
      do 3 sym=1,nsym
         do 2 row=ilevpt(levdiv,sym)+1,
     #                ilevpt(levdiv,sym)+ilevnr(levdiv,sym)
            pt=drtpt(sym)-1+row
            n=max(n,iupnwk(pt)*idnnwk(pt))
    2   continue
    3 continue
c
c     ----- finish core allocation -----
c
c
      if (maxwks.le.0) then
         maxwks=n
      else
         maxwks=min(n,maxwks)
      end if
c
      if (prtflg.ne.'minimum') then
         write (iout,920) maxwks,n
 920     format (/,' m903:',/,5x,'attempting to use ',i7,
     #           ' intermediates at a time',/,
     #           5x,'of a possible     ',i7)
      end if
c
      call getscm(0,z,canget,'m903 available',0)
      top=min(canget,wpadti(nnp*maxwks+10000))
      call getscm(top,z,maxcor,'m903 ec',0)
c
c
c     ----- form the density matrices -----
c
      call fmdm(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,
     #     a,b,sval,levpt,levnr,rowoff,
     #     iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,
     #     ia,ib,isval,ilevpt,ilevnr,inrows,
     #     nrows,nlevs,norbs,
     #     acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #     segsv,levdiv,nnp,c,nwks,h,
     #     g,bcoef,coeffs,maxwks,
     #     nsym,numrow,drtpt,z,z,maxcor,
     #     nnp0,ngsym,numsym,offsym,ijpt)

c
c     ----- un symmetry-block the density matrices -----
c
      call ufmint(orbsym,norbs,nsym,numsym,offsym,ijpt,nnp0,
     #            nnp,ngsym,h,h1,g,g1)
c
c      write(iout,*) '2pdm: b4fixdm'
c      do 123 i=1,nnp
c         write(iout,*) 'kl',i
c         write(iout,*) (g1(j,i),j=1,nnp)
c  123 continue
c     ----- disentangle the one- and two-particle density matrices ----
c
CTEMP
crlm      call fixdm(h1,g1,norbs,nnp)
c
c      write(iout,*) '2pdm: afterfixdm'
c      do 124 i=1,nnp
c         write(iout,*) 'kl',i
c         write(iout,*) (g1(j,i),j=1,nnp)
c  124 continue
c
c     ----- and write out the density matrices to a safe place -----
c
      call iosys('write real "mo 1pdm" to '//dunit,nnp,h1,0,' ')
      call iosys('write real "mo 2pdm" to '//dunit,nnp**2,g1,0,' ')
      write(iout,*) '1pdm:'
      write(iout,*) (h1(j),j=1,nnp)
      write(iout,*) '2pdm:'
      do 122 i=1,nnp
         write(iout,*) 'kl',i
         write(iout,*) (g1(j,i),j=1,nnp)
  122 continue
c
c     --- geminals? ---
c         this feature does not work correctly, because the spin
c         has already been summed over in the 2pdm.
      goto 123
      if(logkey(ops,'geminals',.false.,' ')) then
c        generate the full 2pdm and diagonalize.
         nn=norbs*norbs
         nntri=nn*(nn+1)/2
         eigv=nntri+1
c        should check core available.
         call f2pdm(norbs,nn,nnp,nntri,g1,z(nntri+1),z,iout)
         call rsp(nn,nn,nntri,z,z(eigv),1,z(eigv+nn),z(eigv+nn+nn*nn),
     $            z(eigv+2*nn+nn*nn),
     $            ierr)
         do 100 i=1,nn
            trace=trace+z(eigv+i-1)
  100 continue
         write(iout,*) 'eigenvalue trace',trace
         write(iout,*) 'geminal eigenvalues',(z(eigv+i-1),i=1,nn)
      endif
  123 continue
c
c
      return
      end
      subroutine f2pdm(norbs,nn,nnp,nntri,g1,zsq,z,iout)
      implicit integer(a-z)
c
c     this routine expands the 2pdm into the full array.
c     none of this code has been tested.
      integer norbs,nn,nnp,nntri,iout
      integer ij,kl,i,j,k,l
      real*8 g1(nnp,nnp)
      real*8 zsq(norbs,norbs,norbs,norbs)
      real*8 z(nntri)
      real*8 trace,scale
c     --- the 2pdm is in order ij(electron1),kl(electron2).
c     we want it switched to i(1),k(2),j(1),l(2)
c
      call rzero(zsq,norbs**4)
      trace=0.0d0
      ij=0
      do 100 i=1,norbs
         do 100 j=1,i
            ij=ij+1
            kl=0
            do 90 k=1,norbs
               do 90 l=1,k
                  kl=kl+1
                  scale=1.0d0
                  if((i.eq.j).and.(k.ne.l)) scale=0.5d0
                  if((k.eq.l).and.(i.ne.j)) scale=0.5d0
                  if((i.eq.j).and.(k.eq.l)) scale=0.25d0
                  zsq(i,k,j,l)=zsq(i,k,j,l)+g1(ij,kl)*scale
                  zsq(j,k,i,l)=zsq(j,k,i,l)+g1(ij,kl)*scale
                  zsq(i,l,j,k)=zsq(i,l,j,k)+g1(ij,kl)*scale
                  zsq(j,l,i,k)=zsq(j,l,i,k)+g1(ij,kl)*scale
c                 terms arising from g1(kl,ij)
c                  zsq(k,i,l,j)=g1(ij,kl)
c                  zsq(l,i,k,j)=g1(ij,kl)
c                  zsq(k,j,l,i)=g1(ij,kl)
c                  zsq(l,j,k,i)=g1(ij,kl)
c
c                  ik=(i-1)*norbs+k
c                  jl=(j-1)*norbs+l
c                  if(ik.eq.jl) then
c                     trace=trace+zsq(i,k,j,l)
c                  endif
   90       continue
  100 continue
       write(iout,*) '2pdm trace',trace
       write(iout,*) '2pdm'
       call matout(zsq,nn,nn,nn,nn,iout)
c
c
      call sqtotr(z,zsq,nn,nntri)
         
      return
      end

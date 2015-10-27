*deck @(#)rotcfs.f	4.1  7/7/93
      subroutine rotcfs(c,gamma,t,coeff,done,atprmt,tr,
     #     nfunc,nsymat,nlamda,nrep,lengam,nop,
     #     natoms,angmom,iatom,gampt,relatm,irrep,dump)
c
      implicit integer (a-z)
c
      real*8 c(nfunc,nsymat,nlamda,nrep)
      real*8 t(nfunc,nlamda)
      real*8 gamma(lengam,nop)
      real*8 coeff(nfunc,nlamda,nrep)
      real*8 tr(nfunc,nfunc,nop)
      integer done(nsymat),atprmt(natoms,nop)
      integer relatm(nsymat)
      logical dump
c
      common/io/inp,iout
c
 90   format (/,1x,i5,' vectors of ',i2,'th irrep, degeneracy ',i1,
     #     /,' atoms:',(t10,20i3))
c
c
c     ----- run through the operations, rotating coefficients
c     onto equivalent centres (but not more than once)
c
      call izero(done,nsymat)
      call rzero(c,nfunc*nsymat*nlamda*nrep)
c
      do 500 op=1,nop
         jatom=atprmt(iatom,op)
         do 10 atom=1,nsymat
            if (relatm(atom).eq.jatom) go to 11
 10      continue
         call lnkerr('error with symmetry related atom')
 11      continue
         if (done(atom).eq.1) go to 500
         done(atom)=1
c
         do 400 rep=1,nrep
            if (angmom.eq.0) then
               pt=gampt
               do 2 l=1,nlamda
                  do 1 k=1,nlamda
                     c(1,atom,l,rep)=c(1,atom,l,rep)+
     #                    coeff(1,k,rep)*gamma(pt,op)
                     pt=pt+1
    1             continue
    2          continue
            else
               if (angmom.ge.1) then
                  call ebc(t,tr(1,1,op),coeff(1,1,rep),nfunc,nfunc,
     $                 nlamda)
               end if
               pt=gampt
               do 6 l=1,nlamda
                  do 5 k=1,nlamda
                     do 4 i=1,nfunc
                        c(i,atom,l,rep)=c(i,atom,l,rep)+
     #                       t(i,k)*gamma(pt,op)
    4                continue
                     pt=pt+1
    5             continue
    6          continue
            end if
c
c
 400     continue
 500  continue
c
c     possibly print some debugging information.
      if(dump) then
         write(iout,90) nrep,irrep,nlamda,(relatm(i),i=1,nsymat)
         call matout(c,nfunc*nsymat,nlamda*nrep,nfunc*nsymat,
     #        nlamda*nrep,iout)
      endif
c
c
      return
      end

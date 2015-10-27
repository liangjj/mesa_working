*deck @(#)symprj.f	5.1  11/6/94
      subroutine symprj(nirrep,nbf,nsalc,numso,lambda,lirrep,
     $                  c,nmo,bftoso,csym,lsym,s,t1,t2,t3,t4,
     $                  norm,ops)
c***begin prologue     symprj
c***date written       900410  (yymmdd)
c***revision date      11/30/90
c
c***keywords           symmetry, projection, vectors
c***author             martin, richard(lanl)
c***source             @(#)symprj.f	5.1   11/6/94
c***purpose            projects pure symmetry functions from an arbritrary
c                      input set.
c***description
c
c
c
c***references
c
c***routines called
c     symprj
c       (none)
c
c***end prologue       symprj
      implicit integer(a-z)
      integer numso(nirrep), lambda(nirrep)
      integer lsym(nbf)
      character*(*) lirrep(nirrep), ops
      logical debug, first, logkey, prnt
      real*8 c(nbf,nbf), bftoso(nbf,nsalc), csym(nbf,nbf)
      real*8 s(nbf,nbf), t1(nbf,nbf), t2(nbf,nbf)
      real*8 t3(nbf,nbf), t4(nbf), norm(nbf,nirrep)
      real*8  nrmmax, fac, small, maxerr
      real*8  ovrlap, ovrtst
      parameter (small=1.0d-09, ovrtst=1.0d-06, maxerr=1.0d-06)
      common/io/inp,iout
c
      data first/.true./
      save first
      parameter (debug=.false.)
c
 1000 format(1x,a8,'projected vectors')
 1010 format(1x,a8,'overlap matrix')
 1015 format(1x,'symprj found symmetry contamination')
 1020 format(/6x,'mo:',i4,2x,'irrep:',a8)
c
c     start with a copy of the vectors in csym.
      call scopy(nbf*nbf,c,1,csym,1)
c
      salc=1
      do 100 irrep=1,nirrep
         nsf=numso(irrep)*lambda(irrep)
         if(nsf.ne.0) then
c
c           form the projection operator.
            call ebct(t1,bftoso(1,salc),bftoso(1,salc),nbf,nsf,nbf)
c
c           apply the operator (rho x s x mo's).
            call ebc(t2,t1,s,nbf,nbf,nbf)
            call ebc(t1,t2,csym,nbf,nbf,nbf)
c
            if(debug) then
               write(iout,1000) lirrep(irrep)
               call matout(t1,nbf,nbf,nbf,nbf,iout)
            end if
c
c           find the overlap matrix between projected vectors.
            call ebtc(t2,t1,s,nbf,nbf,nbf)
            call ebc(t3,t2,t1,nbf,nbf,nbf)
            if(debug) then
               write(iout,1010) lirrep(irrep)
               call matout(t3,nbf,nbf,nbf,nbf,iout)
            end if
c
c           save the norms.
            do 60 bf=1,nbf
               norm(bf,irrep)=t3(bf,bf)
   60       continue
c
c           find non-zero projections for this irreducible representation,
c           move them into the output array, and assign a label.
c
c           some of the functions may have been contaminated, and therefore ther
c           will be more functions of representation 1, e.g., than expected. kee
c           the number expected in each group with the largest norms.
            nkept=0
            do 90 sf=1,nsf
c
c              find the largest norm left in the list.
               keep=0
               nrmmax=0.0d+00
               do 80 i=1,nmo
                  if(abs(norm(i,irrep)).gt.nrmmax) then
                     keep=i
                     nrmmax=abs(norm(i,irrep))
                  end if
   80          continue
c
c              make sure this is linearly independent of vectors
c              previously kept. schmidt orthogonalize to the previous vectors.
               do 85 old=1,nkept
                  call ebc(t4,s,t1(1,keep),nbf,nbf,1)
                  call ebtc(ovrlap,t3(1,old),t4,1,nbf,1)
                  call saxpy(nbf,-ovrlap,t3(1,old),1,t1(1,keep),1)
   85          continue
c
c              get the norm of the new orthogonal vector.
               call ebc(t4,s,t1(1,keep),nbf,nbf,1)
               call ebtc(ovrlap,t1(1,keep),t4,1,nbf,1)
               norm(keep,irrep)=ovrlap
c
c              if there is anything left, we wish to keep it.
               if(ovrlap.ge.ovrtst) then
                  nkept=nkept+1
                  fac=1.0d+00/sqrt(abs(ovrlap))
                  call sscal(nbf,fac,t1(1,keep),1,t1(1,keep),1)
                  call scopy(nbf,t1(1,keep),1,t3(1,nkept),1)
                  call scopy(nbf,t1(1,keep),1,csym(1,keep),1)
                  lsym(keep)=irrep
                  norm(keep,irrep)=0.0d+00
               end if
   90       continue
c
c           check to see that things look ok.
            if (nkept.ne.nsf) then
               call lnkerr('symmetry projection cannot '//
     $                     'find all the vectors')
            end if
            salc=salc+nsf
         end if
  100 continue
c
c     report to the user if there are any contaminants left.
      prnt=logkey(ops,'print=m402=contamination-information',
     1               .false.,' ')
      do 200 irrep=1,nirrep
         do 180 mo=1,nbf
            if(abs(norm(mo,irrep)).gt.small) then
c              write message the first time
               if(first) then
                  write(iout,1015)
                  first=.false.
               end if
               if(prnt) then
                  write(iout,1020) mo,lirrep(irrep)
               end if
            end if
  180    continue
  200 continue
c
c     orthonormalize the vectors, once again, just to be sure.
      nnp=nbf*(nbf+1)/2
      call sqtotr(t1,s,nbf,nnp)
      call schmdt(csym,t1,s,t2,t3,nmo,nbf,nnp,maxerr)
c
c
      return
      end

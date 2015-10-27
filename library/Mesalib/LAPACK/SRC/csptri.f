      SUBROUTINE CSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
	implicit none

*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            AP( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CSPTRI compu/k) and not sqrt(1/2k) for photon
C	normalisation.

C	This routine calculates the amplitudes A_3/2 and A_1/2
C	for the process	Delta -> N gamma (for arbitrary k), using 
C	the routine ALLAMP which can find the amplitude for the 
C	operator responsible (at least for the photon attached to 
C	the third quark leg), O_3 between any pair of components of 
C	the wavefunctions. The amplitude for O_2 is found from that
C	of O_3 (in the expanded primed basis) by a transformation
C	of the kind used in the baryon problem, and borrowed directly
C	from there.
C	The wavefunctions for the Nucleon and Delta
C	are read in from files output by the baryon program (and
C	suitably modified), in the form of (eigen)vectors of
C	real coefficients. The corresponding integer 'quantum' 
C	numbers describing the individual components of the wvfn
C	are in exactly the same format as in the baryon program and
C	are read in from the files used there. Likewise for the 
C	transformation coefficients.
	external lnblnk
	integer lnblnk

	character*1 flavf
	integer twojf,paritf,nstatf
	real*8 kmom,amp(8)

	real*8 mq,alpha
	common/setcom/mq,alpha

        character*80 datapath
        character*3 forchar(0:99)
        common/path/datapath, forchar

	real*8 kovera
	common/params/ kovera
	integer	noscli,nosclf,twoji,pariti
	integer nosclfin(7)
	data nosclfin/2,2,1,1,1,1,1/	! number of oscillator levels 
		! for N and Delta states in the I/K model for 
		! 1/2+,3/2+,5/2+,7/2+,1/2-,3/2-,5/2-

	character*1 twicejf(7)
	data twicejf/'1','2','3','4','5','6','7'/
	character*1 parityf(0:1)
	data parityf/'p','m'/
	character*1 nstatef(5)
	data nstatef/'1','2','3','4','5'/
	character*80 filenamef
	character*80 filnoname
	character*80 pfilnoname
	integer filno,pfilno,i,nwi(5),nni,n,qmnoi(120),kwi(5),kki,k,
     &    pqmnoi(240),kzi(120),kfi(120),nwf(5),nnf,qmnof(120),kwf(5),
     &    kkf,pqmnof(240),kzf(120),kff(120),kf,nf,flag,sf,slnr,temp,
     &    lrf,par12f,ki,ni,si,lri,par12i,aind
	real*8 ci(240,120),cf(240,120),eiveci(120),eivecf(120),
     &    o3p(240,240,4),o3(120,120,4),sumi,sumf(4),a(120,120,8)
	common /bigarray/ ci,cf,eiveci,eivecf,
     &    o3p,o3,sumi,sumf,a
	real*8 ap1,ap3,as1,as3,apl,acl,comfac,e1p,m1p,l1p,s1p

	twoji=1
	pariti=0
	noscli=2
	nosclf=nosclfin(4*paritf+(twojf-1)/2+1)
	
	kovera=kmom/alpha

C	Read in the quantum numbers and transformation matrix for the
C	initial baryon.
	
	filno=((twoji-1)/2)*4+11+pariti*2	! see baryon routine 
	pfilno=((twoji-1)/2)*4+12+pariti*2	! DIAGON, these are the
			! quantum number files of the initial baryon

	filnoname = datapath(1:lnblnk(datapath))//'/'
     &       //'for'//forchar(filno)//'.dat'


	pfilnoname = datapath(1:lnblnk(datapath))//'/'
     &       //'for'//forchar(pfilno)//'.dat'

	open(unit=filno,name=filnoname,status='old')
	i=1
	do while(i.le.5)
	  read(filno,*) nwi(i)
	  if(noscli.eq.i) nni=nwi(i)	! nni is the number of 
	  i=i+1				! components in the initial wvfn 
	end do				! to this number of osc. levels
	n=1
	do while(n.le.nni)
	  read(filno,*) qmnoi(n)	! quantum numbers of components
	  n=n+1				! of initial wvfn
	end do
	close(unit=filno,status='keep')

	open(unit=pfilno,name=pfilnoname,status='old')
	i=1
	do while(i.le.5)
	  read(pfilno,*) kwi(i)
	  if(noscli.eq.i) kki=kwi(i)	! kki is the number of 
	  i=i+1				! components of the initial
	end do				! wvfn in the primed basis
	k=1
	do while(k.le.kki)
	  read(pfilno,*) pqmnoi(k)	! quantum numbers of initial
	  k=k+1				! wvfn components in primed
	end do				! basis
	close(unit=pfilno,status='keep')

	filno=((twoji-1)/2)*8! of final wvfn
	  if((lrf/2)*2.ne.lrf) par12f=-par12f
	
	  ki=1
	  ni=1
	  do while(ki.le.kki)
	    si=pqmnoi(ki)/100000	! spin number of initial wvfn
	    slnr=pqmnoi(ki)/1000
	    temp=pqmnoi(ki)-slnr*1000
	    lri=temp/100	! lrho of initial wvfn
	    par12i=1		! (12) exchange parity of this cmpnt
	    if(si.eq.1) par12i=-1		   ! of initial wvfn
	    if((lri/2)*2.ne.lri) par12i=-par12i

	    if(par12i.eq.par12f) then
	      call allamp(twojf,pqmnof(kf),pqmnoi(ki),ap1,ap3,as1,as3,
     &		apl,acl)
		   ! set up an amplitude index, 1=a1,2=a3,3=al,4=as
	      o3p(ki,kf,1)=ap1-as1*kovera	! form these l.c.'s
	      o3p(ki,kf,2)=ap3-as3*kovera	! to save storage,
	      o3p(ki,kf,3)=apl+acl*kovera/2.0d0	! note overall signs
	      o3p(ki,kf,4)=acl			! left until the end
	      if((pqmnoi(ki).eq.qmnoi(ni)).and.(flag.eq.1)) then
	        ni=ni+1			 ! have also found qmnoi(ni) in
		  ! the pqmnoi set so can write	primed basis matrix 
		  ! elts into O_3 matrices; note have already 
		  ! incremented nf and ni
		aind=1
		do while(aind.le.4)
		  o3(ni-1,nf-1,aind)=o3p(ki,kf,aind)
		  aind=aind+1
		end do
 	      endif
	      if(ni.gt.(nni+1)) write(6,500)
500	      format(' error; ni has become larger than nni')
	    endif
	    ki=ki+1
	  end do
	  kf=kf+1
C	Now transform the matrices of the O_3 in the primed basis back 
C	to the usual basis by pre-multiplying by the final baryon 
C	tr'fn matrix and post-multiplying by the initial baryon tr'fn
C	matrix. Note summations for matrix mult are truncated by using
C	the kz and kf's that go with the tr'fn matrices. Note that
C	while we are cycling over nf and ni we may as well also perform
C	the summation of (2x) the O_2 matrix and the O_3 matrix; it is
C	necessary at this point to multiply by the charges e2/e and
C	e3/e. 

	nf=1
	do while(nf.le.nnf)
	  ni=1
	  do while(ni.le.nni)
		! note more efficient to loop over 
		! aind internally, but for now...
	    aind=1	
	    do while(aind.le.4)
	      sumf(aind)=0.0d0
	      kf=kzf(nf)
	      do while(kf.le.kff(nf))
	        ki=kzi(ni)
	        sumi=0.0d0
	        do while(ki.le.kfi(ni))
		  sumi=sumi+o3p(ki,kf,aind)*ci(ki,ni)
		  ki=ki+1
	        end do
	        sumf(aind)=sumf(aind)+sumi*cf(kf,nf)
	        kf=kf+1
	      end do
	      aind=aind+1
	    end do
	    	! set up new amplitude index with two nucleon charge 
		! states, 1=a1p(proton),2=a1n(neutron),3=a3p,4=a3n,
		! 5=alp,6=aln,7=asp,8=asn
	    aind=1
	    do while(aind.le.4)
	      a(ni,nf,2*aind-1)=(4.0d0/3.0d0)*sumf(aind)
     &		-(1.0d0/3.0d0)*o3(ni,nf,aind)
	      a(ni,nf,2*aind)=-(2.0d0/3.0d0)*sumf(aind)
     &		+(2.0d0/3.0d0)*o3(ni,nf,aind)
	      aind=aind+1
	    end do
	    ni=ni+1
	  end do
	  nf=nf+1
	end do

C	Now reduce the resulting matrices to numbers by premultiplying
C	by the eigenvector representing the final baryon (i.e. by its
C	transpose) and post-multiplying by the eigenvector representing
C	the initial baryon.

	aind=1
	do while(aind.le.8)
	  amp(aind)=0.0d0
	  nf=1
	  do while(nf.le.nnf)
	    sumi=0.0d0
	    ni=1
	    do while(ni.le.nni)
	      sumi=sumi+a(ni,nf,aind)*eiveci(ni)
	      ni=ni+1
	    end do
	    amp(aind)=amp(aind)+eivecf(nf)*sumi
	    nf=nf+1
	  end do
	  aind=aind+1
	end do

C	Now multiply the result by the overall factors suppressed
C	above and return the results along with the parameters used.

	comfac=dexp(-kovera**2/6.0d0)
	aind=1
	do while(aind.le.4)
	  amp(aind)=amp(aind)*comfac*alpha/(dsqrt(2.0d0)*mq)
	  aind=aind+1
	end do
	amp(5)=-amp(5)*comfac*alpha/mq
	amp(6)=-amp(6)*comfac*alpha/mq
	amp(7)=amp(7)*comfac
	amp(8)=amp(8)*comfac

C	Need to change the sign of the neutron amplitudes in the case 
C	that the final state is a Delta(^0), because of the sign in 
C	the projection <phi^lambda_n|ddu>.

	if(flavf.eq.'d') then
	  aind=1
	  do while(aind.le.4)
	    amp(2*aind)=-amp(2*aind)
	    aind=aind+1
	  end do
	end if

	return
	end

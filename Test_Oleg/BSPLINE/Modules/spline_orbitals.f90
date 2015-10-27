!=======================================================================
      MODULE spline_orbitals
!=======================================================================
!
!     contains description of atomic orbitals 
!     and their B-spline representation
!
!-----------------------------------------------------------------------
   
      IMPLICIT NONE
      SAVE

! ... LIST OF ONE-ELECTRON ORBITALS

      INTEGER(4) :: mbf = 0        ! max. number of orbitals
      INTEGER(4) :: nbf = 0        ! current number of orbitals
      INTEGER(4) :: ibf = 512      ! initial prediction of mbf
      INTEGER(4) :: jbf = 512      ! incriment for mbf  
      INTEGER(4) :: iscrbf = 99    ! scratch file for reallocation
    
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: nbs   !  n-values
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: lbs   !  l-values
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: kbs   !  set numbers
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: mbs   !  number of splines
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: iech  !  additional pointer
    
      CHARACTER(4), ALLOCATABLE, DIMENSION(:) :: ebs ! spectroscopic notation

! ... B-spline expansion coefficients:

      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: PBS 

! ... convolution with B-overlaps: 

      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: QBS 

! ... orbital orthogonality and AFTER conditions

      INTEGER(4) :: JBORT = 1        
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: IBORT

! ... one-electron overlaps:

      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: OBS


      END MODULE spline_orbitals



!=======================================================================
      SUBROUTINE allocate_bsorb(m)
!=======================================================================
!
!     This program allocates (deallocates) space for list of atomic
!     orbitals or reallocates it if necessary
!
!-----------------------------------------------------------------------
    
      USE spline_orbitals
      USE spline_param

      IMPLICIT NONE

      INTEGER(4), INTENT(in) :: m
      INTEGER(4) :: i,j,met

      if(m.le.0) then
      
       if(Allocated(NBS)) &
	    Deallocate (NBS,LBS,KBS,MBS,iech,EBS,IBORT,PBS,QBS,OBS)
       nbf = 0;  mbf = 0
      
      elseif(m.gt.mbf.or.(nbf.gt.0.and.m.ge.nbf)) then

       met = 0
       if(nbf.gt.0.and.Allocated(NBS)) then
        open(iscrbf,form='UNFORMATTED',status='SCRATCH')
        rewind(iscrbf)
        Do i=1,nbf
         write(iscrbf) nbs(i),lbs(i),kbs(i),ebs(i),mbs(i),iech(i), &
               IBORT(1:nbf,i),PBS(1:ns,i),QBS(1:ns,i),OBS(1:nbf,i)
        End do
        met = 1
       end if
       if(Allocated(nbs)) & 
          Deallocate (nbs,lbs,kbs,mbs,iech,ebs,IBORT,PBS,QBS,OBS)
       
       mbf = m
       Allocate(nbs(mbf),lbs(mbf),kbs(mbf),ebs(mbf),mbs(1:mbf), &
                iech(1:mbf),IBORT(1:mbf,1:mbf),PBS(1:ns,1:mbf), &
                QBS(1:ns,1:mbf),OBS(1:mbf,1:mbf))
       nbs = 0; lbs = 0; kbs = 0; ebs = '****'; mbs = 0; iech = 0
       IBORT = 2; PBS = 0.d0; QBS = 0.d0; OBS = 0.d0

       if(met.gt.0) then
        rewind(iscrbf)
        Do i=1,nbf
         read(iscrbf) nbs(i),lbs(i),kbs(i),ebs(i),mbs(i),iech(i), &
               IBORT(1:nbf,i),PBS(1:ns,i),QBS(1:ns,i),OBS(1:nbf,i)
        End do
        close(iscrbf)
!        write(*,*) 'realoc_BS_orb: mbf=',mbf
       end if

      end if

      END SUBROUTINE allocate_bsorb



!=======================================================================
      INTEGER(4) FUNCTION Ifind_bsorb(n,l,k)
!=======================================================================
      
      USE spline_orbitals

      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: n,l,k
      INTEGER(4) :: i
      
        Ifind_bsorb=0

        Do i=1,nbf
         if(n.eq.nbs(i).and.l.eq.lbs(i).and.k.eq.kbs(i)) then
          Ifind_bsorb = i
          Return
         end if
        End do

      END FUNCTION Ifind_bsorb


!=======================================================================
      INTEGER(4) FUNCTION Jfind_bsorb(n,l,k)
!=======================================================================
     
      USE spline_orbitals

      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: n,l,k
      INTEGER(4), External :: Ifind_bsorb
      INTEGER(4) :: i

      i =  Ifind_bsorb(n,l,k)
      if(i.eq.0) then
       Write(*,*) ' Jfind_bsorb: can not find the orbital: N,L,K=',n,l,k
       Stop
      else
       Jfind_bsorb = i
      end if

      END FUNCTION Jfind_bsorb


!=======================================================================
      Integer(4) Function Iadd_bsorb(n,l,k)
!=======================================================================
!
!     adds the orbital to the list
!
!-----------------------------------------------------------------------

      USE spline_orbitals

      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: n,l,k
      INTEGER(4), External :: Ifind_bsorb
      CHARACTER(4), EXTERNAL :: ELF4
      INTEGER(4) :: i

      i =  Ifind_bsorb(n,l,k)

      if(i.eq.0) then

       if(nbf+1.gt.mbf) Call Allocate_bsorb(mbf+jbf)
       nbf = nbf + 1
       nbs(nbf) = n
       lbs(nbf) = l
       kbs(nbf) = k
       ebs(nbf) = ELF4(n,l,k)
       Iadd_bsorb = nbf

      else

       Iadd_bsorb = i

      end if

      END FUNCTION Iadd_bsorb



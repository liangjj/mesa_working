*deck aaaaaa
      subroutine aaaaaa (ver)
c***begin prologue  aaaaaa
c***purpose  slatec common mathematical library disclaimer and version.
c***library   slatec
c***category  z
c***type      all (aaaaaa-a)
c***keywords  disclaimer, documentation, version
c***author  slatec common mathematical library committee
c***description
c
c   the slatec common mathematical library is issued by the following
c
c           air force weapons laboratory, albuquerque
c           lawrence livermore national laboratory, livermore
c           los alamos national laboratory, los alamos
c           national institute of standards and technology, washington
c           national energy research supercomputer center, livermore
c           oak ridge national laboratory, oak ridge
c           sandia national laboratories, albuquerque
c           sandia national laboratories, livermore
c
c   all questions concerning the distribution of the library should be
c   directed to the national energy software center, 9700 cass ave.,
c   argonne, illinois  60439, and not to the authors of the subprograms.
c
c                    * * * * * notice * * * * *
c
c   this material was prepared as an account of work sponsored by the
c   united states government.  neither the united states, nor the
c   department of energy, nor the department of defense, nor any of
c   their employees, nor any of their contractors, subcontractors, or
c   their employees, makes any warranty, expressed or implied, or
c   assumes any legal liability or responsibility for the accuracy,
c   completeness, or usefulness of any information, apparatus, product,
c   or process disclosed, or represents that its use would not infringe
c   upon privately owned rights.
c
c *usage:
c
c        character * 16 ver
c
c        call aaaaaa (ver)
c
c *arguments:
c
c     ver:out   will contain the version number of the slatec cml.
c
c *description:
c
c   this routine contains the slatec common mathematical library
c   disclaimer and can be used to return the library version number.
c
c***references  kirby w. fong, thomas h. jefferson, tokihiko suyehiro
c                 and lee walton, guide to the slatec common mathema-
c                 tical library, april 10, 1990.
c***routines called  (none)
c***revision history  (yymmdd)
c   800424  date written
c   890414  revision date from version 3.2
c   890713  routine modified to return version number.  (wrb)
c   900330  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c   921215  updated for version 4.0.  (wrb)
c   930701  updated for version 4.1.  (wrb)
c***end prologue  aaaaaa
      character * (*) ver
c***first executable statement  aaaaaa
      ver = ' 4.1'
      return
      end

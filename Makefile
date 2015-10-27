#  @(#)Makefile	5.1 11/6/94
#
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
# If there are other subdirectories that you use they also should be set.
#
#            $MESA_HOME, $MESA_LIB, $MESA_BIN and MESA_TMP
#
#
#
#		$MESA_HOME
#		$MESA_BIN
#		$MESA_LIB
#		$MESA_TMP
#
#			FEDVR_Lib		f90 code
#			Mesalib = UTILITY_LIB now
#				Character_Manipulation_Subroutines
#				Common_los_Alamos_Mathematical_Subroutines
#				Double_Precision_Integral_Library
#				General_Utility_Subroutines
#				Integral_Library
#				IOsys_Subroutines
#				Machine_Dependent_Subroutines
#				Mathematical_Subroutines
#				Symmetry_Library
#				LAPACK
#				LAPACK95
#				slatec
#			Modules			f90 code
#				General_Modules
#				Mathematics_Modules
#				Utility_Modules
#				General_Modules
#			Potential		f90 code
#			Utilities		
#				Same structure as Mesalib
#
#		m0
#		m1
#		...
#		m2001
#		source
#		mesa.dat	
#		mesa.sh
#		resfit
#		Test_Link
#		testmake
#		Test_Oleg
#		tests
#		Time_Propagation
#		bmesa
#			Time_Propagation		
#		DC_DIAG
#		FEDVR_Driver
#		Legendre
#
MFLAGS = $(flags)
LIBRARIES =	$(LIBRARY)


M0 =			$(MESA_HOME)/m0
M1 =			$(MESA_HOME)/m1
M101 =			$(MESA_HOME)/m101
M102 =			$(MESA_HOME)/m102
M104 =			$(MESA_HOME)/m104
M2001 =			$(MESA_HOME)/m2001
M104 =			$(MESA_HOME)/m104
M201 =			$(MESA_HOME)/m201
M202 =			$(MESA_HOME)/m202
M202_BIS =		$(MESA_HOME)/m202bis
M301 =			$(MESA_HOME)/m301
M302 =			$(MESA_HOME)/m302
M312 =			$(MESA_HOME)/m312
M401 =			$(MESA_HOME)/m401
M402 =			$(MESA_HOME)/m402
M403 =			$(MESA_HOME)/m403
M411 =			$(MESA_HOME)/m411
M412 =			$(MESA_HOME)/m412
M501 =			$(MESA_HOME)/m501
M502 =			$(MESA_HOME)/m502
M503 =			$(MESA_HOME)/m503
M601 =			$(MESA_HOME)/m601
M602 =			$(MESA_HOME)/m602
M604 =			$(MESA_HOME)/m604
M330 =			$(MESA_HOME)/m330
M801 =			$(MESA_HOME)/m801
M802 =			$(MESA_HOME)/m802
M805 =			$(MESA_HOME)/m805
M805_NEW =		$(MESA_HOME)/m805_new
M806 =			$(MESA_HOME)/m806
M806_BIS =		$(MESA_HOME)/m806bis
M811 =			$(MESA_HOME)/m811
M811_BIS =		$(MESA_HOME)/m811bis
M811_OLD =		$(MESA_HOME)/m811old
M819 =			$(MESA_HOME)/m819
M820 =			$(MESA_HOME)/m820
M821 =			$(MESA_HOME)/m821
M822 =			$(MESA_HOME)/m822
M901 =			$(MESA_HOME)/m901
M902 =			$(MESA_HOME)/m902
M903 =			$(MESA_HOME)/m903
M930 =			$(MESA_HOME)/m930
M935 =			$(MESA_HOME)/m935
M940 =			$(MESA_HOME)/m940
M940_BIS =		$(MESA_HOME)/m940bis
M103 =			$(MESA_HOME)/m103
M1001 =			$(MESA_HOME)/m1001
M1003 =			$(MESA_HOME)/m1003
M1011 =			$(MESA_HOME)/m1011
M1012 =			$(MESA_HOME)/m1012
M1013 =			$(MESA_HOME)/m1013
M1014 =			$(MESA_HOME)/m1014
M1020 =			$(MESA_HOME)/m1020
M1021 =			$(MESA_HOME)/m1021
M1022 =			$(MESA_HOME)/m1022
M1031 =			$(MESA_HOME)/m1031
M1032 =			$(MESA_HOME)/m1032
M1033 =			$(MESA_HOME)/m1033
M1402 =			$(MESA_HOME)/m1402
M1901 =			$(MESA_HOME)/m1902
M1951 =			$(MESA_HOME)/m1951
M1990 =			$(MESA_HOME)/m1990
M1991 =			$(MESA_HOME)/m1991
M203 =			$(MESA_HOME)/m203
M204 =			$(MESA_HOME)/m204
M205 =			$(MESA_HOME)/m205
M206 =			$(MESA_HOME)/m206
M250 =			$(MESA_HOME)/m250
M303 =			$(MESA_HOME)/m303
M304 =			$(MESA_HOME)/m304
M305 =			$(MESA_HOME)/m305
M309 =			$(MESA_HOME)/m309
M319 =			$(MESA_HOME)/m319
M323 =			$(MESA_HOME)/m323
M390 =			$(MESA_HOME)/m390
M395 =			$(MESA_HOME)/m395
M4001 =			$(MESA_HOME)/m4001
M5000 =			$(MESA_HOME)/m5000
M5001 =			$(MESA_HOME)/m5001
M511 =			$(MESA_HOME)/m511
M515 =			$(MESA_HOME)/m515
M551 =			$(MESA_HOME)/m551
M553 =			$(MESA_HOME)/m553
M6000 =			$(MESA_HOME)/m6000
M6001 =			$(MESA_HOME)/m6001
M6002 =			$(MESA_HOME)/m6002
M6003 =			$(MESA_HOME)/m6003
M6004 =			$(MESA_HOME)/m6004
M6005 =			$(MESA_HOME)/m6005
M6008 =			$(MESA_HOME)/m6008
M6009 =			$(MESA_HOME)/m6009
M6010 =			$(MESA_HOME)/m6010
M6011 =			$(MESA_HOME)/m6011
M6016 =			$(MESA_HOME)/m6016
M6018 =			$(MESA_HOME)/m6018
M6020 =			$(MESA_HOME)/m6020
M6025 =			$(MESA_HOME)/m6025
M6026 =			$(MESA_HOME)/m6026
M6030 =			$(MESA_HOME)/m6030
M6031 =			$(MESA_HOME)/m6031
M6032 =			$(MESA_HOME)/m6032
M6050 =			$(MESA_HOME)/m6050
M6055 =			$(MESA_HOME)/m6055
M6060 =			$(MESA_HOME)/m6060
M6065 =			$(MESA_HOME)/m6065
M6070 =			$(MESA_HOME)/m6070
M6090 =			$(MESA_HOME)/m6090
M6100 =			$(MESA_HOME)/m6100
M611 =			$(MESA_HOME)/m611
M613 =			$(MESA_HOME)/m613
M518 =			$(MESA_HOME)/m618
M619 =			$(MESA_HOME)/m619
M620 =			$(MESA_HOME)/m620
M6200 =			$(MESA_HOME)/m6200
M6201 =			$(MESA_HOME)/m6201
M6202 =			$(MESA_HOME)/m6202
M6202_T =		$(MESA_HOME)/m6202t
M6203 =			$(MESA_HOME)/m6203
M621 =			$(MESA_HOME)/m621
M7003 =			$(MESA_HOME)/m7003
M701 =			$(MESA_HOME)/m701
M702 =			$(MESA_HOME)/m702
M711 =			$(MESA_HOME)/m711
M712 =			$(MESA_HOME)/m712
M715 =			$(MESA_HOME)/m715
M721 =			$(MESA_HOME)/m721
M725 =			$(MESA_HOME)/m725
M731 =			$(MESA_HOME)/m731
M732 =			$(MESA_HOME)/m732
M7777 =			$(MESA_HOME)/m7777
M8000 =			$(MESA_HOME)/m8000
M807 =			$(MESA_HOME)/m807
M812 =			$(MESA_HOME)/m812
M824 =			$(MESA_HOME)/m824
M830 =			$(MESA_HOME)/m830
M840 =			$(MESA_HOME)/m840
M891 =			$(MESA_HOME)/m891
M904 =			$(MESA_HOME)/m904
M905 =			$(MESA_HOME)/m905
M910 =			$(MESA_HOME)/m910
M911 =			$(MESA_HOME)/m911
M914 =			$(MESA_HOME)/m914
M915 =			$(MESA_HOME)/m915
M916 =			$(MESA_HOME)/m916
M921 =			$(MESA_HOME)/m921
M929 =			$(MESA_HOME)/m929
M941 =			$(MESA_HOME)/m941
M942 =			$(MESA_HOME)/m942
M990 =			$(MESA_HOME)/m990
M991 =			$(MESA_HOME)/m991
M992 =			$(MESA_HOME)/m992
M993 =			$(MESA_HOME)/m993
M994 =			$(MESA_HOME)/m994
M999 =			$(MESA_HOME)/m999
RESFIT =		$(MESA_HOME)/resfit
TEST_LINK =		$(MESA_HOME)/Test_Link
TESTMAKE =		$(MESA_HOME)/testmake
TEST_OLEG =		$(MESA_HOME)/Test_Oleg
TESTS =			$(MESA_HOME)/tests
TIME_PROPAGATION_1 =	$(MESA_HOME)/Time_Propagation
BMESA =			$(MESA_HOME)/bmesa
TIME_PROPAGATION_2 =	$(MESA_HOME)/bmesa/Time_Propagation
DC_DIAG =		$(MESA_HOME)/DC_DIAG
FEDVR_DRIVER =		$(MESA_HOME)/FEDVR_Driver
LEGENDRE =		$(MESA_HOME)/Legendre

ALL_LINKS = 		$(M0) \
			$(M1) \
			$(M101) \
			$(M102) \
			$(M202) \
			$(M202_BIS) \
			$(M301) \
			$(M302) \
			$(M312) \
			$(M412) \
			$(M501) \
			$(M502) \
			$(M503) \
			$(M602) \
			$(M806) \
			$(M806_BIS) \
			$(M811_OLD) \
			$(M811) \
			$(M811_BIS) \
			$(M819) \
			$(M820) \
			$(M821) \
			$(M822) \
			$(M930) \
			$(M940) \
			$(M940_BIS) \
			$(M2001) \
			$(M103) \
			$(M1001) \
			$(M1003) \
			$(M1011) \
			$(M1012) \
			$(M1013) \
			$(M1014) \
			$(M1020) \
			$(M1021) \
			$(M1022) \
			$(M1031) \
			$(M1032) \
			$(M1033) \
			$(M1402) \
			$(M1902) \
			$(M1951) \
			$(M1990) \
			$(M1991) \
			$(M203) \
			$(M204) \
			$(M205) \
			$(M206) \
			$(M250) \
			$(M303) \
			$(M304) \
			$(M305) \
			$(M309) \
			$(M319) \
			$(M323) \
			$(M390) \
			$(M395) \
			$(M4001) \
			$(M511) \
			$(M515) \
			$(M551) \
			$(M553) \
			$(M6000) \
			$(M6001) \
			$(M6002) \
			$(M6003) \
			$(M6004) \
			$(M6005) \
			$(M6008) \
			$(M6009) \
			$(M6010) \
			$(M6011) \
			$(M6016) \
			$(M6018) \
			$(M6020) \
			$(M6025) \
			$(M6026) \
			$(M6030) \
			$(M6031) \
			$(M6032) \
			$(M6050) \
			$(M6055) \
			$(M6060) \
			$(M6065) \
			$(M6070) \
			$(M6090) \
			$(M6100) \
			$(M611) \
			$(M613) \
			$(M618) \
			$(M619) \
			$(M620) \
			$(M6200) \
			$(M6201) \
			$(M6202) \
			$(M6202t) \
			$(M6203) \
			$(M621) \
			$(M7003) \
			$(M701) \
			$(M702) \
			$(M711) \
			$(M712) \
			$(M715) \
			$(M721) \
			$(M725) \
			$(M731) \
			$(M732) \
			$(M7777) \
			$(M8000) \
			$(M807) \
			$(M812) \
			$(M824) \
			$(M830) \
			$(M840) \
			$(M891) \
			$(M904) \
			$(M905) \
			$(M910) \
			$(M911) \
			$(M914) \
			$(M915) \
			$(M916) \
			$(M921) \
			$(M929) \
			$(M941) \
			$(M942) \
			$(M990) \
			$(M991) \
			$(M992) \
			$(M993) \
			$(M994) \
			$(M999)

STARTUP =		$(M0) \
			$(M1) \
			$(M103) \
			$(M104) 

GEOM =			$(M101) \
			$(M201) \
			$(M202) \
			$(M202_bis) \
			$(M203) \
			$(M204) \
			$(M205) \
			$(M250)

BASIS =			$(M102) \
			$(M206)

INTEGRALS =		$(M301) \
			$(M302) \
			$(M303) \
			$(M312) \
			$(M330) 

PROPERTIES =		$(M1902) \
			$(M1951) \
			$(M1991)

ORBITALS =		$(M401) \
			$(M402) \
			$(M412)

SCF =	 		$(M501) \
			$(M502) \
			$(M503)

MCSCF =			$(M1001) \
			$(M1003) \
			$(M1013)

LAGRANGIAN =		$(M1011) \
			$(M1031) \
			$(M1032) \
			$(M1033)

HESSIAN =		$(M1012) \
			$(M1014)

CPHF =			$(M1013) \
			$(M1020) \
			$(M1021) \
			$(M1022)

GRIDMO =		$(M1990)

modules: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) 

libraries: FORCE
	cd $(LIBRARY) ; $(MAKE) $(MFLAGS)

precursors: FORCE
	cd $(LIBRARY) ; $(MAKE) $(MFLAGS)
	cd $(FEDVR_HOME) ; $(MAKE) $(MFLAGS)
	cd $(FEDVR_LIB) ; $(MAKE) $(MFLAGS)

links_0: FORCE
	cd $(M0) ; $(MAKE) $(MFLAGS) 
	cd $(M1) ; $(MAKE) $(MFLAGS) 
	cd $(M101) ; $(MAKE) $(MFLAGS) 
	cd $(M102) ; $(MAKE) $(MFLAGS) 
	cd $(M202_BIS) ; $(MAKE) $(MFLAGS) 
	cd $(M301) ; $(MAKE) $(MFLAGS) 
	cd $(M302) ; $(MAKE) $(MFLAGS) 
	cd $(M312) ; $(MAKE) $(MFLAGS) 
	cd $(M330) ; $(MAKE) $(MFLAGS) 
	cd $(M412) ; $(MAKE) $(MFLAGS) 
	cd $(M501) ; $(MAKE) $(MFLAGS) 
	cd $(M502) ; $(MAKE) $(MFLAGS) 
	cd $(M503) ; $(MAKE) $(MFLAGS) 
	cd $(M602) ; $(MAKE) $(MFLAGS) 
	cd $(M806_BIS) ; $(MAKE) $(MFLAGS) 
	cd $(M811_BIS) ; $(MAKE) $(MFLAGS) 
	cd $(M819) ; $(MAKE) $(MFLAGS) 
	cd $(M822) ; $(MAKE) $(MFLAGS) 
	cd $(M921) ; $(MAKE) $(MFLAGS) 
	cd $(M929) ; $(MAKE) $(MFLAGS) 
	cd $(M930) ; $(MAKE) $(MFLAGS) 
	cd $(M935) ; $(MAKE) $(MFLAGS) 
	cd $(M940_BIS) ; $(MAKE) $(MFLAGS) 
	cd $(M2001) ; $(MAKE) $(MFLAGS) 

links_1: FORCE
	cd $(M103) ; $(MAKE) $(MFLAGS)  
	cd $(M1001) ; $(MAKE) $(MFLAGS)  
	cd $(M1003) ; $(MAKE) $(MFLAGS)  
	cd $(M1011) ; $(MAKE) $(MFLAGS)  
	cd $(M1012) ; $(MAKE) $(MFLAGS)  
	cd $(M1013) ; $(MAKE) $(MFLAGS)  
	cd $(M1014) ; $(MAKE) $(MFLAGS) 
	cd $(M1020) ; $(MAKE) $(MFLAGS) 
	cd $(M1021) ; $(MAKE) $(MFLAGS) 
	cd $(M1022) ; $(MAKE) $(MFLAGS) 
	cd $(M1031) ; $(MAKE) $(MFLAGS) 
	cd $(M1032) ; $(MAKE) $(MFLAGS) 
	cd $(M1033) ; $(MAKE) $(MFLAGS) 
	cd $(M1402) ; $(MAKE) $(MFLAGS) 
	cd $(M1902) ; $(MAKE) $(MFLAGS) 
	cd $(M1951) ; $(MAKE) $(MFLAGS) 
	cd $(M1990) ; $(MAKE) $(MFLAGS) 
	cd $(M1991) ; $(MAKE) $(MFLAGS) 
	cd $(M203) ; $(MAKE) $(MFLAGS) 
	cd $(M204) ; $(MAKE) $(MFLAGS) 
	cd $(M205) ; $(MAKE) $(MFLAGS) 
	cd $(M206) ; $(MAKE) $(MFLAGS) 
	cd $(M250) ; $(MAKE) $(MFLAGS) 
	cd $(M303) ; $(MAKE) $(MFLAGS) 
	cd $(M304) ; $(MAKE) $(MFLAGS) 
	cd $(M305) ; $(MAKE) $(MFLAGS) 
	cd $(M309) ; $(MAKE) $(MFLAGS) 
	cd $(M319) ; $(MAKE) $(MFLAGS) 
	cd $(M323) ; $(MAKE) $(MFLAGS) 
	cd $(M390) ; $(MAKE) $(MFLAGS) 
	cd $(M395) ; $(MAKE) $(MFLAGS) 
	cd $(M4001) ; $(MAKE) $(MFLAGS) 
	cd $(M511) ; $(MAKE) $(MFLAGS) 
	cd $(M515) ; $(MAKE) $(MFLAGS) 
	cd $(M551) ; $(MAKE) $(MFLAGS) 
	cd $(M553) ; $(MAKE) $(MFLAGS) 
	cd $(M6000) ; $(MAKE) $(MFLAGS) 
	cd $(M6001) ; $(MAKE) $(MFLAGS) 
	cd $(M6002) ; $(MAKE) $(MFLAGS) 
	cd $(M6003) ; $(MAKE) $(MFLAGS) 
	cd $(M6004) ; $(MAKE) $(MFLAGS) 
	cd $(M6005) ; $(MAKE) $(MFLAGS) 
	cd $(M6008) ; $(MAKE) $(MFLAGS) 
	cd $(M6009) ; $(MAKE) $(MFLAGS) 
	cd $(M6010) ; $(MAKE) $(MFLAGS) 
	cd $(M6011) ; $(MAKE) $(MFLAGS) 
	cd $(M6016) ; $(MAKE) $(MFLAGS) 
	cd $(M6018) ; $(MAKE) $(MFLAGS) 
	cd $(M6020) ; $(MAKE) $(MFLAGS) 
	cd $(M6025) ; $(MAKE) $(MFLAGS) 
	cd $(M6026) ; $(MAKE) $(MFLAGS) 
	cd $(M6030) ; $(MAKE) $(MFLAGS) 
	cd $(M6031) ; $(MAKE) $(MFLAGS) 
	cd $(M6032) ; $(MAKE) $(MFLAGS) 
	cd $(M6050) ; $(MAKE) $(MFLAGS) 
	cd $(M6055) ; $(MAKE) $(MFLAGS) 
	cd $(M6060) ; $(MAKE) $(MFLAGS) 
	cd $(M6065) ; $(MAKE) $(MFLAGS) 
	cd $(M6070) ; $(MAKE) $(MFLAGS) 
	cd $(M6090) ; $(MAKE) $(MFLAGS) 
	cd $(M6100) ; $(MAKE) $(MFLAGS) 
	cd $(M611) ; $(MAKE) $(MFLAGS) 
	cd $(M613) ; $(MAKE) $(MFLAGS) 
	cd $(M618) ; $(MAKE) $(MFLAGS) 
	cd $(M619) ; $(MAKE) $(MFLAGS) 
	cd $(M620) ; $(MAKE) $(MFLAGS) 
	cd $(M6200) ; $(MAKE) $(MFLAGS) 
	cd $(M6201) ; $(MAKE) $(MFLAGS) 
	cd $(M6202) ; $(MAKE) $(MFLAGS) 
	cd $(M6202t) ; $(MAKE) $(MFLAGS) 
	cd $(M621) ; $(MAKE) $(MFLAGS) 
	cd $(M7003) ; $(MAKE) $(MFLAGS) 
	cd $(M701) ; $(MAKE) $(MFLAGS) 
	cd $(M702) ; $(MAKE) $(MFLAGS) 
	cd $(M711) ; $(MAKE) $(MFLAGS) 
	cd $(M712) ; $(MAKE) $(MFLAGS) 
	cd $(M715) ; $(MAKE) $(MFLAGS) 
	cd $(M721) ; $(MAKE) $(MFLAGS) 
	cd $(M725) ; $(MAKE) $(MFLAGS) 
	cd $(M731) ; $(MAKE) $(MFLAGS) 
	cd $(M732) ; $(MAKE) $(MFLAGS) 
	cd $(M7777) ; $(MAKE) $(MFLAGS) 
	cd $(M8000) ; $(MAKE) $(MFLAGS) 
	cd $(M807) ; $(MAKE) $(MFLAGS) 
	cd $(M812) ; $(MAKE) $(MFLAGS) 
	cd $(M824) ; $(MAKE) $(MFLAGS) 
	cd $(M830) ; $(MAKE) $(MFLAGS) 
	cd $(M840) ; $(MAKE) $(MFLAGS) 
	cd $(M891) ; $(MAKE) $(MFLAGS) 
	cd $(M904) ; $(MAKE) $(MFLAGS) 
	cd $(M905) ; $(MAKE) $(MFLAGS) 
	cd $(M910) ; $(MAKE) $(MFLAGS) 
	cd $(M911) ; $(MAKE) $(MFLAGS) 
	cd $(M914) ; $(MAKE) $(MFLAGS) 
	cd $(M915) ; $(MAKE) $(MFLAGS) 
	cd $(M916) ; $(MAKE) $(MFLAGS) 
	cd $(M921) ; $(MAKE) $(MFLAGS) 
	cd $(M929) ; $(MAKE) $(MFLAGS) 
	cd $(M941) ; $(MAKE) $(MFLAGS) 
	cd $(M942) ; $(MAKE) $(MFLAGS) 
	cd $(M990) ; $(MAKE) $(MFLAGS) 
	cd $(M991) ; $(MAKE) $(MFLAGS) 
	cd $(M992) ; $(MAKE) $(MFLAGS) 
	cd $(M993) ; $(MAKE) $(MFLAGS) 
	cd $(M994) ; $(MAKE) $(MFLAGS) 
	cd $(M999) ; $(MAKE) $(MFLAGS) 

svn: FORCE
	cd $(M0) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M101) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M102) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M202) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M202_BIS) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M301) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M302) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M312) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M412) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M501) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M502) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M503) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M602) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M806) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M806_BIS) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M811_OLD) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M811) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M811_BIS) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M819) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M820) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M821) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M822) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M930) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M940) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M940_BIS) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M2001) ; $(MAKE) $(MFLAGS)	clean_svn
	cd $(M103) ; $(MAKE) $(MFLAGS) clean_svn 
	cd $(M1001) ; $(MAKE) $(MFLAGS) clean_svn 
	cd $(M1003) ; $(MAKE) $(MFLAGS) clean_svn 
	cd $(M1011) ; $(MAKE) $(MFLAGS) clean_svn 
	cd $(M1012) ; $(MAKE) $(MFLAGS) clean_svn 
	cd $(M1013) ; $(MAKE) $(MFLAGS) clean_svn 
	cd $(M1014) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1020) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1021) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1022) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1031) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1032) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1033) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1402) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1902) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1951) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1990) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M1991) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M203) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M204) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M205) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M206) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M250) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M303) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M304) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M305) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M309) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M319) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M323) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M390) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M395) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M4001) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M511) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M515) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M551) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M553) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6000) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6001) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6002) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6003) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6004) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6005) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6008) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6009) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6010) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6011) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6016) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6018) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6020) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6025) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6026) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6030) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6031) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6032) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6050) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6055) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6060) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6065) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6070) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6090) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6100) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M611) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M613) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M618) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M619) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M620) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6200) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6201) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6202) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M6202t) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M621) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M7003) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M701) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M702) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M711) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M712) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M715) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M721) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M725) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M731) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M732) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M7777) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M8000) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M807) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M812) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M824) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M830) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M840) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M891) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M904) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M905) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M910) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M911) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M914) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M915) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M916) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M921) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M929) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M941) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M942) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M990) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M991) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M992) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M993) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M994) ; $(MAKE) $(MFLAGS) clean_svn
	cd $(M999) ; $(MAKE) $(MFLAGS) clean_svn

FORCE:


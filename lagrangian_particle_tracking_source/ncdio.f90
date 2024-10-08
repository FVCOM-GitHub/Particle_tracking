!!$/=============================================================================/
!!$Copyright (c) 2007, The University of Massachusetts Dartmouth 
!!$Produced at the School of Marine Science & Technology 
!!$Marine Ecosystem Dynamics Modeling group
!!$All rights reserved.
!!$
!!$The FVCOM Offline Lagrangian Model has been developed by the joint UMASSD-WHOI
!!$research team.   For details of authorship and attribution of credit please see
!!$the FVCOM technical manual or contact the MEDM group.
!!$
!!$ 
!!$This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu/ The
!!$full copyright notice is contained in the file COPYRIGHT located in the root
!!$directory of the FVCOM code. This original header must be maintained in all
!!$distributed versions.
!!$
!!$THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!!$ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
!!$IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
!!$ARE  DISCLAIMED.  
!!$
!!$/-----------------------------------------------------------------------------/
!!$CVS VERSION INFORMATION
!!$$Id: $
!!$$Name: $
!!$$Revision: $
!!$/=============================================================================/
! VERSION 2.0
! created by Martin Huret
! modified by J. Churchill



!    *** VERSION CREATED BY JHC (04/07) TO GET TIME FROM NETCDF FILE ON NCD_READ CALL ***
!    ***   AND TO OUTPUT SURFACE ELEVATION AND BOTTOM DEPTH ON NCD_WRITE CALL         ***
!      *** modified (08/07) to re-try opening the input cdf file on error  ******

SUBROUTINE NCD_READ_GRID(INFILE)
  !---------------------------------------------------------------------
  ! READ DIMENSIONS IN A NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_NCD
  USE LIMS
  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  CHARACTER(LEN=100), INTENT(IN) :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: N_ELEMS,N_NODES,N_SIG_M1,N_SIG
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP
  !----------------------------------------------------------------------------!

  !--Open NetCDF DATA FILE
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR READING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     STOP
  END IF

  !--Get Model Dimensions
  N_ELEMS   = GETDIM(NC_FID,LEN_TRIM('nele'),'nele')
  N_NODES   = GETDIM(NC_FID,LEN_TRIM('node'),'node')
  N_SIG_M1  = GETDIM(NC_FID,LEN_TRIM('siglay'),'siglay')
  N_SIG     = GETDIM(NC_FID,LEN_TRIM('siglev'),'siglev')

  M=N_NODES
  N=N_ELEMS
  KB=N_SIG    
  KBM1=N_SIG_M1
  KBM2=KB-2

  !--close file
  IERR = NF90_CLOSE(NC_FID)

  RETURN
END SUBROUTINE NCD_READ_GRID

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE NCD_READ_SHAPE(INFILE)
  !---------------------------------------------------------------------
  ! READ BATHYMETRY, SIGMA LEVELS AND GRID COEFFICIENTS IN A NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_NCD
  USE ALL_VARS
  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  CHARACTER(LEN=100), INTENT(IN) :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: K,I
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP

  !--OPEN NETCDF DATA FILE
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR READING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     STOP
  END IF

  !--Get Node Coordinates
  ALLOCATE(TEMP(M,1))
# if defined (SPHERICAL)
  CALL GETSVAR(NC_FID,LEN_TRIM('lon'),'lon',M,1,TEMP)
# else
  CALL GETSVAR(NC_FID,LEN_TRIM('x'),'x',M,1,TEMP)
# endif  
  VX(1:M) = TEMP(1:M,1)
# if defined (SPHERICAL)
  WHERE(VX < 0.0_SP) VX=360.0_SP+VX
# endif  
  DEALLOCATE(TEMP)

  ALLOCATE(TEMP(M,1))
# if defined (SPHERICAL)
  CALL GETSVAR(NC_FID,LEN_TRIM('lat'),'lat',M,1,TEMP)
# else
  CALL GETSVAR(NC_FID,LEN_TRIM('y'),'y',M,1,TEMP)
# endif  
  VY(1:M) = TEMP(1:M,1)
  DEALLOCATE(TEMP)

  !--Get Node Numbering
  ALLOCATE(TEMP(N,3))
  CALL GETSVAR(NC_FID,LEN_TRIM('nv'),'nv',N,3,TEMP)
  NV(1:N,1:3) = TEMP(1:N,1:3)
  DEALLOCATE(TEMP)
  NV(:,4) = NV(:,1)

  !--Get Bathymetry
  ALLOCATE(TEMP(N,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('h'),'h',M,1,TEMP)
  H(1:M) = TEMP(1:M,1)
  DEALLOCATE(TEMP)

  !--Get Sigma levels
  ALLOCATE(TEMP(M,KB))
  CALL GETSVAR(NC_FID,LEN_TRIM('siglev'),'siglev',M,KB,TEMP)
  Z(1:M,1:KB) = TEMP(1:M,1:KB)
  DEALLOCATE(TEMP)

  !--Compute derivative and intra-sigma levels
  DO K=1,KBM1
     ZZ(:,K)=0.5_SP*(Z(:,K)+Z(:,K+1))
     DZ(:,K)=Z(:,K)-Z(:,K+1)
  END DO
  ZZ(:,KB)=2.0_SP*ZZ(:,KBM1)-ZZ(:,KBM2)

  DO K=1,KBM2
     DZZ(:,K)=ZZ(:,K)-ZZ(:,K+1)
  END DO
  DZZ(:,KB-1)=0.0
  DZ(:,KB)=0.0

  DO I=1,N
    Z1(I,:)   = (Z(NV(I,1),:)+Z(NV(I,2),:)+Z(NV(I,3),:))/3.0
    ZZ1(I,:)  = (ZZ(NV(I,1),:)+ZZ(NV(I,2),:)+ZZ(NV(I,3),:))/3.0 
    DZ1(I,:)  = (DZ(NV(I,1),:)+DZ(NV(I,2),:)+DZ(NV(I,3),:))/3.0
    DZZ1(I,:) = (DZZ(NV(I,1),:)+DZZ(NV(I,2),:)+DZZ(NV(I,3),:))/3.0
  END DO
  
  !--Get Interpolation Parameters
  ALLOCATE(TEMP(N,4))
  CALL GETSVAR(NC_FID,LEN_TRIM('a1u'),'a1u',N,4,TEMP)
  A1U(1:N,:) = TEMP(1:N,:)
  DEALLOCATE(TEMP)

  ALLOCATE(TEMP(N,4))
  CALL GETSVAR(NC_FID,LEN_TRIM('a2u'),'a2u',N,4,TEMP)
  A2U(1:N,:) = TEMP(1:N,:)
  DEALLOCATE(TEMP)

  ALLOCATE(TEMP(N,3))
  CALL GETSVAR(NC_FID,LEN_TRIM('aw0'),'aw0',N,3,TEMP)
  AW0(1:N,:) = TEMP(1:N,:)
  DEALLOCATE(TEMP)

  ALLOCATE(TEMP(N,3))
  CALL GETSVAR(NC_FID,LEN_TRIM('awx'),'awx',N,3,TEMP)
  AWX(1:N,:) = TEMP(1:N,:)
  DEALLOCATE(TEMP)

  ALLOCATE(TEMP(N,3))
  CALL GETSVAR(NC_FID,LEN_TRIM('awy'),'awy',N,3,TEMP)
  AWY(1:N,:) = TEMP(1:N,:)
  DEALLOCATE(TEMP)

  !--Close file
  IERR = NF90_CLOSE(NC_FID)

  RETURN
END SUBROUTINE NCD_READ_SHAPE

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|


! **************MODIFIED BY JCH (04/07) TO ACQUIRE TIME FROM NETCDF FILE AT "LEVEL' OF READ

SUBROUTINE NCD_READ(INFILE,UL,VL,WWL,KHL,ELL,time,HO)
  !---------------------------------------------------------------------
  ! READ DATA FROM DAILY NETCDF FILES
  !---------------------------------------------------------------------

  USE MOD_NCD
  USE ALL_VARS, ONLY : NV
  USE LIMS
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  REAL(SP), DIMENSION(0:N,KB),INTENT(OUT)   :: UL,VL,WWL 
  REAL(SP), DIMENSION(0:M,KB),INTENT(OUT)   :: KHL 
  ! REAL(SP), DIMENSION(0:M,KB),INTENT(OUT) :: T1L
  REAL(SP), DIMENSION(0:M),INTENT(OUT)      :: ELL
  REAL(SP), INTENT(OUT)                     :: time
  INTEGER, INTENT(IN)                         :: HO
  CHARACTER(LEN=100), INTENT(IN)              :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: HT,I
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP,TEMP2
  INTEGER            :: N_TIMES


  !--Adjustement to read in Netcdf file
  HT=HO+1

  !--Open NetCDF Datafile 

!  **** MODIFIED 08/07 BY JHC TO RETRY OPENING THE INPUT FILE TWICE ON ERROR

  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR READING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     WRITE(*,*)' Second Try'
     IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
     IF(IERR /=NF90_NOERR)THEN
        WRITE(*,*)' 2nd ERROR READING ',TRIM(INFILE)
        WRITE(*,*)TRIM(NF90_STRERROR(IERR))
		WRITE(*,*)' Third Try'
        IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
        IF(IERR /=NF90_NOERR)THEN
           WRITE(*,*)' 3rd ERROR READING ',TRIM(INFILE)
           WRITE(*,*)TRIM(NF90_STRERROR(IERR))
		   WRITE(*,*)' Yer Out!'
		   STOP
		END IF
	  END IF
  END IF

  !---------------------------------------------------------------------
  ! Read Data from file INFILE at time level ht
  !---------------------------------------------------------------------

  ! ----- ADDED BY JHC (04/07) TO GET time (SECONDS AFTER START TIME) OF EACH READ------

  N_TIMES   = GETDIM(NC_FID,LEN_TRIM('time'),'time') ! number of times in file

  ALLOCATE(TEMP(HT,1))
  
  CALL GETSVAR(NC_FID,LEN_TRIM('time'),'time',HT,1,TEMP)
  
  
  time = TEMP(HT,1)
  DEALLOCATE(TEMP)

  !--free surface elevation
  ALLOCATE(TEMP(M,1))
  CALL GETDVAR(NC_FID,LEN_TRIM('zeta'),'zeta',M,1,TEMP,HT)
  ELL(1:M) = TEMP(1:M,1)
  DEALLOCATE(TEMP)

  !--salinity
  !ALLOCATE(TEMP(M,KBM1))
  !CALL GETDVAR(NC_FID,LEN_TRIM('salinity'),'salinity',M,KBM1,TEMP,HT)
  !S1L(1:M,1:KBM1) = TEMP(1:M,1:KBM1)
  !DEALLOCATE(TEMP)

  ! temperature
  ! ALLOCATE(TEMP(M,KBM1))
  ! CALL GETDVAR(NC_FID,LEN_TRIM('temp'),'temp',M,KBM1,TEMP,HT)
  ! T1L(1:M,1:KBM1) = TEMP(1:M,1:KBM1)
  ! DEALLOCATE(TEMP)


  !--U velocity 
  ALLOCATE(TEMP(N,KBM1))
  CALL GETDVAR(NC_FID,LEN_TRIM('u'),'u',N,KBM1,TEMP,HT)
  UL(1:N,1:KBM1) = TEMP(1:N,1:KBM1)
  DEALLOCATE(TEMP)

  !--V velocity
  ALLOCATE(TEMP(N,KBM1))
  CALL GETDVAR(NC_FID,LEN_TRIM('v'),'v',N,KBM1,TEMP,HT)
  VL(1:N,1:KBM1) = TEMP(1:N,1:KBM1)
  DEALLOCATE(TEMP)
!----------------------------here you have to check what name and format(node or cell) you output
  !--WW velocity
!  ALLOCATE(TEMP(N,KBM1))
!  CALL GETDVAR(NC_FID,LEN_TRIM('ww'),'ww',N,KBM1,TEMP,HT)
!  WWL(1:N,1:KBM1) = TEMP(1:N,1:KBM1)
!  DEALLOCATE(TEMP)
  !--W velocity (omega)
  ALLOCATE(TEMP(M,KBM1))
  ALLOCATE(TEMP2(N,KBM1))
!  CALL GETDVAR(NC_FID,LEN_TRIM('omega'),'omega',M,KBM1,TEMP,HT)  !original
 CALL GETDVAR(NC_FID,LEN_TRIM('wts'),'wts',M,KBM1,TEMP,HT) !pengfei some case people output as wts, maybe from cell(N) or from node(M)
  DO I = 1,N
    TEMP2(I,:) = ((TEMP(NV(I,1),:))+(TEMP(NV(I,2),:))+(TEMP(NV(I,3),:)))/3.0
  END DO  
  WWL(1:N,1:KBM1) = TEMP2(1:N,1:KBM1)
  DEALLOCATE(TEMP,TEMP2)

  !--KH
  ALLOCATE(TEMP(M,KBM1))
  CALL GETDVAR(NC_FID,LEN_TRIM('km'),'km',M,KBM1,TEMP,HT)
  KHL(1:M,1:KBM1) = TEMP(1:M,1:KBM1)
  DEALLOCATE(TEMP)

  !--Close file
  IERR = NF90_CLOSE(NC_FID)
  
  RETURN
END SUBROUTINE NCD_READ

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
! MODIFIED BY JHC 04/07 TO OUTPUT SURFACE ELEVATION (EP) AND BOTTOM DEPTH (HP)
!          Modified by JHC 07/07 to include the inwater variable in the output NCD file

!SUBROUTINE NCD_WRITE(INFILE,NPTS,TIME,LABEL,INDOMAIN,XP,YP,ZP,UP,VP,WP,NT)
SUBROUTINE NCD_WRITE(INFILE,NPTS,TIME,LABEL,INDOMAIN,XP,YP,ZP,UP,VP,WP,EP,HP,INWATER,NT)

  USE MOD_NCD
  IMPLICIT NONE
  !----------------------------------------------------------------------------! 
  INTEGER, INTENT(IN)                   :: NPTS, NT
  INTEGER, DIMENSION(NPTS),INTENT(IN)   :: LABEL,INDOMAIN,INWATER   ! 07/07 ADDED INWATER
  REAL(SP), DIMENSION(NPTS),INTENT(IN)  :: XP,YP,ZP,UP,VP,WP,EP,HP  ! 04/07 ADDED EP & HP
  REAL(SP), INTENT(IN)                  :: TIME
  CHARACTER(LEN=100), INTENT(IN)        :: INFILE
  !----------------------------------------------------------------------------! 
  CHARACTER(LEN=100)                  :: TSTRING,NETCDF_CONVENTION
  INTEGER                             :: IERR
  INTEGER                             :: TIME_DID,TIME_VID,NLAG_DID,IND_VID,LAB_VID
  INTEGER                             :: X_VID,Y_VID,Z_VID,U_VID,V_VID,W_VID
  INTEGER                             :: DYNMTIME(1)
  INTEGER                             :: STAT1D(1),DYNM1D(2)
  REAL(SP), ALLOCATABLE, DIMENSION(:) :: TEMP
  !----------------------------------------------------------------------------!

  IF(NT == 1)THEN

     !--Create File --> nc_fid (modified by JHC 4/07 to allow for input of file into MATLAB)
	 IERR = NF90_CREATE(TRIM(INFILE) ,NF90_CLOBBER, NC_FID)
     IF(IERR /=NF90_NOERR)THEN
        WRITE(*,*)'error creating',TRIM(INFILE)
        WRITE(*,*)TRIM(NF90_STRERROR(IERR))
        STOP
     END IF

     !--Netcdf Convention String
     NETCDF_CONVENTION = 'CF-1.0'

     !Global Attributes     
     TSTRING = "FVCOM Offline Lagrangian Particle Data"
     IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"title",TSTRING)
     IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"institution","SMAST") 
     IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"source"     ,"OFFLINE_FVCOM") 
     IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"modeler"    ,"PHIL MCCRACKEN") 
     IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"Conventions",TRIM(NETCDF_CONVENTION))

     !Dimensioning
     IERR = NF90_DEF_DIM(NC_FID,"nlag" ,NPTS,NLAG_DID)
     IERR = NF90_DEF_DIM(NC_FID,"time" ,NF90_UNLIMITED,TIME_DID)

     DYNMTIME     = (/TIME_DID/)                 !!Time
     DYNM1D       = (/NLAG_DID,TIME_DID/)        !!Dynamic 1d var 
     STAT1D       = (/NLAG_DID/)                 !!Static 1d var 
     !--Variable Definitions

     !!====time===============================================!
     IERR = NF90_DEF_VAR(NC_FID,"time",NF90_FLOAT,DYNMTIME,TIME_VID)
     IERR = NF90_PUT_ATT(NC_FID,TIME_VID,"long_name","time")
     IERR = NF90_PUT_ATT(NC_FID,TIME_VID,"units","seconds")

     !!=====label==============================================!
     IERR = NF90_DEF_VAR(NC_FID,"label",NF90_FLOAT,STAT1D,LAB_VID)
     IERR = NF90_PUT_ATT(NC_FID,LAB_VID,"long_name","particle label")
     IERR = NF90_PUT_ATT(NC_FID,LAB_VID,"units","")

     !!=====indomain=================================================!
     IERR = NF90_DEF_VAR(NC_FID,"indomain",NF90_FLOAT,DYNM1D,IND_VID)
     IERR = NF90_PUT_ATT(NC_FID,IND_VID,"long_name","particle indomain (1)")
     IERR = NF90_PUT_ATT(NC_FID,IND_VID,"units","")
     !!=====x=================================================!
     IERR = NF90_DEF_VAR(NC_FID,"x",NF90_FLOAT,DYNM1D,X_VID)
     IERR = NF90_PUT_ATT(NC_FID,X_VID,"long_name","particle x position")
     IERR = NF90_PUT_ATT(NC_FID,X_VID,"units","m")

     !!=====y=================================================!
     IERR = NF90_DEF_VAR(NC_FID,"y",NF90_FLOAT,DYNM1D,Y_VID)
     IERR = NF90_PUT_ATT(NC_FID,Y_VID,"long_name","particle y position")
     IERR = NF90_PUT_ATT(NC_FID,Y_VID,"units","m")

     !!=====z=================================================!
     IERR = NF90_DEF_VAR(NC_FID,"z",NF90_FLOAT,DYNM1D,Z_VID)
     IERR = NF90_PUT_ATT(NC_FID,Z_VID,"long_name","particle z position")
     IERR = NF90_PUT_ATT(NC_FID,Z_VID,"units","m")

     !!=====u=================================================!
     IERR = NF90_DEF_VAR(NC_FID,"u",NF90_FLOAT,DYNM1D,U_VID)
     IERR = NF90_PUT_ATT(NC_FID,U_VID,"long_name","particle u velocity")
     IERR = NF90_PUT_ATT(NC_FID,U_VID,"units","cm/s")

     !!=====v=================================================!
     IERR = NF90_DEF_VAR(NC_FID,"v",NF90_FLOAT,DYNM1D,V_VID)
     IERR = NF90_PUT_ATT(NC_FID,V_VID,"long_name","particle v velocity")
     IERR = NF90_PUT_ATT(NC_FID,V_VID,"units","cm/s")

     !!=====w=================================================!
     IERR = NF90_DEF_VAR(NC_FID,"omega",NF90_FLOAT,DYNM1D,W_VID)
     IERR = NF90_PUT_ATT(NC_FID,W_VID,"long_name","particle w velocity")
     IERR = NF90_PUT_ATT(NC_FID,W_VID,"units","mm/s")

 !   ADDED BY JCH 04/07 TO OUTPUT SURFACE ELEVATION AND BOTTOM DEPTH
     !!=====elev=============================================!
     IERR = NF90_DEF_VAR(NC_FID,"elev",NF90_FLOAT,DYNM1D,W_VID)
     IERR = NF90_PUT_ATT(NC_FID,W_VID,"long_name","surface elevation above particle")
     IERR = NF90_PUT_ATT(NC_FID,W_VID,"units","m")
      !!=====depth============================================!
     IERR = NF90_DEF_VAR(NC_FID,"depth",NF90_FLOAT,DYNM1D,W_VID)
     IERR = NF90_PUT_ATT(NC_FID,W_VID,"long_name","bottom depth at particle")
     IERR = NF90_PUT_ATT(NC_FID,W_VID,"units","m")

 !   ADDED BY JCH 07/07 TO OUTPUT INWATER	
     IERR = NF90_DEF_VAR(NC_FID,"inwater",NF90_FLOAT,DYNM1D,IND_VID)
     IERR = NF90_PUT_ATT(NC_FID,IND_VID,"long_name","particle inwater (1)")
     IERR = NF90_PUT_ATT(NC_FID,IND_VID,"units","")
  


 !   END ADDITION
     !--End definition section
     IERR = NF90_ENDDEF(NC_FID)
    
     !--Write Particle label
     ALLOCATE(TEMP(NPTS))
     TEMP(:) = FLOAT(LABEL(:))
     CALL PUTSVAR(NC_FID,LEN_TRIM('label'),'label',NPTS,TEMP)
     DEALLOCATE(TEMP)

     !--Close file
     IERR = NF90_CLOSE(NC_FID)

  ENDIF
  
  IERR = NF90_OPEN(TRIM(INFILE),NF90_WRITE,NC_FID)
  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR OPENING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     STOP
  END IF

  !write time to file
  ALLOCATE(TEMP(1))
  TEMP(1)=TIME
  CALL PUTDVAR(NC_FID,LEN_TRIM('time'),'time',1,TEMP,NT)
  DEALLOCATE(TEMP)

  !write indomain to file
  ALLOCATE(TEMP(NPTS))
  TEMP(:) = FLOAT(INDOMAIN(:))
  CALL PUTDVAR(NC_FID,LEN_TRIM('indomain'),'indomain',NPTS,TEMP,NT)
  DEALLOCATE(TEMP)

  !write position and velocity to file
  CALL PUTDVAR(NC_FID,LEN_TRIM('x'),'x',NPTS,XP,NT)
  CALL PUTDVAR(NC_FID,LEN_TRIM('y'),'y',NPTS,YP,NT)
  CALL PUTDVAR(NC_FID,LEN_TRIM('z'),'z',NPTS,ZP,NT)
  CALL PUTDVAR(NC_FID,LEN_TRIM('u'),'u',NPTS,UP,NT)
  CALL PUTDVAR(NC_FID,LEN_TRIM('v'),'v',NPTS,VP,NT)
  CALL PUTDVAR(NC_FID,LEN_TRIM('omega'),'omega',NPTS,WP,NT)

  ! ADDED BY JHC 04/07 TO WRITE SURFACE ELEVATION AND BOTTOM DEPTH TO FILE
  CALL PUTDVAR(NC_FID,LEN_TRIM('elev'),'elev',NPTS,EP,NT)
  CALL PUTDVAR(NC_FID,LEN_TRIM('depth'),'depth',NPTS,HP,NT)
  
  !   ADDED BY JCH 07/07 TO OUTPUT INWATER	
  ALLOCATE(TEMP(NPTS))
  TEMP(:) = FLOAT(INWATER(:))
  CALL PUTDVAR(NC_FID,LEN_TRIM('inwater'),'inwater',NPTS,TEMP,NT)
  DEALLOCATE(TEMP)

  ! END ADDITION

  !close file
  IERR = NF90_CLOSE(NC_FID)
!  write(*,'(i4)') IERR
END SUBROUTINE NCD_WRITE

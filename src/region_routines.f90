!> \file
!> \author Chris Bradley
!> \brief This module contains all region routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module contains all region routines.
MODULE REGION_ROUTINES

  USE BASE_ROUTINES
  USE COORDINATE_ROUTINES
  USE CMISS_CELLML
  USE DATA_POINT_ROUTINES
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES
  USE LISTS
  USE COMP_ENVIRONMENT
  USE CMISS_MPI
  USE MATHS




#ifndef NOMPIMOD
  USE MPI
#endif




#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  TYPE(REGIONS_TYPE) :: REGIONS
   !Interfaces

  INTERFACE REGION_LABEL_GET
    MODULE PROCEDURE REGION_LABEL_GET_C
    MODULE PROCEDURE REGION_LABEL_GET_VS
  END INTERFACE !REGION_LABEL_GET

  INTERFACE REGION_LABEL_SET
    MODULE PROCEDURE REGION_LABEL_SET_C
    MODULE PROCEDURE REGION_LABEL_SET_VS
  END INTERFACE !REGION_LABEL_SET

  PUBLIC REGION_COORDINATE_SYSTEM_GET,REGION_COORDINATE_SYSTEM_SET

  PUBLIC REGION_CREATE_START,REGION_CREATE_FINISH

  PUBLIC REGION_DATA_POINTS_GET

  PUBLIC REGION_DESTROY

  PUBLIC REGION_INITIALISE,REGION_FINALISE

  PUBLIC REGION_LABEL_GET,REGION_LABEL_SET

  PUBLIC REGION_NODES_GET

  PUBLIC REGION_USER_NUMBER_FIND, REGION_USER_NUMBER_TO_REGION

  PUBLIC REGIONS_INITIALISE,REGIONS_FINALISE

  PUBLIC COUPLED_DECOMPOSITION_CREATE_START, COUPLED_DECOMPOSITION_ADD_COUPLED_MESH, &
    & COUPLED_DECOMPOSITION_ADD_INTERFACE, COUPLED_DECOMPOSITION_CREATE_FINISH, &
    & DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD, COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION, &
    & COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION

CONTAINS

  !  ================================================================================================================================
  !

  !>Returns the coordinate system of region. \see OPENCMISS::CMISSRegionCoordinateSystemGet
  SUBROUTINE REGION_COORDINATE_SYSTEM_GET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<On exit, the coordinate system for the specified region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_COORDINATE_SYSTEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          CALL FlagError("Coordinate system is already associated.",ERR,ERROR,*999)
        ELSE
          COORDINATE_SYSTEM=>REGION%COORDINATE_SYSTEM
        ENDIF
      ELSE
        CALL FlagError("Region has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_COORDINATE_SYSTEM_GET")
    RETURN
999 ERRORSEXITS("REGION_COORDINATE_SYSTEM_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_COORDINATE_SYSTEM_GET

  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of region.  \see OPENCMISS::CMISSRegionCoordinateSystemSet
  SUBROUTINE REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to set the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_COORDINATE_SYSTEM_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FlagError("Region has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
            REGION%COORDINATE_SYSTEM=>COORDINATE_SYSTEM
          ELSE
            CALL FlagError("Coordinate system has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Coordinate system is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_COORDINATE_SYSTEM_SET")
    RETURN
999 ERRORSEXITS("REGION_COORDINATE_SYSTEM_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_COORDINATE_SYSTEM_SET

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a region. \see OPENCMISS::CMISSRegionCreateFinish
  SUBROUTINE REGION_CREATE_FINISH(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FlagError("Region has already been finished.",ERR,ERROR,*999)
      ELSE
        REGION%REGION_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Region : ",REGION%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",REGION%LABEL,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%PARENT_REGION)) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",REGION%PARENT_REGION%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",REGION%PARENT_REGION%LABEL, &
          & ERR,ERROR,*999)
      ENDIF
    ENDIF

    EXITS("REGION_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("REGION_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation a new region number USER_NUMBER as a sub region to the given PARENT_REGION, initialises all
  !>variables and inherits the PARENT_REGIONS coordinate system. \see OPENCMISS::CMISSRegionCreateFinish
  !>Default values set for the REGION's attributes are:
  !>- COORDINATE_SYSTEM: parent coordinate system. See \ref COORDINATE_SYSTEM_TYPE
  !>- DATA_POINTS: null
  !>- NODES: null
  !>- MESHES: 0 mesh
  !>- FIELDS: 0 field
  !>- EQUATIONS_SETS: 0 equation set
  !>- PARENT_REGION: global region
  !>- NUMBER_OF_SUB_REGIONS: 0
  !>- SUB_REGIONS: 0 region
  SUBROUTINE REGION_CREATE_START(USER_NUMBER,PARENT_REGION,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to create
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer to the parent region
    TYPE(REGION_TYPE), POINTER :: REGION !<On exit, a pointer to the created region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,region_idx
    TYPE(REGION_TYPE), POINTER :: NEW_REGION
    TYPE(REGION_PTR_TYPE), POINTER :: NEW_SUB_REGIONS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_REGION)
    NULLIFY(NEW_SUB_REGIONS)

    ENTERS("REGION_CREATE_START",ERR,ERROR,*997)

    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,NEW_REGION,ERR,ERROR,*997)
    IF(ASSOCIATED(NEW_REGION)) THEN
      LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
        & " has already been created."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*997)
    ELSE
      IF(ASSOCIATED(REGION)) THEN
        CALL FlagError("Region is already associated.",ERR,ERROR,*997)
      ELSE
        NULLIFY(REGION)
        IF(ASSOCIATED(PARENT_REGION)) THEN
          IF(PARENT_REGION%REGION_FINISHED) THEN
            IF(ASSOCIATED(PARENT_REGION%COORDINATE_SYSTEM)) THEN
              !Initialise the region
              CALL REGION_INITIALISE(REGION,ERR,ERROR,*999)
              !Set the user number
              REGION%USER_NUMBER=USER_NUMBER
              !CPB 21/02/07 The vstring operation crashes the AIX compiler so put a CHAR() etc. around it.
              !REGION%LABEL="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
              LOCAL_STRING="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
              REGION%LABEL=CHAR(LOCAL_STRING)
              IF(ERR/=0) GOTO 999
              REGION%COORDINATE_SYSTEM=>PARENT_REGION%COORDINATE_SYSTEM
              !Adjust the parent region to include this new daughter
              ALLOCATE(NEW_SUB_REGIONS(PARENT_REGION%NUMBER_OF_SUB_REGIONS+1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new sub-regions.",ERR,ERROR,*999)
              DO region_idx=1,PARENT_REGION%NUMBER_OF_SUB_REGIONS
                NEW_SUB_REGIONS(region_idx)%PTR=>PARENT_REGION%SUB_REGIONS(region_idx)%PTR
              ENDDO !region_no
              PARENT_REGION%NUMBER_OF_SUB_REGIONS=PARENT_REGION%NUMBER_OF_SUB_REGIONS+1
              NEW_SUB_REGIONS(PARENT_REGION%NUMBER_OF_SUB_REGIONS)%PTR=>REGION
              IF(ASSOCIATED(PARENT_REGION%SUB_REGIONS)) DEALLOCATE(PARENT_REGION%SUB_REGIONS)
              PARENT_REGION%SUB_REGIONS=>NEW_SUB_REGIONS
              !Set the new regions parent region to the parent region
              REGION%PARENT_REGION=>PARENT_REGION
            ELSE
              CALL FlagError("Parent region does not have an associated coordinate system.",ERR,ERROR,*997)
            ENDIF
          ELSE
            CALL FlagError("Parent region has not been finished.",ERR,ERROR,*997)
          ENDIF
        ELSE
          CALL FlagError("Parent region is not associated.",ERR,ERROR,*997)
        ENDIF
      ENDIF
    ENDIF

    EXITS("REGION_CREATE_START")
    RETURN
999 CALL REGION_FINALISE(REGION,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_SUB_REGIONS)) DEALLOCATE(NEW_SUB_REGIONS)
997 ERRORSEXITS("REGION_CREATE_START",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data points for a region. \see OPENCMISS::CMISSRegionDataPointsGet
  SUBROUTINE REGION_DATA_POINTS_GET(REGION,DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the data points for
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<On exit, a pointer to the data points for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_DATA_POINTS_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        IF(ASSOCIATED(DATA_POINTS)) THEN
          CALL FlagError("Data points is already associated.",ERR,ERROR,*998)
        ELSE
          DATA_POINTS=>REGION%DATA_POINTS
          IF(.NOT.ASSOCIATED(DATA_POINTS)) CALL FlagError("Data points is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Region has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*998)
    ENDIF

    EXITS("REGION_DATA_POINTS_GET")
    RETURN
999 NULLIFY(DATA_POINTS)
998 ERRORSEXITS("REGION_DATA_POINTS_GET",ERR,ERROR)
    RETURN 1

  END SUBROUTINE REGION_DATA_POINTS_GET


  !
  !================================================================================================================================
  !

  !>Destroys a region given by USER_NUMBER and all sub-regions under it. \todo create destroy by pointer method. \see OPENCMISS::CMISSRegionDestroy
  RECURSIVE SUBROUTINE REGION_DESTROY_NUMBER(USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: count,nr
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(REGION_PTR_TYPE), POINTER :: NEW_SUB_REGIONS(:)

    ENTERS("REGION_DESTROY_NUMBER",ERR,ERROR,*999)

    NULLIFY(REGION)
    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN

!!NOTE: We have to find a pointer to the region to destroy within this routine rather than passing in a pointer to a
!!DESTROY_REGION_PTR type routine because we need to change REGION%SUB_REGIONS of the PARENT region and this would violate section
!!12.4.1.6 of the Fortran standard if the dummy REGION pointer argument was associated with the SUB_REGIONS(x)%PTR actual
!!argument.

      IF(REGION%NUMBER_OF_SUB_REGIONS==0) THEN
        !No more daughter sub regions so delete this instance
        IF(ASSOCIATED(REGION%PARENT_REGION)) THEN
          NULLIFY(NEW_SUB_REGIONS)
          IF(REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS>1) THEN
            !If the parent region has more than one sub regions then remove this instance from its sub-regions list
            ALLOCATE(NEW_SUB_REGIONS(REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new sub-regions.",ERR,ERROR,*999)
            count=0
            DO nr=1,REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS
              IF(REGION%PARENT_REGION%SUB_REGIONS(nr)%PTR%USER_NUMBER/=REGION%USER_NUMBER) THEN
                count=count+1
                NEW_SUB_REGIONS(count)%PTR=>REGION%PARENT_REGION%SUB_REGIONS(nr)%PTR
              ENDIF
            ENDDO !nr
          ENDIF
          REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS=REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS-1
          IF(ASSOCIATED(REGION%PARENT_REGION%SUB_REGIONS)) DEALLOCATE(REGION%PARENT_REGION%SUB_REGIONS)
          REGION%PARENT_REGION%SUB_REGIONS=>NEW_SUB_REGIONS
          !Finalise the region
          CALL REGION_FINALISE(REGION,ERR,ERROR,*999)
        ELSE
          CALL FlagError("Parent region is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !Recursively delete sub regions first
        DO WHILE(REGION%NUMBER_OF_SUB_REGIONS>0)
          CALL REGION_DESTROY_NUMBER(REGION%SUB_REGIONS(1)%PTR%USER_NUMBER,ERR,ERROR,*999)
        ENDDO
        !Now delete this instance
        CALL REGION_DESTROY_NUMBER(REGION%USER_NUMBER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Region number does not exist.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_DESTROY_NUMBER")
    RETURN
999 ERRORSEXITS("REGION_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a region identified by a pointer and all sub-regions under it. \see OPENCMISS::CMISSRegionDestroy
  SUBROUTINE REGION_DESTROY(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: USER_NUMBER

    ENTERS("REGION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      USER_NUMBER=REGION%USER_NUMBER
      CALL REGION_DESTROY_NUMBER(USER_NUMBER,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_DESTROY")
    RETURN
999 ERRORSEXITS("REGION_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a region and deallocates all memory
  SUBROUTINE REGION_FINALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      REGION%LABEL=""
      CALL CELLML_ENVIRONMENTS_FINALISE(REGION%CELLML_ENVIRONMENTS,ERR,ERROR,*999)
      CALL EQUATIONS_SETS_FINALISE(REGION,ERR,ERROR,*999)
      CALL FIELDS_FINALISE(REGION%FIELDS,ERR,ERROR,*999)
      CALL MESHES_FINALISE(REGION%MESHES,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%DATA_POINTS)) CALL DATA_POINTS_DESTROY(REGION%DATA_POINTS,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%NODES)) CALL NODES_DESTROY(REGION%NODES,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%SUB_REGIONS)) DEALLOCATE(REGION%SUB_REGIONS)
      IF(ASSOCIATED(REGION%INTERFACES)) CALL INTERFACES_FINALISE(REGION%INTERFACES,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%GENERATED_MESHES)) CALL GENERATED_MESHES_FINALISE(REGION%GENERATED_MESHES,ERR,ERROR,*999)
      DEALLOCATE(REGION)
    ENDIF

    EXITS("REGION_FINALISE")
    RETURN
999 ERRORSEXITS("REGION_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a region.
  SUBROUTINE REGION_INITIALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("REGION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      CALL FlagError("Region is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(REGION,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate region.",ERR,ERROR,*999)
      REGION%USER_NUMBER=0
      REGION%REGION_FINISHED=.FALSE.
      REGION%LABEL=""
      NULLIFY(REGION%COORDINATE_SYSTEM)
      NULLIFY(REGION%DATA_POINTS)
      NULLIFY(REGION%NODES)
      NULLIFY(REGION%MESHES)
      NULLIFY(REGION%GENERATED_MESHES)
      NULLIFY(REGION%FIELDS)
      NULLIFY(REGION%EQUATIONS_SETS)
      NULLIFY(REGION%CELLML_ENVIRONMENTS)
      NULLIFY(REGION%PARENT_REGION)
      REGION%NUMBER_OF_SUB_REGIONS=0
      NULLIFY(REGION%SUB_REGIONS)
      NULLIFY(REGION%INTERFACES)
      CALL MESHES_INITIALISE(REGION,ERR,ERROR,*999)
      CALL GENERATED_MESHES_INITIALISE(REGION,ERR,ERROR,*999)
      CALL FIELDS_INITIALISE(REGION,ERR,ERROR,*999)
      CALL EQUATIONS_SETS_INITIALISE(REGION,ERR,ERROR,*999)
      CALL CELLML_ENVIRONMENTS_INITIALISE(REGION,ERR,ERROR,*999)
      CALL INTERFACES_INITIALISE(REGION,ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_INITIALISE")
    RETURN
999 CALL REGION_FINALISE(REGION,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("REGION_INITIALISE",ERR,ERROR)
    RETURN 1

  END SUBROUTINE REGION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the label of a region. \see OPENCMISS::CMISSRegionLabelGet
  SUBROUTINE REGION_LABEL_GET_C(REGION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On return the region label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: C_LENGTH,VS_LENGTH

    ENTERS("REGION_LABEL_GET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      C_LENGTH=LEN(LABEL)
      VS_LENGTH=LEN_TRIM(REGION%LABEL)
      IF(C_LENGTH>VS_LENGTH) THEN
        LABEL=CHAR(REGION%LABEL,VS_LENGTH)
      ELSE
        LABEL=CHAR(REGION%LABEL,C_LENGTH)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_LABEL_GET_C")
    RETURN
999 ERRORSEXITS("REGION_LABEL_GET_C",ERR,ERROR)
    RETURN 1

  END SUBROUTINE REGION_LABEL_GET_C

   !
  !================================================================================================================================
  !

  !>Returns the label of a region. \see OPENCMISS::CMISSRegionLabelGet
  SUBROUTINE REGION_LABEL_GET_VS(REGION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On return the region label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_LABEL_GET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      !CPB 20/2/07 The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
      LABEL=VAR_STR(CHAR(REGION%LABEL))
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_LABEL_GET_VS")
    RETURN
999 ERRORSEXITS("REGION_LABEL_GET_VS",ERR,ERROR)
    RETURN 1

  END SUBROUTINE REGION_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Sets the label of a region. \see OPENCMISS::CMISSRegionLabelSet
  SUBROUTINE REGION_LABEL_SET_C(REGION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_LABEL_SET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FlagError("Region has been finished.",ERR,ERROR,*999)
      ELSE
        REGION%LABEL=LABEL
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_LABEL_SET_C")
    RETURN
999 ERRORSEXITS("REGION_LABEL_SET_C",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_LABEL_SET_C

  !
  !================================================================================================================================
  !

  !>Sets the label of a region. \see OPENCMISS::CMISSRegionLabelSet
  SUBROUTINE REGION_LABEL_SET_VS(REGION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_LABEL_SET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FlagError("Region has been finished.",ERR,ERROR,*999)
      ELSE
        REGION%LABEL=LABEL
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_LABEL_SET_VS")
    RETURN
999 ERRORSEXITS("REGION_LABEL_SET_VS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_LABEL_SET_VS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the nodes for a region. \see OPENCMISS::CMISSRegionNodesGet
  SUBROUTINE REGION_NODES_GET(REGION,NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the nodes for
    TYPE(NODES_TYPE), POINTER :: NODES !<On exit, a pointer to the nodes for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGION_NODES_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        IF(ASSOCIATED(NODES)) THEN
          CALL FlagError("Nodes is already associated.",ERR,ERROR,*998)
        ELSE
          NODES=>REGION%NODES
          IF(.NOT.ASSOCIATED(NODES)) CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Region has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*998)
    ENDIF

    EXITS("REGION_NODES_GET")
    RETURN
999 NULLIFY(NODES)
998 ERRORSEXITS("REGION_NODES_GET",ERR,ERROR)
    RETURN 1

  END SUBROUTINE REGION_NODES_GET

  !
  !================================================================================================================================
  !

  !>Finds and returns in REGION a pointer to the region with the number given in USER_NUMBER. If no region with that number
  !>exits REGION is left nullified.
  SUBROUTINE REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*)

     !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to find
    TYPE(REGION_TYPE), POINTER :: REGION !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nr
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION

    ENTERS("REGION_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      CALL FlagError("Region is already associated.",ERR,ERROR,*999)
    ELSE
      NULLIFY(REGION)
      WORLD_REGION=>REGIONS%WORLD_REGION
      IF(ASSOCIATED(WORLD_REGION)) THEN
        IF(USER_NUMBER==0) THEN
          REGION=>WORLD_REGION
        ELSE
          nr=1
          DO WHILE(nr<=WORLD_REGION%NUMBER_OF_SUB_REGIONS.AND..NOT.ASSOCIATED(REGION))
            CALL REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,WORLD_REGION%SUB_REGIONS(nr)%PTR,REGION,ERR,ERROR,*999)
            IF(.NOT.ASSOCIATED(REGION)) nr=nr+1
          END DO
        ENDIF
      ELSE
        CALL FlagError("World region is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF

    EXITS("REGION_USER_NUMBER_FIND")
    RETURN
999 ERRORSEXITS("REGION_USER_NUMBER_FIND",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finds and returns in REGION a pointer to the region with the number given in USER_NUMBER starting from the
  !>START_REGION and searching all sub-regions under the START_REGION. If no region with that number exit REGION is
  !>left nullified.
  RECURSIVE SUBROUTINE REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,START_REGION,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find
    TYPE(REGION_TYPE), POINTER :: START_REGION !<A pointer to the region to start the search from
    TYPE(REGION_TYPE), POINTER :: REGION !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nr

    ENTERS("REGION_USER_NUMBER_FIND_PTR",ERR,ERROR,*999)

    NULLIFY(REGION)
    IF(ASSOCIATED(START_REGION)) THEN
      IF(START_REGION%USER_NUMBER==USER_NUMBER) THEN
        REGION=>START_REGION
      ELSE
        nr=1
        DO WHILE(nr<=START_REGION%NUMBER_OF_SUB_REGIONS.AND..NOT.ASSOCIATED(REGION))
          CALL REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,REGION,START_REGION%SUB_REGIONS(nr)%PTR,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(REGION)) nr=nr+1
        END DO
      ENDIF
    ELSE
      CALL FlagError("Start region is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("REGION_USER_NUMBER_FIND_PTR")
    RETURN
999 ERRORSEXITS("REGION_USER_NUMBER_FIND_PTR",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGION_USER_NUMBER_FIND_PTR

  !
  !================================================================================================================================
  !

  !>Finalises the regions and destroys any current regions.
  SUBROUTINE REGIONS_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("REGIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGIONS%WORLD_REGION)) THEN
      !Destroy any global region daughter regions first
      DO WHILE(REGIONS%WORLD_REGION%NUMBER_OF_SUB_REGIONS>0)
        CALL REGION_DESTROY_NUMBER(REGIONS%WORLD_REGION%SUB_REGIONS(1)%PTR%USER_NUMBER,ERR,ERROR,*999)
      ENDDO !region
      !Destroy global region and deallocated any memory allocated in the global region
      CALL REGION_FINALISE(REGIONS%WORLD_REGION,ERR,ERROR,*999)
      NULLIFY(REGIONS%WORLD_REGION)
    ENDIF

    EXITS("REGIONS_FINALISE")
    RETURN
999 ERRORSEXITS("REGIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the regions and creates the global world region.
  SUBROUTINE REGIONS_INITIALISE(WORLD_REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION !<On exit, a pointer to the world region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: WORLD_COORDINATE_SYSTEM

    NULLIFY(WORLD_COORDINATE_SYSTEM)

    ENTERS("REGIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(WORLD_REGION)) THEN
      CALL FlagError("World region is already associated.",ERR,ERROR,*999)
    ELSE
      CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(0,WORLD_COORDINATE_SYSTEM,ERR,ERROR,*999)
      IF(ASSOCIATED(WORLD_COORDINATE_SYSTEM)) THEN
        CALL REGION_INITIALISE(REGIONS%WORLD_REGION,ERR,ERROR,*999)
        REGIONS%WORLD_REGION%USER_NUMBER=0
        REGIONS%WORLD_REGION%LABEL="World Region"
        REGIONS%WORLD_REGION%COORDINATE_SYSTEM=>WORLD_COORDINATE_SYSTEM
        REGIONS%WORLD_REGION%REGION_FINISHED=.TRUE.
        !Return the pointer
        WORLD_REGION=>REGIONS%WORLD_REGION
      ELSE
        CALL FlagError("World coordinate system has not been created.",ERR,ERROR,*999)
      ENDIF
    ENDIF

    EXITS("REGIONS_INITIALISE")
    RETURN
999 ERRORSEXITS("REGIONS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE REGIONS_INITIALISE

  !
  !================================================================================================================================
  !

  !> Find the region with the given user number, or throw an error if it does not exist.
  SUBROUTINE REGION_USER_NUMBER_TO_REGION( USER_NUMBER, REGION, ERR, ERROR, * )
    !Arguments
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to find
    TYPE(REGION_TYPE), POINTER :: REGION !<On return, a pointer to the region with the specified user number.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("REGION_USER_NUMBER_TO_REGION", ERR, ERROR, *999 )

    NULLIFY( REGION )
    CALL REGION_USER_NUMBER_FIND( USER_NUMBER, REGION, ERR, ERROR, *999 )
    IF( .NOT.ASSOCIATED( REGION ) ) THEN
      LOCAL_ERROR = "A region with an user number of "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*", ERR, ERROR ) )//" does not exist."
      CALL FlagError( LOCAL_ERROR, ERR, ERROR, *999 )
    ENDIF

    EXITS( "REGION_USER_NUMBER_TO_REGION" )
    RETURN
999 ERRORSEXITS( "REGION_USER_NUMBER_TO_REGION", ERR, ERROR )
    RETURN 1

  END SUBROUTINE REGION_USER_NUMBER_TO_REGION

 !=====================================================================================================
 !> The following subroutine calculate the maximum edge length in the interface graph G_I
  SUBROUTINE GET_MAXIMUM_EDGE_LENGTH(COUPLED_DECOMPOSITION, MAXIMUM_EDGE_LENGTH, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN)     :: COUPLED_DECOMPOSITION !<A pointer to the decomposition to calculate the node domains for.
    REAL(RP), INTENT(OUT)                     :: MAXIMUM_EDGE_LENGTH   !Maximum edge length in a mesh given by object MESH.
    INTEGER(INTG), INTENT(OUT)                :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)         :: ERROR !<The error string

    !Local Variables
    REAL(RP), ALLOCATABLE                     :: EDGE_LENGTHS(:)
    INTEGER(INTG)                             :: ADJACENCY_NODE, adjacency_idx, &
      & edge_length_idx, MY_COMPUTATIONAL_NODE, node_idx
    REAL(RP)                                  :: VECTOR_ALONG_THE_EDGE(3)
    TYPE(MESH_TYPE), POINTER                  :: INTERFACE_MESH

    ENTERS("GET_MAXIMUM_EDGE_LENGTH",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN

      INTERFACE_MESH=>&
        COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH

      MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      IF(ERR/=0) GOTO 999


      edge_length_idx = 0
      ALLOCATE(EDGE_LENGTHS&
        &(INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes*INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes), STAT=ERR)

      IF(ERR/=0) CALL FlagError("Could not allocate EDGE_LENGTHS array",ERR,ERROR,*999)
      EDGE_LENGTHS =  0

        DO node_idx = 1, SIZE(INTERFACE_MESH%TOPOLOGY(1)%ptr%Nodes%Nodes,1)

           DO adjacency_idx = 1, SIZE(INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNodes)


             ADJACENCY_NODE = INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNodes(adjacency_idx)

             VECTOR_ALONG_THE_EDGE(1) =  COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(node_idx,1) - &
               & COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(ADJACENCY_NODE,1)
             VECTOR_ALONG_THE_EDGE(2) =  COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(node_idx,2) - &
               & COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(ADJACENCY_NODE,2)
             VECTOR_ALONG_THE_EDGE(3) =  COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(node_idx,3) - &
               & COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(ADJACENCY_NODE,3)

             edge_length_idx          =  edge_length_idx + 1

             EDGE_LENGTHS(edge_length_idx)   = NORM2(VECTOR_ALONG_THE_EDGE)

           END DO !adjncy_idx
         END DO !node_idx

       MAXIMUM_EDGE_LENGTH = MAXVAL(EDGE_LENGTHS)

    ELSE
     CALL FlagError("COupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    RETURN
999 RETURN 1

  END SUBROUTINE GET_MAXIMUM_EDGE_LENGTH

 !================================================================================================================
   !> The following subroutine builds interedges I_{iI} between coupled mesh graphs G_i and the interface graphs G_I

  SUBROUTINE COUPLED_DECOMPOSITION_GET_INTER_EDGES(COUPLED_DECOMPOSITION, MAXIMUM_INTERFACE_EDGE_LENGTH, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)   :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    REAL(RP), INTENT(IN)                                       :: MAXIMUM_INTERFACE_EDGE_LENGTH !<Maximum edge length in interface graph mesh G^I.
    INTEGER(INTG), INTENT(OUT)                                 :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                          :: ERROR !<The error string.

    !Local Variables
    TYPE(LIST_TYPE),  POINTER :: BOUNDARY_NODE_LIST, INTER_EDGE_LIST_INTERFACE_MESH, INTER_EDGE_LIST_COUPLED_MESH
    TYPE(LIST_PTR_TYPE), ALLOCATABLE:: ADJNCY_RESTRICTED_TO_INTERFACE_LIST(:)
    TYPE(MESH_TYPE), POINTER:: COUPLED_MESH, INTERFACE_MESH
    INTEGER(INTG) :: boundary_node_idx,coupled_mesh_normal_vector_idx, interface_mesh_normal_vector_idx, &
      & coupled_mesh_node_idx, distance_idx,INTERFACE_MESH_NODE_TO_ADD,interface_node_idx, &
      & inter_edge_idx,MY_COMPUTATIONAL_NODE,NUMBER_OF_BOUNDARY_NODES,NUMBER_OF_INTER_EDGES, &
      & NUMBER_OF_NODES,surrounding_node_idx
    REAL(RP) :: DIFFERENCE_IN_ANGLE
    REAL(RP), ALLOCATABLE :: DISTANCE(:), NORMAL_VECTOR_COUPLED_MESH(:,:), NORMAL_VECTOR_INTERFACE_MESH(:,:)
    INTEGER(INTG), ALLOCATABLE  :: ADJNCY_RESTRICTED_TO_INTERFACE(:,:),BOUNDARY_NODES(:), &
      & COUPLED_MESH_INTER_EDGE_NODES(:),INTERFACE_MESH_INTER_EDGE_NODES(:),LIST_OF_NODES(:), &
      & NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE(:)
    LOGICAL :: INTERFACE_NODE_DETECTED

   ENTERS("COUPLED_DECOMPOSITION_GET_INTER_EDGES",ERR,ERROR,*999)

   IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
     COUPLED_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%MESH
     IF(ASSOCIATED(COUPLED_MESH)) THEN
       INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH
       IF(ASSOCIATED(INTERFACE_MESH)) THEN

         MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
         IF(ERR/=0) GOTO 999

         NULLIFY(INTER_EDGE_LIST_INTERFACE_MESH)
         CALL LIST_CREATE_START(INTER_EDGE_LIST_INTERFACE_MESH,ERR,ERROR,*999)
         CALL LIST_DATA_TYPE_SET(INTER_EDGE_LIST_INTERFACE_MESH,LIST_INTG_TYPE,ERR,ERROR,*999)
         CALL LIST_INITIAL_SIZE_SET(INTER_EDGE_LIST_INTERFACE_MESH,COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes,ERR,ERROR,*999)
         CALL LIST_CREATE_FINISH(INTER_EDGE_LIST_INTERFACE_MESH,ERR,ERROR,*999)

         NULLIFY(INTER_EDGE_LIST_COUPLED_MESH)
         CALL LIST_CREATE_START(INTER_EDGE_LIST_COUPLED_MESH,ERR,ERROR,*999)
         CALL LIST_DATA_TYPE_SET(INTER_EDGE_LIST_COUPLED_MESH,LIST_INTG_TYPE,ERR,ERROR,*999)
         CALL LIST_INITIAL_SIZE_SET(INTER_EDGE_LIST_COUPLED_MESH,COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes,ERR,ERROR,*999)
         CALL LIST_CREATE_FINISH(INTER_EDGE_LIST_COUPLED_MESH,ERR,ERROR,*999)

         NULLIFY(BOUNDARY_NODE_LIST)
         CALL LIST_CREATE_START(BOUNDARY_NODE_LIST,ERR,ERROR,*999)
         CALL LIST_DATA_TYPE_SET(BOUNDARY_NODE_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
         CALL LIST_INITIAL_SIZE_SET(BOUNDARY_NODE_LIST,COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes,ERR,ERROR,*999)
         CALL LIST_CREATE_FINISH(BOUNDARY_NODE_LIST,ERR,ERROR,*999)

         !populate boundary nodes
         DO coupled_mesh_node_idx = 1 , COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes
           IF(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(coupled_mesh_node_idx)%boundaryNode) THEN
             CALL LIST_ITEM_ADD(BOUNDARY_NODE_LIST, coupled_mesh_node_idx, ERR, ERROR, *999)
           END IF
         END DO ! coupled_mesh_node_idx

         CALL LIST_REMOVE_DUPLICATES(BOUNDARY_NODE_LIST, ERR, ERROR, *999)
         CALL LIST_DETACH_AND_DESTROY(BOUNDARY_NODE_LIST, NUMBER_OF_BOUNDARY_NODES, BOUNDARY_NODES, ERR,ERROR,*999)

         ALLOCATE(NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE(NUMBER_OF_BOUNDARY_NODES), STAT=ERR)
         IF(ERR/=0) CALL FlagError("Could not allocate ADJNCY_RESTRICTED_TO_INTERFACE array",ERR,ERROR,*999)

         ALLOCATE(ADJNCY_RESTRICTED_TO_INTERFACE(NUMBER_OF_BOUNDARY_NODES, &
           & COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes),STAT=ERR)
         IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE array",ERR,ERROR,*999)

         ALLOCATE(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(NUMBER_OF_BOUNDARY_NODES), STAT=ERR)
         IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE_LIST array",ERR,ERROR,*999)


         DO coupled_mesh_node_idx = 1 , NUMBER_OF_BOUNDARY_NODES

           boundary_node_idx      =   BOUNDARY_NODES(coupled_mesh_node_idx)

           INTERFACE_NODE_DETECTED = .FALSE.
           NULLIFY(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR)

           CALL LIST_CREATE_START(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR,ERR,ERROR,*999)

           CALL LIST_DATA_TYPE_SET(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR,LIST_INTG_TYPE, ERR, ERROR, *999)

           CALL LIST_INITIAL_SIZE_SET(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR, &
             & SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(boundary_node_idx)%surroundingNodes),ERR,ERROR,*999)

           CALL LIST_CREATE_FINISH(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR,ERR, ERROR, *999)

           ! To determine the size of COupledMeshNodeAdjancyRestrictedToBoundary
           DO surrounding_node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(boundary_node_idx)%surroundingNodes)


             IF(ANY(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(boundary_node_idx)%surroundingNodes(surrounding_node_idx)== &
               BOUNDARY_NODES)) THEN
               ! Node adjancy containing only the nodes that belong to the boundary.

               CALL LIST_ITEM_ADD(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR, &
                 & COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(boundary_node_idx)%surroundingNodes(surrounding_node_idx), &
                   & ERR,ERROR,*999)

             END  IF  !  IF( any(CoupledMesh%adjancy(COupledMesh%BOundaryNodes ...


           END DO !surrounding_node_idx

           CALL LIST_REMOVE_DUPLICATES(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR,ERR,ERROR,*999)

           CALL LIST_DETACH_AND_DESTROY(ADJNCY_RESTRICTED_TO_INTERFACE_LIST(coupled_mesh_node_idx)%PTR, &
             & NUMBER_OF_NODES,LIST_OF_NODES,ERR,ERROR,*999)


           NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE(coupled_mesh_node_idx:coupled_mesh_node_idx)=NUMBER_OF_NODES

           ADJNCY_RESTRICTED_TO_INTERFACE( &
             & coupled_mesh_node_idx,1:NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE(coupled_mesh_node_idx)) = &
               & LIST_OF_NODES(1:NUMBER_OF_NODES)
           DEALLOCATE(LIST_OF_NODES)

         END DO !COUPLED_mesh_idx

         ALLOCATE(DISTANCE(INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes), STAT=ERR)
         IF(ERR /= 0) CALL FlagError(" Unable to allocate DISTANCE array.",ERR,ERROR,*999)


         DO coupled_mesh_node_idx = 1, NUMBER_OF_BOUNDARY_NODES

           distance_idx  = 0

           boundary_node_idx      =   BOUNDARY_NODES(coupled_mesh_node_idx)

           DISTANCE =  -1

           DO interface_node_idx = 1, INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes

             DISTANCE(interface_node_idx) = NORM2(COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(interface_node_idx,:) &
               & - COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES(boundary_node_idx,:))

           END DO ! interface_node_idx

           IF(MINVAL(PACK(DISTANCE,DISTANCE >= -1E-5)) < MAXIMUM_INTERFACE_EDGE_LENGTH ) THEN

             INTERFACE_MESH_NODE_TO_ADD= INT(SQRT(REAL(DOT_PRODUCT(&
               & MINLOC(PACK(DISTANCE,DISTANCE >= -1E-5)-MINVAL(PACK(DISTANCE,DISTANCE >= -1E-5))),&
                 & MINLOC(PACK(DISTANCE,DISTANCE >= -1E-5)-MINVAL(PACK(DISTANCE,DISTANCE >= -1E-5)))))),INTG)

             ! Calculate normal at vertex coupled_mesh_node_idx.
             CALL COUPLED_MESH_CALCULATE_NORMAL(COUPLED_DECOMPOSITION,  &
               & boundary_node_idx, NORMAL_VECTOR_COUPLED_MESH, ERR, ERROR, *999)

             ! Calculate normal at vertex interface_mesh_node_idx.
             CALL INTERFACE_MESH_CALCULATE_NORMAL(COUPLED_DECOMPOSITION,  &
               & INTERFACE_MESH_NODE_TO_ADD, NORMAL_VECTOR_INTERFACE_MESH, ERR, ERROR, *999)

                      DO coupled_mesh_normal_vector_idx = 1, SIZE(NORMAL_VECTOR_COUPLED_MESH,2)
               DO interface_mesh_normal_vector_idx = 1, SIZE(NORMAL_VECTOR_INTERFACE_MESH,2)

                ! COS(\Theeta)  =  (U.V)/( |U|.|V|)
                 DIFFERENCE_IN_ANGLE =  ACOS(DOT_PRODUCT(NORMAL_VECTOR_COUPLED_MESH(:,coupled_mesh_normal_vector_idx), &
                   & NORMAL_VECTOR_INTERFACE_MESH(:,interface_mesh_normal_vector_idx))/ &
                     & ((NORM2(NORMAL_VECTOR_COUPLED_MESH(:,coupled_mesh_normal_vector_idx))* &
                       & NORM2(NORMAL_VECTOR_INTERFACE_MESH(:,interface_mesh_normal_vector_idx)))))

                 DIFFERENCE_IN_ANGLE = DIFFERENCE_IN_ANGLE*180/3.14159 ! Calculating the difference of angle in degrees

                 IF((DIFFERENCE_IN_ANGLE .GE. 0 .AND. DIFFERENCE_IN_ANGLE .LE. 5.0001) .OR. &
                   & (DIFFERENCE_IN_ANGLE .GE. 174.9999 .AND. DIFFERENCE_IN_ANGLE .LE. 181.0001)) THEN !With tolerance of 1E-4


                   CALL LIST_ITEM_ADD(INTER_EDGE_LIST_INTERFACE_MESH,INTERFACE_MESH_NODE_TO_ADD, ERR, ERROR, *999)

                   CALL LIST_ITEM_ADD(INTER_EDGE_LIST_COUPLED_MESH,boundary_node_idx, ERR, ERROR, *999)

                   EXIT

                 END IF

               END DO !interface_mesh_normal_vector_idx

                 IF((DIFFERENCE_IN_ANGLE .GE. 0 .AND. DIFFERENCE_IN_ANGLE .LE. 5.0001) .OR. &
                   & (DIFFERENCE_IN_ANGLE .GE. 174.9999 .AND. DIFFERENCE_IN_ANGLE .LE. 181.0001)) EXIT !With tolerance of 1E-4

             END DO ! coupled_mesh_normal_vector_idx

             IF(ALLOCATED(NORMAL_VECTOR_COUPLED_MESH)) DEALLOCATE(NORMAL_VECTOR_COUPLED_MESH)
             IF(ALLOCATED(NORMAL_VECTOR_INTERFACE_MESH)) DEALLOCATE(NORMAL_VECTOR_INTERFACE_MESH)

           END IF


         END DO !coupled_mesh_node_idx

         CALL LIST_DETACH_AND_DESTROY(INTER_EDGE_LIST_INTERFACE_MESH, &
           & NUMBER_OF_INTER_EDGES ,INTERFACE_MESH_INTER_EDGE_NODES, ERR,ERROR,*999)

         CALL LIST_DETACH_AND_DESTROY(INTER_EDGE_LIST_COUPLED_MESH, &
           & NUMBER_OF_INTER_EDGES ,COUPLED_MESH_INTER_EDGE_NODES, ERR,ERROR,*999)


         ALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES(NUMBER_OF_INTER_EDGES,2), STAT=ERR)
         IF(ERR /=0) CALL FlagError(" Unable to allocate INTER_EDGES array.",ERR,ERROR,*999)

         COUPLED_DECOMPOSITION%INTER_EDGES = 0

         DO inter_edge_idx = 1,   NUMBER_OF_INTER_EDGES

           coupled_mesh_node_idx =  COUPLED_MESH_INTER_EDGE_NODES(inter_edge_idx)
           interface_node_idx = INTERFACE_MESH_INTER_EDGE_NODES(inter_edge_idx)

           COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2) = coupled_mesh_node_idx
           COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,1) = interface_node_idx


         END DO  !inter_edge_idx

         IF(ALLOCATED(DISTANCE)) DEALLOCATE(DISTANCE)
         IF(ALLOCATED(ADJNCY_RESTRICTED_TO_INTERFACE)) DEALLOCATE(ADJNCY_RESTRICTED_TO_INTERFACE)
         IF(ALLOCATED(BOUNDARY_NODES)) DEALLOCATE(BOUNDARY_NODES)
         IF(ALLOCATED(COUPLED_MESH_INTER_EDGE_NODES)) DEALLOCATE(COUPLED_MESH_INTER_EDGE_NODES)
         IF(ALLOCATED(LIST_OF_NODES)) DEALLOCATE(LIST_OF_NODES)
         IF(ALLOCATED(NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE)) DEALLOCATE(NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE)
         IF(ALLOCATED(INTERFACE_MESH_INTER_EDGE_NODES)) DEALLOCATE(INTERFACE_MESH_INTER_EDGE_NODES)

       ELSE
         CALL FlagError(" Interface mesh is not associated.",ERR,ERROR,*999)
       ENDIF
     ELSE
       CALL FlagError(" COupled mesh is not associated.",ERR,ERROR,*999)
     ENDIF
   ELSE
     CALL FlagError(" Coupled mesh Field is not associated.",ERR,ERROR,*999)
   ENDIF
   EXITS("COUPLED_DECOMPOSITION_GET_INTER_EDGES.")
   RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_GET_INTER_EDGES.",ERR,ERROR)
   RETURN 1

  END SUBROUTINE COUPLED_DECOMPOSITION_GET_INTER_EDGES

!==================================================================================================================
!> Calculate all possible normal angles on node NODE_ID of the coupled mesh graph G_i

  SUBROUTINE COUPLED_MESH_CALCULATE_NORMAL(COUPLED_DECOMPOSITION, NODE_ID, NORMAL_VECTOR, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER ,INTENT(INOUT) :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    INTEGER(INTG), INTENT(OUT)                               :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                        :: ERROR !<The error string.
    INTEGER(INTG)                                            :: NODE_ID !<The global node id to find normal at.
    REAL(RP), ALLOCATABLE                                    :: NORMAL_VECTOR(:,:) !< 2D array containing all possible normal vectors at NODE_ID.
    !Local variables
    INTEGER(INTG)  :: coordinate_idx, counter, MY_COMPUTATIONAL_NODE, node_idx, normal_vector_idx, &
      & normal_vector_idx_1,normal_vector_idx_2, NUMBER_OF_SURROUNDING_NODES, &
      & surrouding_node, surrounding_node_idx, vector_idx
    REAL(RP), ALLOCATABLE                                    :: TANGENT_VECTOR(:,:)
    TYPE(MESH_TYPE), POINTER                                 :: COUPLED_MESH
                                                                ! In the given mesh, for vertex_1, the tangent and normal vectors are calculated as follow:

                                                                !                              ^
                                                                !                              |NORMAL_VECTOR(1,:)
                                                                !                              |
                                                                !                              |
                                                                !                              |        
                                                                !           vertex_1  o---------------->o vertex_2
                                                                !                     | TANGENT_VECTOR(1,:)
                                                                !                     |   
                                                                !                     |
                                                                !                     |
                                                                !                     |
                                                                !NORMAL_VECTOR(2,:)   |TANGENT_VECTOR(2,:)
                                                                !       <-------------|
                                                                !                     |
                                                                !                     V
                                                                !                     o
                                                                !                     vertex_3

                                                                ! Description of the labels:
                                                                ! TANGENT_VECTOR(1,:) is a vector drawn between vertex 1 and vertex 2.
                                                                ! TANGENT_VECTOR(2,:) is a vector drawn between vertex 2 and vertex 3.
                                                                ! TANGENT_VECTOR(3,:) is a vector which is normalized average of  TANGENT_VECTOR(1,:) and TANGENT_VECTOR(2,:).
                                                                ! NORMAL_VECTOR(1,:) is perpendicular to TANGENT_VECTOR(1,:).
                                                                ! NORMAL_VECTOR(2,:) is perpendicular to TANGENT_VECTOR(2,:).
                                                                ! NORMAL_VECTOR(3,:) is perpendicular to TANGENT_VECTOR(3,:).


    ENTERS(" COUPLED_MESH_CALCULATE_NORMAL",ERR,ERROR,*999)
    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      COUPLED_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%MESH
      IF(ASSOCIATED(COUPLED_MESH)) THEN

        MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
        IF(ERR/=0) GOTO 999

        NUMBER_OF_SURROUNDING_NODES =  0
        DO node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODE_ID)%surroundingNodes)

          surrouding_node = COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODE_ID)%surroundingNodes(node_idx)
          IF(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(surrouding_node)%BoundaryNode) &
            & NUMBER_OF_SURROUNDING_NODES = NUMBER_OF_SURROUNDING_NODES +1

        END DO

        ALLOCATE(TANGENT_VECTOR(3,NUMBER_OF_SURROUNDING_NODES), STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate TANGENT_VECTOR array.",ERR,ERROR,*999)
        counter=  0
        DO surrounding_node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODE_ID)%surroundingNodes)

          surrouding_node = COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODE_ID)%surroundingNodes(surrounding_node_idx)
          IF(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(surrouding_node)%BoundaryNode) THEN
            counter= counter+ 1
            TANGENT_VECTOR(:,counter)  = COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES(surrouding_node,:) - &
              &  COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES(NODE_ID,:)

          END IF
        END DO ! surrounding_node_idx


        IF(COUPLED_MESH%NUMBER_OF_DIMENSIONS == 2) THEN


          DO coordinate_idx = 1, 3
            counter= 0
            ! The following do-loop figures out which one of the coordinates in the coupled mesh is all zero.
            ! In other words, it will help algorithm to judge which plane the 2D mesh lie on.             
            DO node_idx = 1, SIZE(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES,1)

              IF(ABS(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES(node_idx, coordinate_idx))<=1E-12) &
                & counter= counter+ 1

            END DO !node_idx
            IF(counter== SIZE(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES,1)) EXIT
          END DO !coordinate_idx

          IF(ALLOCATED(NORMAL_VECTOR)) DEALLOCATE(NORMAL_VECTOR)
          ALLOCATE(NORMAL_VECTOR(3,NUMBER_OF_SURROUNDING_NODES+1), STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate NORMAL_VECTOR array.",ERR,ERROR,*999)

          NORMAL_VECTOR = 0

          DO normal_vector_idx = 1, NUMBER_OF_SURROUNDING_NODES

            IF(coordinate_idx == 1) THEN

              NORMAL_VECTOR(1,normal_vector_idx) = 0
              NORMAL_VECTOR(2,normal_vector_idx) = TANGENT_VECTOR(3, normal_vector_idx)
              NORMAL_VECTOR(3,normal_vector_idx) = -TANGENT_VECTOR(2, normal_vector_idx)

            ELSE IF(coordinate_idx == 2) THEN

              NORMAL_VECTOR(1,normal_vector_idx) = TANGENT_VECTOR(3, normal_vector_idx)
              NORMAL_VECTOR(2,normal_vector_idx) = 0
              NORMAL_VECTOR(3,normal_vector_idx) = -TANGENT_VECTOR(1, normal_vector_idx)


            ELSE

              NORMAL_VECTOR(1,normal_vector_idx) = TANGENT_VECTOR(2, normal_vector_idx)
              NORMAL_VECTOR(2,normal_vector_idx) = -TANGENT_VECTOR(1, normal_vector_idx)
              NORMAL_VECTOR(3,normal_vector_idx) = 0

            END IF

            NORMAL_VECTOR(:,normal_vector_idx) = NORMAL_VECTOR(:,normal_vector_idx)/ &
              & NORM2(NORMAL_VECTOR(:,normal_vector_idx))

          END DO   ! normal_vector_idx

          DO normal_vector_idx = 1, NUMBER_OF_SURROUNDING_NODES

            NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) = NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) + &
              & NORMAL_VECTOR(:,normal_vector_idx)

          END DO

            NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) = NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1)/ &   ! Averaging out the vector
              & NUMBER_OF_SURROUNDING_NODES

            NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) = NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1)/ &
              & NORM2(NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1))

        ELSE IF(COUPLED_MESH%NUMBER_OF_DIMENSIONS == 3) THEN

          IF(ALLOCATED(NORMAL_VECTOR)) DEALLOCATE(NORMAL_VECTOR)
          ALLOCATE( &
            & NORMAL_VECTOR(3,INT(0.5*(NUMBER_OF_SURROUNDING_NODES-1)*(NUMBER_OF_SURROUNDING_NODES)+1)), STAT=ERR) ! (n-1)n/2 + 1
          IF(ERR/=0) CALL FlagError("Could not allocate NORMAL_VECTOR array.",ERR,ERROR,*999)

          NORMAL_VECTOR = 0
          vector_idx = 0

          DO normal_vector_idx_1 = 1, NUMBER_OF_SURROUNDING_NODES-1
            DO normal_vector_idx_2 = normal_vector_idx_1+1, NUMBER_OF_SURROUNDING_NODES
             vector_idx = vector_idx + 1
             !  cross = axb
             !  cross(1) = a(2) * b(3) - a(3) * b(2)
             !  cross(2) = a(3) * b(1) - a(1) * b(3)
             !  cross(3) = a(1) * b(2) - a(2) * b(1)

              NORMAL_VECTOR(1,vector_idx)=TANGENT_VECTOR(2,normal_vector_idx_1)*TANGENT_VECTOR(3,normal_vector_idx_2)- &
                & TANGENT_VECTOR(3,normal_vector_idx_1)*TANGENT_VECTOR(2,normal_vector_idx_2)

              NORMAL_VECTOR(2,vector_idx)=TANGENT_VECTOR(3,normal_vector_idx_1)*TANGENT_VECTOR(1,normal_vector_idx_2)- &
                & TANGENT_VECTOR(1,normal_vector_idx_1)*TANGENT_VECTOR(3,normal_vector_idx_2)

              NORMAL_VECTOR(3,vector_idx)=TANGENT_VECTOR(1,normal_vector_idx_1)*TANGENT_VECTOR(2,normal_vector_idx_2)- &
                & TANGENT_VECTOR(2,normal_vector_idx_1)*TANGENT_VECTOR(1,normal_vector_idx_2)

              NORMAL_VECTOR(:,vector_idx) = NORMAL_VECTOR(:,vector_idx)/NORM2(NORMAL_VECTOR(:,vector_idx))
            END DO   ! normal_vector_idx_1
          END DO     ! normal_vector_idx_2

          DO normal_vector_idx = 1, SIZE(NORMAL_VECTOR,2)-1

            NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) = NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) + &
              & NORMAL_VECTOR(:,normal_vector_idx)

          END DO

          NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) = NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2))/ &   ! Averaging out the vector
            & SIZE(NORMAL_VECTOR,2)-1

          NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) = NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2))/ &
            & NORM2(NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)))

       END IF !(COUPLED_MESH%NUMBER_OF_DIMENSIONS == 2)

       IF(ALLOCATED(TANGENT_VECTOR)) DEALLOCATE(TANGENT_VECTOR)

     ELSE
       CALL FlagError(" Coupled mesh is not associated.",ERR,ERROR,*999)
     ENDIF
   ELSE
     CALL FlagError(" Coupled decomposition is not associated.",ERR,ERROR,*999)
   ENDIF
   EXITS(" COUPLED_MESH_CALCULATE_NORMAL.")
   RETURN
999 ERRORSEXITS(" COUPLED_MESH_CALCULATE_NORMAL.",ERR,ERROR)
   RETURN 1
  END SUBROUTINE  COUPLED_MESH_CALCULATE_NORMAL

!=================================================================================================================
!> Calculate all possible normal angles on node NODE_ID of the interface mesh graph G_I

  SUBROUTINE INTERFACE_MESH_CALCULATE_NORMAL(COUPLED_DECOMPOSITION, NODE_ID, NORMAL_VECTOR, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER ,INTENT(INOUT) :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    INTEGER(INTG), INTENT(OUT)                               :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                        :: ERROR !<The error string.
    INTEGER(INTG)                                            :: NODE_ID !<The global node id to find normal at.
    REAL(RP), ALLOCATABLE                                    :: NORMAL_VECTOR(:,:) !< 2D array containing all possible normal vectors at NODE_ID.

                                                                ! In the given mesh, for vertex_1, the tangent and normal vectors are calculated as follow:

                                                                !                              ^
                                                                !                              |NORMAL_VECTOR(1,:)
                                                                !                              |
                                                                !                              |
                                                                !                              |        
                                                                !           vertex_1  o---------------->o vertex_2
                                                                !                     | TANGENT_VECTOR(1,:)
                                                                !                     |   
                                                                !                     |
                                                                !                     |
                                                                !                     |
                                                                !NORMAL_VECTOR(2,:)   |TANGENT_VECTOR(2,:)
                                                                !       <-------------|
                                                                !                     |
                                                                !                     V
                                                                !                     o
                                                                !                     vertex_3

                                                                ! Description of the labels:
                                                                ! TANGENT_VECTOR(1,:) is a vector drawn between vertex 1 and vertex 2
                                                                ! TANGENT_VECTOR(2,:) is a vector drawn between vertex 2 and vertex 3
                                                                ! TANGENT_VECTOR(3,:) is a vector which is normalized average of  TANGENT_VECTOR(1,:) and TANGENT_VECTOR(2,:)
                                                                ! NORMAL_VECTOR(1,:) is perpendicular to TANGENT_VECTOR(1,:)
                                                                ! NORMAL_VECTOR(2,:) is perpendicular to TANGENT_VECTOR(2,:)
                                                                ! NORMAL_VECTOR(3,:) is perpendicular to TANGENT_VECTOR(3,:)
                                                                                                                                    
    !Local variables
    INTEGER(INTG)  :: coordinate_idx, counter, MY_COMPUTATIONAL_NODE, node_idx, normal_vector_idx, &
      & normal_vector_idx_1,normal_vector_idx_2, NUMBER_OF_SURROUNDING_NODES, &
      & surrouding_node, surrounding_node_idx, vector_idx
    REAL(RP), ALLOCATABLE                                    :: TANGENT_VECTOR(:,:)!                          TANGENT_VECTOR(1,:)         TANGENT_VECTOR(2,:)
                                                                                   ! For instance 1D mesh: o-<--------------------o------------------------->-o
    TYPE(MESH_TYPE), POINTER                                 :: INTERFACE_MESH, COUPLED_MESH              !2                      1                           3 


    ENTERS(" INTERFACE_MESH_CALCULATE_NORMAL",ERR,ERROR,*999)
    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH
      IF(ASSOCIATED(INTERFACE_MESH)) THEN
        COUPLED_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(1)%PTR%DECOMPOSITION%MESH
        MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
        IF(ERR/=0) GOTO 999


        CALL GET_SURROUNDING_NODES(INTERFACE_MESH, ERR, error, *999)
        NUMBER_OF_SURROUNDING_NODES =  SIZE(INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODE_ID)%surroundingNodes)

        ALLOCATE(TANGENT_VECTOR(3,NUMBER_OF_SURROUNDING_NODES), STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate TANGENT_VECTOR array.",ERR,ERROR,*999)

        DO surrounding_node_idx = 1, NUMBER_OF_SURROUNDING_NODES

          surrouding_node = INTERFACE_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODE_ID)%surroundingNodes(surrounding_node_idx)

          TANGENT_VECTOR(:,surrounding_node_idx)=COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(surrouding_node,:)-&
            &  COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(NODE_ID,:)

        END DO ! surrounding_node_idx


        IF(COUPLED_MESH%NUMBER_OF_DIMENSIONS == 2) THEN

            ! The following do-loop figures out which one of the coordinates in the coupled mesh is all zero.
            ! In other words, it will help algorithm to judge which plane the 2D mesh lie on.
          DO coordinate_idx = 1, 3
            counter= 0
            DO node_idx = 1, SIZE(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES,1)

              IF(ABS(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES(node_idx, coordinate_idx))<=1E-12) &
                & counter= counter+ 1

            END DO !node_idx
            IF(counter== SIZE(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES,1)) EXIT
          END DO !coordinate_idx

          IF(ALLOCATED(NORMAL_VECTOR)) DEALLOCATE(NORMAL_VECTOR)
          ALLOCATE(NORMAL_VECTOR(3,NUMBER_OF_SURROUNDING_NODES+1), STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate NORMAL_VECTOR array.",ERR,ERROR,*999)

          NORMAL_VECTOR = 0

          DO normal_vector_idx = 1, NUMBER_OF_SURROUNDING_NODES

            IF(coordinate_idx == 1) THEN

              NORMAL_VECTOR(1,normal_vector_idx) = 0
              NORMAL_VECTOR(2,normal_vector_idx) = TANGENT_VECTOR(3, normal_vector_idx)
              NORMAL_VECTOR(3,normal_vector_idx) = -TANGENT_VECTOR(2, normal_vector_idx)

            ELSE IF(coordinate_idx == 2) THEN

              NORMAL_VECTOR(1,normal_vector_idx) = TANGENT_VECTOR(3, normal_vector_idx)
              NORMAL_VECTOR(2,normal_vector_idx) = 0
              NORMAL_VECTOR(3,normal_vector_idx) = -TANGENT_VECTOR(1, normal_vector_idx)


            ELSE

              NORMAL_VECTOR(1,normal_vector_idx) = TANGENT_VECTOR(2, normal_vector_idx)
              NORMAL_VECTOR(2,normal_vector_idx) = -TANGENT_VECTOR(1, normal_vector_idx)
              NORMAL_VECTOR(3,normal_vector_idx) = 0

            END IF

            NORMAL_VECTOR(:,normal_vector_idx) = NORMAL_VECTOR(:,normal_vector_idx)/ &
              & NORM2(NORMAL_VECTOR(:,normal_vector_idx))

          END DO   ! normal_vector_idx

          DO normal_vector_idx = 1, NUMBER_OF_SURROUNDING_NODES

            NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) = NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) + &
              & NORMAL_VECTOR(:,normal_vector_idx)

          END DO

          NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) = NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1)/ &   ! Averaging out the vector
            & NUMBER_OF_SURROUNDING_NODES

          NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1) = NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1)/ &
            & NORM2(NORMAL_VECTOR(:,NUMBER_OF_SURROUNDING_NODES+1))

        ELSE IF(COUPLED_MESH%NUMBER_OF_DIMENSIONS == 3) THEN

          IF(ALLOCATED(NORMAL_VECTOR)) DEALLOCATE(NORMAL_VECTOR)
          ALLOCATE(NORMAL_VECTOR(3,INT(0.5*(NUMBER_OF_SURROUNDING_NODES-1)*(NUMBER_OF_SURROUNDING_NODES)+1)), &
            & STAT=ERR) ! (n-1)n/2 + 1
          IF(ERR/=0) CALL FlagError("Could not allocate NORMAL_VECTOR array.",ERR,ERROR,*999)

          NORMAL_VECTOR = 0
          vector_idx = 0

          DO normal_vector_idx_1 = 1, NUMBER_OF_SURROUNDING_NODES-1
            DO normal_vector_idx_2 = normal_vector_idx_1+1, NUMBER_OF_SURROUNDING_NODES
             vector_idx = vector_idx + 1
             !  cross = axb
             !  cross(1) = a(2) * b(3) - a(3) * b(2)
             !  cross(2) = a(3) * b(1) - a(1) * b(3)
             !  cross(3) = a(1) * b(2) - a(2) * b(1)

              NORMAL_VECTOR(1,vector_idx)=TANGENT_VECTOR(2,normal_vector_idx_1)*TANGENT_VECTOR(3,normal_vector_idx_2)- &
                & TANGENT_VECTOR(3,normal_vector_idx_1)*TANGENT_VECTOR(2,normal_vector_idx_2)

              NORMAL_VECTOR(2,vector_idx)=TANGENT_VECTOR(3,normal_vector_idx_1)*TANGENT_VECTOR(1,normal_vector_idx_2)- &
                & TANGENT_VECTOR(1,normal_vector_idx_1)*TANGENT_VECTOR(3,normal_vector_idx_2)

              NORMAL_VECTOR(3,vector_idx)=TANGENT_VECTOR(1,normal_vector_idx_1)*TANGENT_VECTOR(2,normal_vector_idx_2)- &
                & TANGENT_VECTOR(2,normal_vector_idx_1)*TANGENT_VECTOR(1,normal_vector_idx_2)

              NORMAL_VECTOR(:,vector_idx) = NORMAL_VECTOR(:,vector_idx)/NORM2(NORMAL_VECTOR(:,vector_idx))
            END DO   ! normal_vector_idx_1
          END DO     ! normal_vector_idx_2

          DO normal_vector_idx = 1, SIZE(NORMAL_VECTOR,2)-1

            NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) = NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) + &
              & NORMAL_VECTOR(:,normal_vector_idx)

          END DO

          NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) = NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2))/ &   ! Averaging out the vector
            & SIZE(NORMAL_VECTOR,2)-1

          NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)) = NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2))/ &
            & NORM2(NORMAL_VECTOR(:,SIZE(NORMAL_VECTOR,2)))

       END IF !(INTERFACE_MESH%NUMBER_OF_DIMENSIONS == 2)

     ELSE
       CALL FlagError(" Coupled mesh is not associated.",ERR,ERROR,*999)
     ENDIF
   ELSE
     CALL FlagError(" Coupled decomposition is not associated.",ERR,ERROR,*999)
   ENDIF
   EXITS(" INTERFACE_MESH_CALCULATE_NORMAL.")
   RETURN
999 ERRORSEXITS(" INTERFACE_MESH_CALCULATE_NORMAL.",ERR,ERROR)
   RETURN 1
  END SUBROUTINE  INTERFACE_MESH_CALCULATE_NORMAL

!=================================================================================================================
!> This subroutine implements the coupling-aware partitioning of multidomain problems arising from
!> interface-coupled problems using vertex-based partitioning as described in Waleed Mirza's M.Sc. thesis (referred to as Scheme B)

  SUBROUTINE COUPLED_DECOMPOSITION_FIXED_VERTEX_PARTITIONING(COUPLED_DECOMPOSITION, ERR, ERROR, *)

  !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER ,INTENT(INOUT) :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    INTEGER(INTG), INTENT(OUT)                               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)                        :: ERROR !<The error string

    !Local Variables
    REAL(RP)                                   :: INTERFACE_MESH_MAXIMUM_EDGE_LENGTH
    TYPE(MESH_TYPE), POINTER                   :: COUPLED_MESH, INTERFACE_MESH, NEW_COUPLED_MESH
    TYPE(FIELD_TYPE), POINTER                  :: FIELD_COUPLED_MESH, FIELD_INTERFACE_MESH
    REAL(RP), ALLOCATABLE                      :: NODE_WEIGHTS_TO_BE_PROJECTED(:)
    TYPE(DECOMPOSITION_TYPE), POINTER          :: COUPLED_MESH_DECOMPOSITION, INTERFACE_MESH_DECOMPOSITION
    INTEGER(INTG)                              :: COORDINATE_SYSTEM_USER_NUMBER, element_idx, MESH_ELEMENTS, &
      & NUM, NUMBER_OF_COMPUTATIONAL_NODES, PARMETIS_OPTIONS(3), proc_idx_send,  &
      & MY_COMPUTATIONAL_NODE, node_idx, STATUS(MPI_STATUS_SIZE), REGION_SYSTEM_USER_NUMBER, &
      & vtx_dist_idx, proc_idx_receive, proc_idx
    INTEGER(INTG), ALLOCATABLE                 :: FLIP_PARTITION(:,:), NEW_TO_OLD_INDEX_MAPPING(:,:), &
      & NODES_TO_IMPOSE_WEIGHTS_ON(:), TEMP_ARRAY(:)
    INTEGER(INTG)                              :: xadj_counter, adjncy_size, counter, new_node_idx
    TYPE(REGION_TYPE), POINTER                 :: COUPLED_MESH_REGION
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER      :: COUPLED_MESH_COORDINATE_SYSTEM
    TYPE(NODES_TYPE), POINTER                  :: COUPLED_MESH_NODE
    TYPE(meshElementsType), POINTER            :: COUPLED_MESH_ELEMENTS
    LOGICAL                                    :: FLAG


    ENTERS("COUPLED_DECOMPOSITION_FIXED_VERTEX_PARTITIONING",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
        COUPLED_MESH=>&
          COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%MESH

      IF(ASSOCIATED(COUPLED_MESH)) THEN
        FIELD_COUPLED_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR

        IF(ASSOCIATED(FIELD_COUPLED_MESH)) THEN
          INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH

          IF(ASSOCIATED(INTERFACE_MESH)) THEN
            FIELD_INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR

            IF(ASSOCIATED(FIELD_INTERFACE_MESH)) THEN
              COUPLED_MESH_DECOMPOSITION=> &
                & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION

              IF(ASSOCIATED(COUPLED_MESH_DECOMPOSITION)) THEN
                INTERFACE_MESH_DECOMPOSITION=> COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION
                IF(ASSOCIATED(INTERFACE_MESH_DECOMPOSITION)) THEN

                  MY_COMPUTATIONAL_NODE = COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  NUMBER_OF_COMPUTATIONAL_NODES = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999

                  CALL GET_SURROUNDING_NODES(COUPLED_MESH, ERR, error, *999)

                  CALL GATHER_MESH_COORDINATES(FIELD_COUPLED_MESH, COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES,&
                    & ERR, Error, *999) ! Gather the geometric field of coupled mesh graph denoted by mesh_idx

                  CALL GET_MAXIMUM_EDGE_LENGTH(COUPLED_DECOMPOSITION, &
                    & INTERFACE_MESH_MAXIMUM_EDGE_LENGTH, ERR, ERROR, *999) ! Get the maximum edge length of the interface graph.

                  IF(ALLOCATED(COUPLED_DECOMPOSITION%INTER_EDGES)) DEALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES) !WOuld have meaning for mesh_idx=2

                  CALL COUPLED_DECOMPOSITION_GET_INTER_EDGES(COUPLED_DECOMPOSITION, &
                    & INTERFACE_MESH_MAXIMUM_EDGE_LENGTH, ERR, ERROR, *999) !Build interedges between the coupled mesh graph and the interface graph.

                  CALL COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING(COUPLED_DECOMPOSITION,ERR,ERROR,*999) !Gather the vertices of the coupled mesh graph where fixed partitioning is to be imposed.
                  CALL COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING( &
                    & COUPLED_DECOMPOSITION,COUPLED_DECOMPOSITION%OLD_TO_NEW_VERTEX_MAPPING,ERR,ERROR,*999) !Build a 2D array OLD_TO_NEW_VERTEX_MAPPING(: , :) that stores a one-to-one relationship between old vertices and the new vertices.

                 ! From here onward initialize a mesh object which is an identical copy of coupled_mesh object.

                 ! The following extracts the coordindate system user number that has been been used before in create a coordinate system.
                  NULLIFY(COUPLED_MESH_COORDINATE_SYSTEM)
                  FLAG = .TRUE.
                  DO WHILE (FLAG)
                    COORDINATE_SYSTEM_USER_NUMBER = IRAND()
                    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(COORDINATE_SYSTEM_USER_NUMBER, &
                      & COUPLED_MESH_COORDINATE_SYSTEM,ERR,ERROR,*999)
                    IF(ASSOCIATED(COUPLED_MESH_COORDINATE_SYSTEM)) THEN
                      FLAG = .TRUE.
                      NULLIFY(COUPLED_MESH_COORDINATE_SYSTEM)
                    ELSE
                      FLAG = .FALSE.
                    END IF
                  END DO
                 ! Initialite the coordinate system.
                  CALL COORDINATE_SYSTEM_CREATE_START(&
                    & COORDINATE_SYSTEM_USER_NUMBER,COUPLED_MESH_COORDINATE_SYSTEM,ERR,ERROR,*999)
                  CALL COORDINATE_SYSTEM_DIMENSION_SET(&
                    & COUPLED_MESH_COORDINATE_SYSTEM,COUPLED_MESH%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  CALL COORDINATE_SYSTEM_CREATE_FINISH(COUPLED_MESH_COORDINATE_SYSTEM,ERR,ERROR,*999)


                 ! The following extracts the region user number that has been been used before in create a region.

                  FLAG = .TRUE.
                  DO WHILE (FLAG)
                    REGION_SYSTEM_USER_NUMBER = IRAND()
                    NULLIFY(COUPLED_MESH_REGION)
                    CALL REGION_USER_NUMBER_FIND(REGION_SYSTEM_USER_NUMBER,COUPLED_MESH_REGION,ERR,ERROR,*999)
                    IF(ASSOCIATED(COUPLED_MESH_REGION)) THEN
                      FLAG = .TRUE.
                    ELSE
                      FLAG = .FALSE.
                    END IF
                  END DO
                  ! Initialize the mesh region.
                  CALL REGION_INITIALISE(COUPLED_MESH_REGION,ERR,ERROR,*999)
                  NULLIFY(COUPLED_MESH_REGION)
                  CALL REGION_CREATE_START(REGION_SYSTEM_USER_NUMBER, &
                    & COUPLED_MESH%REGION%PARENT_REGION,COUPLED_MESH_REGION,ERR,ERROR,*999)
                  CALL REGION_COORDINATE_SYSTEM_SET(COUPLED_MESH_REGION,COUPLED_MESH_COORDINATE_SYSTEM,ERR,ERROR,*999)
                  CALL REGION_CREATE_FINISH(COUPLED_MESH_REGION,ERR,ERROR,*999)
                  ! Initialize the mesh coupled mesh object from here onward.
                  CALL MESH_NUMBER_OF_ELEMENTS_GET(COUPLED_MESH,MESH_ELEMENTS,ERR,ERROR,*999)
                  NULLIFY(COUPLED_MESH_NODE)
                  CALL NODES_CREATE_START(COUPLED_MESH_REGION,MAXVAL(COUPLED_DECOMPOSITION%OLD_TO_NEW_VERTEX_MAPPING(:,1)), &
                    & COUPLED_MESH_NODE,ERR,ERROR,*999)
                  CALL NODES_CREATE_FINISH(COUPLED_MESH_NODE,ERR,ERROR,*999)
                  NULLIFY(NEW_COUPLED_MESH)
                  CALL MESH_INITIALISE(NEW_COUPLED_MESH,ERR,ERROR,*999)
                  NULLIFY(NEW_COUPLED_MESH)
                  CALL MESH_CREATE_START(COUPLED_DECOMPOSITION%mesh_idx,COUPLED_MESH_REGION, &
                    & COUPLED_MESH%NUMBER_OF_DIMENSIONS,NEW_COUPLED_MESH,ERR,ERROR,*999)
                  CALL MESH_NUMBER_OF_ELEMENTS_SET(NEW_COUPLED_MESH,MESH_ELEMENTS,ERR,ERROR,*999)
                  CALL MESH_NUMBER_OF_COMPONENTS_SET(NEW_COUPLED_MESH,COUPLED_MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                  NULLIFY(COUPLED_MESH_ELEMENTS)
                  CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(NEW_COUPLED_MESH,1, &
                    & COUPLED_MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(1)%BASIS, COUPLED_MESH_ELEMENTS,ERR,ERROR,*999)
                  DO element_idx = 1, MESH_ELEMENTS
                    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(element_idx,COUPLED_MESH_ELEMENTS, &
                      & COUPLED_MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(element_idx)%MESH_ELEMENT_NODES,ERR,ERROR,*999)
                  END DO !element_idx
                  CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(COUPLED_MESH_ELEMENTS,ERR,ERROR,*999)
                  CALL MESH_CREATE_FINISH(NEW_COUPLED_MESH,ERR,ERROR,*999)

                  ! Calculate surrounding nodes.
                  CALL GET_SURROUNDING_NODES(NEW_COUPLED_MESH, ERR, error, *999)
                  ! Alter the graph adjacencies for mesh fixed partitioning.
                  CALL COUPLED_DECOMPOSITION_GET_NEW_GRAPH(COUPLED_DECOMPOSITION,NEW_COUPLED_MESH, &
                    & COUPLED_DECOMPOSITION%OLD_TO_NEW_VERTEX_MAPPING,ERR,ERROR,*999)

                  NEW_COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNOdes = SIZE(NEW_COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)

                  ! Calculate the decomposition parameters.
                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN)
                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN &
                    & (0:NEW_COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNOdes-1), STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate NODE_DOMAIN array.",ERR,ERROR,*999)


                  COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN = 0

                  COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS = NUMBER_OF_COMPUTATIONAL_NODES

                  COUPLED_MESH_DECOMPOSITION%WEIGHT_FLAG = 2_INTG ! edge and vertex eights activated
                  COUPLED_MESH_DECOMPOSITION%NUM_FLAG    = 0_INTG ! node numbering starts from 0

                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%VTX_DIST)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%VTX_DIST)
                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%VTX_DIST(0:COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate VTX_DIST array.",ERR,ERROR,*999)

                  ! Calculate VTX_DIST
                  COUPLED_MESH_DECOMPOSITION%VTX_DIST(0)  = 0
                  COUPLED_MESH_DECOMPOSITION%VTX_DIST(COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS) =  &
                    & NEW_COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes
                  ! Calculate COUPLED_MESH_DECOMPOSITION%VTX_DIST
                  DO vtx_dist_idx = 1, COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS-1

                   COUPLED_MESH_DECOMPOSITION%VTX_DIST(vtx_dist_idx) =  COUPLED_MESH_DECOMPOSITION%VTX_DIST(vtx_dist_idx-1) + &
                    & INT(NEW_COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes/COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS, INTG)

                  END DO !vtx_dist_idx

                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%XADJ)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%XADJ)
                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%XADJ(0:(COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1)-&
                    & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE))), STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate XADJ array.",ERR,ERROR,*999)

                  COUPLED_MESH_DECOMPOSITION%XADJ(0) = 0
                  xadj_counter= 0
                  adjncy_size = 0

                ! Calculate COUPLED_MESH_DECOMPOSITION%XADJ
                  DO vtx_dist_idx = COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE), &
                    & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1)-1

                    xadj_counter= xadj_counter+ 1

                    adjncy_size = adjncy_size + &
                      & SIZE(NEW_COUPLED_MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes,1)
                    COUPLED_MESH_DECOMPOSITION%XADJ(xadj_counter) = adjncy_size

                  END DO !vtx_dist_idx

                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%ADJNCY)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%ADJNCY)
                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%ADJNCY(0:adjncy_size-1), STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ADJNCY array.",ERR,ERROR,*999)

                  counter= 0
                  DO vtx_dist_idx = COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE), &
                    & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1)-1


                    COUPLED_MESH_DECOMPOSITION%ADJNCY( &
                      & counter: counter+SIZE(NEW_COUPLED_MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes)-1)= &
                      & NEW_COUPLED_MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes

                    counter= counter+ &
                      & SIZE(NEW_COUPLED_MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes,1)

                  END DO !vtx_dist_idx.
                  COUPLED_MESH_DECOMPOSITION%ADJNCY = COUPLED_MESH_DECOMPOSITION%ADJNCY - 1 ! to make vertex numbering start from 0.

                  COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS = 2

                  ! Define array ADJWT with edge weights.
                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%ADJWT)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%ADJWT)
                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%ADJWT(SIZE(COUPLED_MESH_DECOMPOSITION%ADJNCY,1)), STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ADJWT array.",ERR,ERROR,*999)
                  COUPLED_MESH_DECOMPOSITION%ADJWT = 1  ! Set to default value

                  ALLOCATE(NODES_TO_IMPOSE_WEIGHTS_ON(NUMBER_OF_COMPUTATIONAL_NODES), STAT=ERR) ! The array stores the vertices of the new coupled mesh graph G_{i,merged} where vertex weights are supposed to be imposed.
                  IF(ERR/=0) CALL FlagError("Could not allocate NODES_TO_IMPOSE_WEIGHTS_ON array.",ERR,ERROR,*999)

                  DO proc_idx = 1, NUMBER_OF_COMPUTATIONAL_NODES
                    node_idx = COUPLED_DECOMPOSITION%INTER_EDGES(1,proc_idx)
                    NODES_TO_IMPOSE_WEIGHTS_ON(proc_idx) = COUPLED_DECOMPOSITION%OLD_TO_NEW_VERTEX_MAPPING(node_idx, 2)
                  END DO

        ! initialize  NODE_WEIGHT_SET and impose node weights
                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET)
                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET( &
                    & 1:COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS*( &
                    & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1) &
                    & - COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE))), STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate NODE_WEIGHT_SET array.",ERR,ERROR,*999)
                  COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET   = 1

                  DO node_idx = 1, SIZE(NODES_TO_IMPOSE_WEIGHTS_ON,1)

                    IF(NODES_TO_IMPOSE_WEIGHTS_ON(node_idx) <= &
                      & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1) .AND. &
                      & NODES_TO_IMPOSE_WEIGHTS_ON(node_idx) > &
                      & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE)) THEN

                      COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET&
                        & ((NODES_TO_IMPOSE_WEIGHTS_ON(node_idx)-COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE))*&
                        & COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS) =  &
                        & SIZE(COUPLED_MESH%TOPOLOGY(1)%ptr%Nodes%Nodes,1)

                      COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET&
                        & ((NODES_TO_IMPOSE_WEIGHTS_ON(node_idx)-COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE))*&
                        & COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS-1) = &
                        & SIZE(PACK(COUPLED_DECOMPOSITION%INTER_EDGES(:,node_idx), &
                        & COUPLED_DECOMPOSITION%INTER_EDGES(:,node_idx) /= 0))

                    END IF
                  END DO !node_idx


                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%TPWGT)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%TPWGT)
                  IF(ALLOCATED(COUPLED_MESH_DECOMPOSITION%UBVEC)) DEALLOCATE(COUPLED_MESH_DECOMPOSITION%UBVEC)

                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%TPWGT &
                    & (COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS*COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
                  IF(ERR/=0)  CALL FlagError("Could not allocate TPWGT array",ERR,ERROR,*999)

                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%UBVEC(COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS), STAT=ERR)
                  IF(ERR/=0)  CALL FlagError("Could not allocate UBVEC array",ERR,ERROR,*999)

                  COUPLED_MESH_DECOMPOSITION%TPWGT = REAL(1./COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS,RP)
                  COUPLED_MESH_DECOMPOSITION%UBVEC = 1.00001_RP

                  PARMETIS_OPTIONS(1) = 1
                  PARMETIS_OPTIONS(2) = 7
                  PARMETIS_OPTIONS(3) = 99999

                  ! Calculate node/vertex domains.
                  CALL ParMETIS_V3_PartKway(COUPLED_MESH_DECOMPOSITION%VTX_DIST,COUPLED_MESH_DECOMPOSITION%XADJ, &
                    & COUPLED_MESH_DECOMPOSITION%ADJNCY,COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET, &
                    & COUPLED_MESH_DECOMPOSITION%ADJWT, COUPLED_MESH_DECOMPOSITION%WEIGHT_FLAG, &
                    & COUPLED_MESH_DECOMPOSITION%NUM_FLAG, 2_INTG, &
                    & COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS, COUPLED_MESH_DECOMPOSITION%TPWGT, &
                    & COUPLED_MESH_DECOMPOSITION%UBVEC, PARMETIS_OPTIONS, &
                    & COUPLED_MESH_DECOMPOSITION%NUMBER_OF_EDGES_CUT, &
                    & COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN( &
                    & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE): &
                    & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1)-1), &
                    & COMPUTATIONAL_ENVIRONMENT%MPI_COMM)

                  ! Store all local arrays in one array.
                  DO proc_idx_receive = 0, COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS-1

                    IF(MY_COMPUTATIONAL_NODE /= proc_idx_receive) THEN

                      CALL MPI_SEND(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE):&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1)-1),&
                        & SIZE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE):&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1)-1)),&
                        & MPI_INT, proc_idx_receive, MY_COMPUTATIONAL_NODE ,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, Err)


                    END IF
                  END DO !!ROC_idx_RECEIVE

                  DO proc_idx_send = 0, COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS-1


                    IF(proc_idx_send /= MY_COMPUTATIONAL_NODE) THEN

                      CALL MPI_RECV( &
                        & COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(proc_idx_send):&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(proc_idx_send+1)-1),&
                        & SIZE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(proc_idx_send):&
                        & COUPLED_MESH_DECOMPOSITION%VTX_DIST(proc_idx_send+1)-1)),&
                        & MPI_INT, proc_idx_send,proc_idx_send, COMPUTATIONAL_ENVIRONMENT%MPI_COMM,STATUS, Err)


                     END IF !proc_idx_send

                  END DO !procidx_send = 0, number_parts-1


                  ! The following flip the sub-domain ids if necessary.
                  ALLOCATE(FLIP_PARTITION(COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS,2), STAT=ERR)
                  IF(ERR/=0)  CALL FlagError("Could not allocate FLIP_PARTITION array",ERR,ERROR,*999)

                  DO proc_idx = 1, NUMBER_OF_COMPUTATIONAL_NODES

                    node_idx = COUPLED_DECOMPOSITION%OLD_TO_NEW_VERTEX_MAPPING(COUPLED_DECOMPOSITION%INTER_EDGES(1,proc_idx),2)
                    FLIP_PARTITION(proc_idx,1) = proc_idx-1
                    FLIP_PARTITION(proc_idx,2) = COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_idx-1)

                  END DO

                  DO node_idx = 1, SIZE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN,1)

                    NUM = INT(SQRT(REAL(DOT_PRODUCT(FLIP_PARTITION &
                      & (MINLOC(ABS(FLIP_PARTITION(:,2)-COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_idx-1))),1), &
                      & FLIP_PARTITION(MINLOC(ABS(FLIP_PARTITION(:,2)-COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_idx-1))) &
                      & ,1)))),INTG)

                    COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_idx-1) = NUM

                  END DO

                  ALLOCATE(TEMP_ARRAY(0:SIZE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN)-1), STAT=ERR)
                  IF(ERR/=0)  CALL FlagError("Could not allocate TEMP_ARRAY array.",ERR,ERROR,*999)

                  TEMP_ARRAY = COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN

                  DEALLOCATE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN)
                  ALLOCATE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(0:COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNodes-1), STAT=ERR)
                  IF(ERR/=0)  CALL FlagError("Could not allocate NODE_DOMAIN array.",ERR,ERROR,*999)

                  !Project vetex domains from now coupled mesh graph G_{i,merged} to the old coupled mesh graph G_i.
                  DO node_idx= 1, COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNodes

                    new_node_idx = COUPLED_DECOMPOSITION%OLD_TO_NEW_VERTEX_MAPPING(node_idx,2)

                    COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_idx-1) = TEMP_ARRAY(new_node_idx-1)

                  END DO !node_idx


                  COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION=>COUPLED_MESH_DECOMPOSITION
                  COUPLED_DECOMPOSITION%mesh_idx = COUPLED_DECOMPOSITION%mesh_idx + 1

                  NULLIFY(COUPLED_MESH_COORDINATE_SYSTEM)
                  NULLIFY(COUPLED_MESH_REGION)
                  NULLIFY(COUPLED_MESH_NODE)
                  NULLIFY(COUPLED_MESH_ELEMENTS)
                  IF(ALLOCATED(COUPLED_DECOMPOSITION%INTER_EDGES)) DEALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES)
                  IF(ALLOCATED(FLIP_PARTITION)) DEALLOCATE(FLIP_PARTITION)
                  IF(ALLOCATED(NEW_TO_OLD_INDEX_MAPPING)) DEALLOCATE(NEW_TO_OLD_INDEX_MAPPING)
                  IF(ALLOCATED(TEMP_ARRAY)) DEALLOCATE(TEMP_ARRAY)
                  IF(ALLOCATED(NODE_WEIGHTS_TO_BE_PROJECTED)) DEALLOCATE( NODE_WEIGHTS_TO_BE_PROJECTED)

                ELSE
                  CALL FlagError("Interface mesh decomposition is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface mesh decomposition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Coupled mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Field coupled mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Field Interface mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
       CALL FlagError("Coupled mesh decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
     CALL FlagError("Coupled decomposition  is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("COUPLED_DECOMPOSITION_FIXED_VERTEX_PARTITIONING.")

    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_FIXED_VERTEX_PARTITIO",err,error)
    RETURN 1

  END SUBROUTINE COUPLED_DECOMPOSITION_FIXED_VERTEX_PARTITIONING


!========================================================================================
 !> The following subroutine trivially decomposes the interface graph G_I.

  SUBROUTINE DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMPOSITION_SET(INTERFACE_MESH, ERR, ERROR, *)

  !Argument variables
    TYPE(MESH_TYPE), POINTER, INTENT(IN)      :: INTERFACE_MESH !<A pointer to the interface mesh object.
    INTEGER(INTG), INTENT(OUT)                :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)         :: ERROR !<The error string

    ! Local argument
    TYPE(DECOMPOSITION_TYPE), POINTER         :: INTERFACE_DECOMPOSITION
    INTEGER(INTG)                             :: NUMBER_OF_COMPUTATIONAL_NODES

    ENTERS("DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMPOSITION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MESH)) THEN

      INTERFACE_DECOMPOSITION=>INTERFACE_MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR

      IF(ASSOCIATED(INTERFACE_DECOMPOSITION)) THEN

        NULLIFY(INTERFACE_MESH%DECOMPOSITIONS)

        NUMBER_OF_COMPUTATIONAL_NODES = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
        IF(ERR/=0) GOTO 999

        ! Initialize and create a decomposition object.
        CALL DECOMPOSITIONS_INITIALISE(INTERFACE_MESH,ERR,ERROR,*999)

        CALL DECOMPOSITION_CREATE_START(INTERFACE_DECOMPOSITION%USER_NUMBER,INTERFACE_MESH,INTERFACE_DECOMPOSITION,ERR,ERROR,*999)

        CALL DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET(INTERFACE_DECOMPOSITION, .TRUE. , ERR, error, *999)

        CALL DECOMPOSITION_TYPE_SET(INTERFACE_DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)

        CALL DECOMPOSITION_NUMBER_OF_CONSTRAINTS_SET(INTERFACE_DECOMPOSITION, 2_INTG, err, error, *999)

        CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(INTERFACE_DECOMPOSITION,NUMBER_OF_COMPUTATIONAL_NODES,err,error,*999)

        CALL DECOMPOSITION_CREATE_FINISH(INTERFACE_DECOMPOSITION,ERR,ERROR,*999)

        INTERFACE_MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR=>INTERFACE_DECOMPOSITION
      ELSE
       CALL FlagError(" Interface decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
     CALL FlagError(" Interface Mesh  is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMPOSITION_SET")
   RETURN
999 ERRORSEXITS("DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMP",err,error)
   RETURN 1

  END SUBROUTINE DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMPOSITION_SET
!===================================================================================================================!
 !> The following subroutine adds the geometric field information of the interface graph G_{i}  in COUPLED_DECOMPOSITION object.

   SUBROUTINE COUPLED_DECOMPOSITION_ADD_INTERFACE(COUPLED_DECOMPOSITION,FIELD,ERR,ERROR,*)

  !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN)      :: COUPLED_DECOMPOSITION !<A pointer to the coupled decomposition object.
    TYPE(FIELD_TYPE), POINTER                                  :: FIELD !<A pointer to the field, representing geometric field of the coupled mesh graph.
    INTEGER(INTG), INTENT(OUT)                                 :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)                          :: ERROR !<The error string

    ! Local argument
    TYPE(VARYING_STRING)                                       :: LOCAL_ERROR

    ENTERS("COUPLED_DECOMPOSITION_ADD_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      IF(COUPLED_DECOMPOSITION%mesh_idx > 3) THEN
        LOCAL_ERROR="Coupled decomposition of user number "// &
          & TRIM(NUMBER_TO_VSTRING(COUPLED_DECOMPOSITION%USER_NUMBER,"*",ERR,ERROR))//&
          & " has number of assogned coupled meshes greater than 3."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR => FIELD
      END IF
    ELSE
     CALL FlagError(" Coupled decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("COUPLED_DECOMPOSITION_ADD_INTERFACE_MESH")
   RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_ADD_INTERFACE_MESH",err,error)
   RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_ADD_INTERFACE

!===================================================================================================================!
 !> The following subroutine adds the geometric field information of the coupled mesh graph G_{i} in COUPLED_DECOMPOSITION object.

   SUBROUTINE COUPLED_DECOMPOSITION_ADD_COUPLED_MESH(COUPLED_DECOMPOSITION,FIELD, ERR, ERROR,*)

  !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN)      :: COUPLED_DECOMPOSITION !<A pointer to the coupled decomposition object.
    TYPE(FIELD_TYPE), POINTER                                  :: FIELD !<A pointer to the field, representing geometric field of the coupled mesh graph.
    INTEGER(INTG), INTENT(OUT)                                 :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)                          :: ERROR !<The error string

    ! Local argument
    TYPE(VARYING_STRING)                                       :: LOCAL_ERROR


    ENTERS("COUPLED_DECOMPOSITION_ADD_COUPLED_MESH",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN

      IF(COUPLED_DECOMPOSITION%mesh_idx > 3) THEN
        LOCAL_ERROR="Coupled decomposition of user number "// &
          & TRIM(NUMBER_TO_VSTRING(COUPLED_DECOMPOSITION%USER_NUMBER,"*",ERR,ERROR))//&
          & " has number of assigned coupled meshes greater than 3."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR => FIELD
        COUPLED_DECOMPOSITION%mesh_idx = COUPLED_DECOMPOSITION%mesh_idx + 1
      END IF
    ELSE
     CALL FlagError(" Coupled decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("COUPLED_DECOMPOSITION_ADD_COUPLED_MESH")
   RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_ADD_COUPLED_MESH",err,error)
   RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_ADD_COUPLED_MESH
!===================================================================================================================!
!> The following subroutine initialize members of COUPLED_DECOMPOSITION object.

  SUBROUTINE COUPLED_DECOMPOSITION_CREATE_START(COUPLED_DECOMPOSITION, &
    & COUPLED_DECOMSPOSITION_USER_NUMBER, ERR, Error, *)

  !Argument variables
   TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)  :: COUPLED_DECOMPOSITION !<A pointer to the coupled decomposition object.
   INTEGER(INTG)                                             :: COUPLED_DECOMSPOSITION_USER_NUMBER !<A unique user number to indetify the  COUPLED_DECOMPOSITION object
   INTEGER(INTG), INTENT(OUT)                                :: ERR !<The error code
   TYPE(VARYING_STRING), INTENT(OUT)                         :: ERROR !<The error string
  !LOCAL Variables
   TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER                 :: NEW_COUPLED_DECOMPOSITION
   TYPE(VARYING_STRING)                                      :: LOCAL_ERROR

   ENTERS("COUPLED_DECOMPOSITION_CREATE_START",ERR,ERROR,*999)

   IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
     LOCAL_ERROR="Coupled Decomposition number "//&
       & TRIM(NUMBER_TO_VSTRING(COUPLED_DECOMSPOSITION_USER_NUMBER,"*",ERR,ERROR))// &
       & " has already been created."
     CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
   ELSE

     ALLOCATE(NEW_COUPLED_DECOMPOSITION, STAT=ERR)
     IF(ERR /=0) CALL FlagError(" Cannot allocate NEW_COUPLED_DECOMPOSITION.",ERR,ERROR,*999)

     NEW_COUPLED_DECOMPOSITION%mesh_idx = 1_INTG ! This member acts as an index for NEW_COUPLED_DECOMPOSITION%COUPLED_DECOMPOSITION(:).

     ALLOCATE(NEW_COUPLED_DECOMPOSITION%COUPLED_FIELDS(3),STAT=ERR) ! the first two indices store the geometric field of the coupled mesh objects and the 3rd index store the geometric field of the interface mesh object.
     IF(ERR /=0) CALL FlagError(" Cannot allocate COUPLED_MESH array.",ERR,ERROR,*999)

     NEW_COUPLED_DECOMPOSITION%USER_NUMBER=COUPLED_DECOMSPOSITION_USER_NUMBER

     COUPLED_DECOMPOSITION=>NEW_COUPLED_DECOMPOSITION

   ENDIF
   EXITS("COUPLED_DECOMPOSITION_CREATE_START")
   RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_CREATE_START",err,error)
   RETURN 1

  END SUBROUTINE COUPLED_DECOMPOSITION_CREATE_START
!====================================================================================!
 !> In the following subroutine the "coupling aware" mesh partitioning is implemented.


  SUBROUTINE COUPLED_DECOMPOSITION_CREATE_FINISH(COUPLED_DECOMPOSITION, ERR, Error, *)

   !Argument variables
   TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)  :: COUPLED_DECOMPOSITION !<Coupled decomposition object.
   INTEGER(INTG), INTENT(OUT)                                :: ERR !<The error code.
   TYPE(VARYING_STRING), INTENT(OUT)                         :: ERROR !<The error string.
   !LOCAL Variables
   TYPE(MESH_TYPE), POINTER   :: COUPLED_MESH_1,COUPLED_MESH_2,INTERFACE_MESH
   TYPE(FIELD_TYPE), POINTER  :: FIELD_COUPLED_MESH_1,FIELD_COUPLED_MESH_2,FIELD_INTERFACE_MESH
   TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION_COUPLED_MESH_1, DECOMPOSITION_COUPLED_MESH_2, &
     & DECOMPOSITION_INTERFACE_MESH
   INTEGER(INTG)  ::  MY_COMPUTATIONAL_NODE


   ENTERS("COUPLED_DECOMPOSITION_CREATE_START",ERR,ERROR,*999)


   IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
     FIELD_COUPLED_MESH_1=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(1)%PTR

     IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(1)%PTR)) THEN
       FIELD_COUPLED_MESH_2=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(2)%PTR

       IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(2)%PTR)) THEN
          FIELD_INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR

          IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR)) THEN
            DECOMPOSITION_COUPLED_MESH_1=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(1)%PTR%DECOMPOSITION

            IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(1)%PTR%DECOMPOSITION)) THEN
              DECOMPOSITION_COUPLED_MESH_2=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(2)%PTR%DECOMPOSITION

              IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(2)%PTR%DECOMPOSITION)) THEN
                DECOMPOSITION_INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION

                IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION)) THEN
                  COUPLED_MESH_1=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(1)%PTR%DECOMPOSITION%MESH

                  IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(1)%PTR%DECOMPOSITION%MESH)) THEN
                    COUPLED_MESH_2=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(2)%PTR%DECOMPOSITION%MESH

                    IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(2)%PTR%DECOMPOSITION%MESH)) THEN
                      INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH

                      IF(ASSOCIATED(COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH)) THEN


                        MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     ! Get the processor Id.

                        CALL GATHER_MESH_COORDINATES(FIELD_INTERFACE_MESH, COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES, &
                          & ERR, Error, *999) ! Gather interface mesh coordinates

                        ! The following subroutine trivially decomposes the interface mesh.
                        CALL DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMPOSITION_SET(INTERFACE_MESH, ERR, ERROR, *999)

                        COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION=> &
                          & INTERFACE_MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR

                        COUPLED_DECOMPOSITION%mesh_idx = 1 !resetting the index to 1.

                        CALL COUPLED_DECOMPOSITION_FIXED_VERTEX_PARTITIONING(COUPLED_DECOMPOSITION,ERR,ERROR,*999)
                        CALL COUPLED_DECOMPOSITION_FIXED_VERTEX_PARTITIONING(COUPLED_DECOMPOSITION,ERR,ERROR,*999)

                        COUPLED_DECOMPOSITION%mesh_idx = 1  !Reset the value

                        ELSE
                          CALL FlagError("Interface mesh is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Coupled mesh 2 is not associated.",ERR,ERROR,*999)
                      END IF
                   ELSE
                     CALL FlagError("Coupled mesh 1 is not associated.",ERR,ERROR,*999)
                   END IF
                 ELSE
                   CALL FlagError("Decomposition of interface mesh is not associated.",ERR,ERROR,*999)
                 END IF
              ELSE
                 CALL FlagError("Decomposition of coupled mesh 1 is not associated.",ERR,ERROR,*999)
              END IF
            ELSE
              CALL FlagError("Decomposition of coupled mesh 2 is not associated.",ERR,ERROR,*999)
            END IF
         ELSE
           CALL FlagError("Geometric Field of Interface mesh is not associated.",ERR,ERROR,*999)
         END IF
       ELSE
         CALL FlagError("Geometric Field of coupled mesh 2 is not associated.",ERR,ERROR,*999)
       END IF
     ELSE
       CALL FlagError("Geometric Field of coupled mesh 1 is not associated.",ERR,ERROR,*999)
     END IF
   ELSE
     CALL FlagError("Coupled decomposition is not associated.",ERR,ERROR,*999)
   END IF

   EXITS("COUPLED_DECOMPOSITION_CREATE_FINISH.")
   RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_CREATE_FINISH.",ERR,ERROR)
   RETURN 1

  END SUBROUTINE COUPLED_DECOMPOSITION_CREATE_FINISH

!=========================================================================================================================
 !> The following subroutine gathers the geometric field information form all processors in one data structure MESH_COORDINATES(node_idx,:),
 !> i.e. x = MESH_COORDINATES(node_idx,1),  y = MESH_COORDINATES(node_idx,2), z = MESH_COORDINATES(node_idx,3).

  SUBROUTINE GATHER_MESH_COORDINATES(GEOMETRIC_FIELD, MESH_COORDINATES,  ERR, Error, * )

   !Argument variables
   TYPE(FIELD_TYPE), POINTER, INTENT(INOUT)             :: GEOMETRIC_FIELD ! A pointer to the geometric field of the mesh whose coordinates are to be gathered.
   REAL(RP), INTENT(OUT), ALLOCATABLE                   :: MESH_COORDINATES(:,:) ! Data structure where the geometric field coordinates are gathered.
   INTEGER(INTG), INTENT(OUT)                           :: ERR !<The error code
   TYPE(VARYING_STRING), INTENT(OUT)                    :: ERROR !<The error string

   !Local variables
   INTEGER(INTG)                                        :: component_idx, DOMAIN_NUMBER, node_idx, NUMBER_OF_FIELD_COMPONENTS, &
     & MY_COMPUTATIONAL_NODE
   TYPE(MESH_TYPE) , POINTER                            :: MESH
   TYPE(DECOMPOSITION_TYPE) , POINTER                   :: DECOMPOSITION
   TYPE(FIELD_VARIABLE_TYPE), POINTER                   :: FIELD_VARIABLE_MESH
   REAL(RP), ALLOCATABLE                                :: MESH_COORDINATES_NEW(:)

   ENTERS("GATHER_MESH_COORDINATES(MESH, MESH_COORDIANTES",ERR,ERROR,*999)

   IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
     DECOMPOSITION=>GEOMETRIC_FIELD%DECOMPOSITION
     IF(ASSOCIATED(DECOMPOSITION)) THEN
       MESH=>DECOMPOSITION%MESH
       IF(ASSOCIATED(MESH)) THEN
         FIELD_VARIABLE_MESH=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(1)%PTR ! one because every geometric field has one variable
         IF(ASSOCIATED(FIELD_VARIABLE_MESH)) THEN

           MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
           NUMBER_OF_FIELD_COMPONENTS = GEOMETRIC_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS !NUmber of components of a geometric field.
           ALLOCATE(MESH_COORDINATES(MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes,3), STAT=ERR)
           IF(ERR /= 0)   CALL FlagError(" Cannot allocate MESH_COORDINATES.",ERR,ERROR,*999)
           MESH_COORDINATES=0
           DO node_idx = 1, MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes

             CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,node_idx,1_INTG, &
               & DOMAIN_NUMBER,ERR,ERROR,*999)

             IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE)  THEN

               DO component_idx=1,NUMBER_OF_FIELD_COMPONENTS
                 CALL FIELD_PARAMETER_SET_GET_NODE(GEOMETRIC_FIELD,FIELD_VARIABLE_MESH%VARIABLE_TYPE, &
                   & 1_INTG, 1_INTG, 1_INTG, &
                   & node_idx,component_idx,MESH_COORDINATES(node_idx,component_idx),ERR,ERROR,*999)
               END DO ! component_idx
             END IF ! IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE)
           END DO ! node_idx

           ALLOCATE(MESH_COORDINATES_NEW(MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes*3), STAT=ERR)
           IF(ERR /= 0)   CALL FlagError(" Cannot allocate MESH_COORDINATES_NEW.",ERR,ERROR,*999)
           MESH_COORDINATES_NEW = 0.


           CALL MPI_ALLREDUCE(MESH_COORDINATES, &
             & MESH_COORDINATES_NEW, MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes*3, &
             & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,Err)

           MESH_COORDINATES = 0.
           MESH_COORDINATES = RESHAPE(MESH_COORDINATES_NEW, (/MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes,3/))

           DEALLOCATE(MESH_COORDINATES_NEW)

           ELSE
             CALL FlagError(" Field variable is not associated.",ERR,ERROR,*999)
           ENDIF
         ELSE
           CALL FlagError(" Mesh is not associated.",ERR,ERROR,*999)
         ENDIF
       ELSE
         CALL FlagError(" Decomposition is not associated.",ERR,ERROR,*999)
       ENDIF
     ELSE
       CALL FlagError(" Geometric field is not associated.",ERR,ERROR,*999)
     ENDIF
     RETURN
999 ERRORSEXITS("GATHER_MESH_COORDINATES.",ERR,ERROR)
    RETURN 1

  END SUBROUTINE GATHER_MESH_COORDINATES

!========================================================================================================
  !>The following subroutine updates the new geometric field.

  SUBROUTINE DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD(FIELD,DECOMPOSITION,VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER, INTENT(IN)         :: FIELD !<A pointer to the field to update the geometric parameters for
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: DECOMPOSITION !<The mesh which is generated by the generated mesh \todo is this necessary???
    INTEGER(INTG), INTENT(IN)  :: VARIABLE_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER           :: MESH
    INTEGER(INTG)                      :: DOMAIN,MY_COMPUTATIONAL_NODE, node_idx, TOTAL_NODES

    ENTERS("DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
          TOTAL_NODES = SIZE(MESH%TOPOLOGY(1)%PTR%NODES%NODES)
          DO node_idx = 1, TOTAL_NODES
            CALL DECOMPOSITION_NODE_DOMAIN_GET( DECOMPOSITION, node_idx, 1_INTG, DOMAIN, ERR, ERROR, *999)
            IF(DOMAIN==MY_COMPUTATIONAL_NODE) THEN
              CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,VARIABLE_TYPE,1_INTG,1_INTG,1_INTG, &
                & node_idx,1, DECOMPOSITION%NODE_DOMAIN(node_idx-1), ERR,Error,*999)
            END IF
          END DO !node_idx
        ELSE
          CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD")

    RETURN
999 ERRORSEXITS("DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD",ERR,ERROR)
    RETURN 1

  END SUBROUTINE DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD

!========================================================================================================
  !>The following subroutine updates the decomposition field of the coupled mesh graph G_i.


  SUBROUTINE COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION(COUPLED_DECOMPOSITION,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)         :: DECOMPOSITION !< A decomposition type object to be updated.
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !<Coupled decomposition type object that stores information about the coupled mesh partitioning.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER        :: MESH
    INTEGER(INTG)                   :: mapping_idx, NUMBER_OF_COMPUTATIONAL_NODES
    TYPE(DECOMPOSITION_TYPE), POINTER  :: NEW_DECOMPOSITION

    ENTERS("COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          NUMBER_OF_COMPUTATIONAL_NODES = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
          IF(ERR/=0) GOTO 999
          NULLIFY(MESH%DECOMPOSITIONS)

          CALL GET_SURROUNDING_NODES(MESH, ERR, error, *999)

          CALL DECOMPOSITIONS_INITIALISE(MESH,ERR,ERROR,*999)

          CALL DECOMPOSITION_CREATE_START(DECOMPOSITION%USER_NUMBER,MESH,NEW_DECOMPOSITION,ERR,ERROR,*999)

          CALL DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET(NEW_DECOMPOSITION, .TRUE. , ERR, error, *999)

          CALL DECOMPOSITION_TYPE_SET(NEW_DECOMPOSITION,DECOMPOSITION_USER_DEFINED_TYPE,ERR,ERROR,*999)

          CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(NEW_DECOMPOSITION,NUMBER_OF_COMPUTATIONAL_NODES,err,error,*999)

          DO mapping_idx = 1, MESH%TOPOLOGY(1)%PTR%NODES%numberOfNodes

            NEW_DECOMPOSITION%NODE_DOMAIN(mapping_idx-1) = &
              & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%NODE_DOMAIN &
              & (mapping_idx-1)

          END DO !mapping_idx

          CALL DECOMPOSITION_CREATE_FINISH(NEW_DECOMPOSITION, ERR, Error, *999)

          COUPLED_DECOMPOSITION%mesh_idx = COUPLED_DECOMPOSITION%mesh_idx + 1

          DECOMPOSITION=>NEW_DECOMPOSITION
        ELSE
          CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Coupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION")
    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION


!========================================================================================================
 !>The following subroutine updates the decomposition field of the interface mesh graph G_I.

  SUBROUTINE COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION(COUPLED_DECOMPOSITION,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)         :: DECOMPOSITION !<A pointer to the field to update the geometric parameters for
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !<The mesh which is generated by the generated mesh \todo is this necessary???
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        NULLIFY(DECOMPOSITION)
        DECOMPOSITION=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION
      ELSE
        CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Coupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION")
    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DEC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION


!====================================================================================================
 !> The following subroutine identifies the vertices to impose the fixed vertex partitioning o.

 SUBROUTINE COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING(COUPLED_DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !< Coupled decomposition type objects.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG)                      :: inter_edge_idx, MY_COMPUTATIONAL_NODES, proc_idx, PROC_ID, TOTAL_COMPUTATIONAL_NODES
    INTEGER(INTG), ALLOCATABLE         :: sub_domain_counter(:), TEMP_ARRAY(:,:)
    TYPE(DECOMPOSITION_TYPE), POINTER  :: INTERFACE_MESH_DECOMPOSITION

    ENTERS("COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN

      INTERFACE_MESH_DECOMPOSITION=> &
        & COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION

      IF(ASSOCIATED(INTERFACE_MESH_DECOMPOSITION)) THEN

        TOTAL_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR) ! Get rank.
        MY_COMPUTATIONAL_NODES=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) ! Get id of the processor

        ALLOCATE(sub_domain_counter(TOTAL_COMPUTATIONAL_NODES), STAT=ERR)
        IF(ERR/=0)   CALL FlagError("Sub_domain_countercannot be allocated.",ERR,ERROR,*999)

        ALLOCATE(TEMP_ARRAY(SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1),TOTAL_COMPUTATIONAL_NODES), STAT=ERR) !TEMP_ARRAY(:,PROC_idx) contains all the vertices that are assigned sub-domain PROC_idx
        IF(ERR/=0)   CALL FlagError("TEMP_ARRAY cannot be allocated.",ERR,ERROR,*999)

        ! This subroutine changes the values of COUPLED_DECOMPOSITION%INTER_EDGES such that the new...
        ! structure of COUPLED_DECOMPOSITION%INTER_EDGES will consist of NUMBER_OF_COMPUTATIONAL_NODES...
        ! number of columns and COUPLED_DECOMPOSITION%INTER_EDGES(:,subdomain_idx) will consist of all the vertices of .. 
        ! the coupled mesh graph G_i that are suppose to belong to subdomain_idx.
  
        TEMP_ARRAY =  0
        sub_domain_counter= 1

        DO proc_idx = 1, TOTAL_COMPUTATIONAL_NODES

          DO inter_edge_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)

            PROC_ID = INTERFACE_MESH_DECOMPOSITION%NODE_DOMAIN(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,1)-1)

            IF(PROC_ID == proc_idx-1) THEN

              TEMP_ARRAY(sub_domain_counter(proc_idx),proc_idx)=COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2)

              sub_domain_counter(proc_idx) = sub_domain_counter(proc_idx) + 1

            END IF

          END DO !inter_edge_idx

        END DO !proc_idx

        DEALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES) ! Store the vertices in COUPLED_DECOMPOSITION%INTER_EDGES data structure ...
                                                      !... such that !COUPLED_DECOMPOSITION%INTER_EDGES(:,PROC_idx) contains all  ... 
                                                      !.... the vertices that are assigned sub-domain PROC_idx
        ALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES(SIZE(TEMP_ARRAY,1),TOTAL_COMPUTATIONAL_NODES))
        COUPLED_DECOMPOSITION%INTER_EDGES = TEMP_ARRAY
        DEALLOCATE(TEMP_ARRAY)
        DEALLOCATE(sub_domain_counter)
      ELSE
        CALL FlagError("Interface mesh decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Coupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING")
    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DEC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING

!=====================================================================================================================
  !>This subroutine builds mapping between numbering of a vertices of the new graph G_{i,merged} and the orginal graph G_i.

 SUBROUTINE COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING(COUPLED_DECOMPOSITION, NEW_TO_OLD_INDEX_MAPPING, ERR,ERROR,*)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !<Coupled decomposition type objects.
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: NEW_TO_OLD_INDEX_MAPPING(:,:) !<NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,:). Contains one-to-one mapping between new vertex ids and old vertex ids
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG)    :: inter_edge_idx,NUMBER_OF_COMPUTATIONAL_NODES,new_vertex_idx, old_vertex_idx, proc_idx
    LOGICAL          :: FLAG1, FLAG2
    TYPE(DECOMPOSITION_TYPE), POINTER  :: COUPLED_MESH_DECOMPOSITION
    TYPE(MESH_TYPE), POINTER  :: COUPLED_MESH

    ENTERS("COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN

      COUPLED_MESH_DECOMPOSITION=> &
        & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION
      IF(ASSOCIATED(COUPLED_MESH_DECOMPOSITION)) THEN

        COUPLED_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%MESH

        IF(ASSOCIATED(COUPLED_MESH)) THEN

          NUMBER_OF_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)

          ALLOCATE(NEW_TO_OLD_INDEX_MAPPING(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNOdes,2), STAT=ERR)
          new_vertex_idx = 1

          DO old_vertex_idx = 1, COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNOdes

            FLAG1= .FALSE.
            FLAG2= .FALSE.
            NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,1) = old_vertex_idx
  
            ! proc_idx represents a sub-domain. 
            DO proc_idx = 1, NUMBER_OF_COMPUTATIONAL_NODES 

              DO inter_edge_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)

                IF(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,proc_idx)==old_vertex_idx) THEN !If old_vertex_idx is one of the vertices to be merged
                  FLAG2= .TRUE.
                  IF(inter_edge_idx==1) THEN ! In this case the Id would of old_vertex_idx  would be the same as the id of the first merged vertex that is supposed to go in ...
                    FLAG1=.FALSE.  !.... sub-domain proc_idx. 
                  ELSE
                    FLAG1=.TRUE.
                  END IF
                  EXIT
                END IF
              END DO

              IF(FLAG2) EXIT ! Exit outer loop if old_vertex_idx is detected as one of the vertices to be merged.

            END DO

            IF(.NOT. FLAG1) THEN ! If old_vertex_idx is not the merged vertex then assign old_vertex_idx a new Id.

              NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,2)=new_vertex_idx
              new_vertex_idx = new_vertex_idx + 1

            ELSE ! If old_vertex_idx  is suppose to be a merged vertex belonging to sub-domain proc_idx then .... 
                 ! ... assign old_vertex_idx the Id of the new id of the merged vertex.
  
              NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,2)= &
                & NEW_TO_OLD_INDEX_MAPPING(COUPLED_DECOMPOSITION%INTER_EDGES(1,proc_idx),2)

            END IF
          END DO ! OLD_VERTEX_idx

        ELSE
          CALL FlagError("Coupled mesh decomposition is not associated.",ERR,ERROR,*999)
        ENDIF  
      ELSE
        CALL FlagError("Coupled mesh decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Coupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING")
    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING


! =========================================================================================================!
  !>The following subroutine build adjacencies of the new coupled mesh graph G_{i,merged}.

  SUBROUTINE COUPLED_DECOMPOSITION_GET_NEW_GRAPH(COUPLED_DECOMPOSITION, COUPLED_MESH, &
   & NEW_TO_OLD_INDEX_MAPPING, ERR, ERROR,* )

    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER :: COUPLED_DECOMPOSITION !<Coupled decomposition type objects.
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT) :: COUPLED_MESH  !<Mesh type object.
    INTEGER(INTG), INTENT(IN)  :: NEW_TO_OLD_INDEX_MAPPING(:,:) !<NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,:). Cotains one-to-one mapping between new vertex ids and old vertex ids
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG)    :: index_start, index_end, inter_edge_idx, idx, LOC, MY_COMPUTATIONAL_NODE, &
      & NUMBER_OF_COMPUTATIONAL_NODES, node_idx, NODES_TO_RETAIN, NODES_TO_COLLAPSE, NUMBER_OF_ROWS_TO_DELETE, &
      & proc_idx, surrounding_node_idx, SURROUNDING_NODE
    INTEGER(INTG), ALLOCATABLE    :: TEMP_ARRAY(:), TEMP_ARRAY_2D(:,:)


    ENTERS("COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      IF(ASSOCIATED(COUPLED_MESH)) THEN
        NUMBER_OF_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
        MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)

        DO proc_idx = 1,  NUMBER_OF_COMPUTATIONAL_NODES

          NODES_TO_RETAIN = COUPLED_DECOMPOSITION%INTER_EDGES(1,proc_idx)

          ALLOCATE(TEMP_ARRAY(SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)* &
            & SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)), STAT=ERR)

          TEMP_ARRAY = 0
          index_start = 0

          TEMP_ARRAY(1:SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_RETAIN)%surroundingNOdes,1)) = &
            &  COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_RETAIN)%surroundingNOdes
          ! In the Step 1 add adjacencies of the nodes with common sub-domain.
          DO inter_edge_idx = 2 , SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)

            NODES_TO_COLLAPSE = COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,proc_idx)

            IF(NODES_TO_COLLAPSE/=0) THEN

              index_start=index_start + SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_RETAIN)%surroundingNOdes,1)+1

              index_end = index_start + SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_COLLAPSE)%surroundingNOdes,1)-1

              TEMP_ARRAY(index_start:index_end) = COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_COLLAPSE)%surroundingNOdes

            END IF

          END DO  !inter_edge_idx

          DEALLOCATE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_RETAIN)%surroundingNOdes)

          ALLOCATE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_RETAIN)%surroundingNOdes &
            & (SIZE(PACK(TEMP_ARRAY, TEMP_ARRAY /= 0 ),1)), STAT=ERR)

          COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(NODES_TO_RETAIN)%surroundingNOdes = PACK(TEMP_ARRAY, TEMP_ARRAY/=0)

          DEALLOCATE(TEMP_ARRAY)

        END DO !proc_idx
        ! In step 2, try to nullify nodes that are supposed to be merged.
        DO inter_edge_idx = 2 , SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)

          DO proc_idx = 1, NUMBER_OF_COMPUTATIONAL_NODES

            node_idx = COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,proc_idx)
            IF(node_idx/=0) THEN
              COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes = -1
            END IF
          END DO !proc_idx
        END DO !inter_edge_idx

        NUMBER_OF_ROWS_TO_DELETE = 0

        DO node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)

          IF(SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes) > 0) THEN
            IF(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(1) == -1) THEN

              NUMBER_OF_ROWS_TO_DELETE = NUMBER_OF_ROWS_TO_DELETE + 1     !No. of rows less in the new coupled mesh graph G_{i,merged} compared to the orginal coupled mesh graph G_i.

            END IF
          END IF

        END DO    !node_idx

        IF(ALLOCATED(TEMP_ARRAY_2D)) DEALLOCATE(TEMP_ARRAY_2D)
        ALLOCATE(TEMP_ARRAY_2D(SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)-NUMBER_OF_ROWS_TO_DELETE, &
          & SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)*10), STAT=ERR)
        IF(ERR/=0) CALL FlagError("Unable to allocate TEMP_ARRAY_2D array.",ERR,ERROR,*999)

        TEMP_ARRAY_2D = 0 ! THis array temporarily stores the new node adjacencies.
        idx = 1

        DO node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)

          IF(SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes) > 0) THEN
            IF(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(1) /= -1) THEN

              TEMP_ARRAY_2D(idx,1:SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes,1)) = &
                & COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes

              idx = idx + 1
            END IF
          END IF
        END DO  !node_idx

        DEALLOCATE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES)
        ALLOCATE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(SIZE(TEMP_ARRAY_2D,1)))

        !  Move the new vertex/node adjacencies to the surroundingNOdes data structure.
        DO node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)

          ALLOCATE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes( &
            & SIZE(PACK(TEMP_ARRAY_2D(node_idx,:),TEMP_ARRAY_2D(node_idx,:) /=0))))

          COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(1: &
            & SIZE(PACK(TEMP_ARRAY_2D(node_idx,:),TEMP_ARRAY_2D(node_idx,:) /=0))) = &
            & TEMP_ARRAY_2D(node_idx,1:SIZE(PACK(TEMP_ARRAY_2D(node_idx,:),TEMP_ARRAY_2D(node_idx,:) /=0)))

        END DO  !node_idx

       ! IN step 3 change the vertex ids according to the information stored in NEW_TO_OLD_INDEX_MAPPING(:,:) array.
        DO node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)

          DO surrounding_node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes)

            LOC = COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(surrounding_node_idx)

            COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(surrounding_node_idx)= &
              & NEW_TO_OLD_INDEX_MAPPING(LOC,2)

          END DO !surrounding_node_idx

        END DO ! node_idx

        ! In step 4 remove the node ids that are adjacent to themselves. For instance if adjacency of node 1 is [1,2,3,4] then remove node 1 in the adjacency vector as there is no edge between a node and itself. Therefore the new node adjacency will be [2,3,4].
        DO node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)

          DO surrounding_node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes)

            SURROUNDING_NODE = COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(surrounding_node_idx)

            IF(node_idx == SURROUNDING_NODE) THEN

              COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(surrounding_node_idx)= 0

            END IF

          END DO !surrounding_node_idx

        END DO ! node_idx

        DO node_idx = 1, SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1)

          ALLOCATE(TEMP_ARRAY(SIZE(PACK(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes,&
            & COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes /=0))))

          TEMP_ARRAY = PACK(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes,&
            & COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes /=0)

          DEALLOCATE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes)

          ALLOCATE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes(SIZE(TEMP_ARRAY,1)))

          COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(node_idx)%surroundingNOdes = TEMP_ARRAY

          IF(ALLOCATED(TEMP_ARRAY))  DEALLOCATE(TEMP_ARRAY)
          IF(ALLOCATED(TEMP_ARRAY_2D)) DEALLOCATE(TEMP_ARRAY_2D)
        END DO
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    ELSE
      CALL FlagError("Coupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_DECOMPOSITION_GET_NEW_GRAPH")
    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_GET_NEW_GRAPH",ERR,ERROR)
    RETURN 1

  END SUBROUTINE COUPLED_DECOMPOSITION_GET_NEW_GRAPH

END MODULE REGION_ROUTINES

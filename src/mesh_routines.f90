!> \file
!> \author Chris Bradley
!> \brief This module handles all mesh (node and element) routines.
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

!> This module handles all mesh (node and element) routines.
MODULE MESH_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS_MPI
  USE CMISS_PARMETIS
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE DOMAIN_MAPPINGS
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE LISTS
#ifndef NOMPIMOD
  USE MPI
#endif
  USE NODE_ROUTINES
  USE STRINGS
  USE TREES
  USE TYPES
  USE FIELD_ROUTINES !MIRZAWD
  USE REGION_ROUTINES
  USE COORDINATE_ROUTINES
  USE NODE_ROUTINES

#include "macros.h"  

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PRIVATE
  !Module parameters

  !> \addtogroup MESH_ROUTINES_DecompositionTypes MESH_ROUTINES::DecompositionTypes
  !> \brief The Decomposition types parameters
  !> \see MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_ALL_TYPE=1 !<The decomposition contains all elements. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_CALCULATED_TYPE=2 !<The element decomposition is calculated by graph partitioning. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_USER_DEFINED_TYPE=3 !<The user will set the element decomposition. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Starts the process of creating a mesh
  INTERFACE MESH_CREATE_START
    MODULE PROCEDURE MESH_CREATE_START_INTERFACE
    MODULE PROCEDURE MESH_CREATE_START_REGION
  END INTERFACE !MESH_CREATE_START

  !>Initialises the meshes for a region or interface.
  INTERFACE MESHES_INITIALISE
    MODULE PROCEDURE MESHES_INITIALISE_INTERFACE
    MODULE PROCEDURE MESHES_INITIALISE_REGION
  END INTERFACE !MESHES_INITIALISE

  INTERFACE MESH_USER_NUMBER_FIND
    MODULE PROCEDURE MESH_USER_NUMBER_FIND_INTERFACE
    MODULE PROCEDURE MESH_USER_NUMBER_FIND_REGION
  END INTERFACE !MESH_USER_NUMBER_FIND

  INTERFACE MeshTopologyNodeCheckExists
    MODULE PROCEDURE MeshTopologyNodeCheckExistsMesh
    MODULE PROCEDURE MeshTopologyNodeCheckExistsMeshNodes    
  END INTERFACE MeshTopologyNodeCheckExists  
  
  INTERFACE MeshTopologyElementCheckExists
    MODULE PROCEDURE MeshTopologyElementCheckExistsMesh
    MODULE PROCEDURE MeshTopologyElementCheckExistsMeshElements
  END INTERFACE MeshTopologyElementCheckExists
  
  PUBLIC DECOMPOSITION_ALL_TYPE,DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE

  PUBLIC DECOMPOSITIONS_INITIALISE,DECOMPOSITIONS_FINALISE

  PUBLIC DECOMPOSITION_CREATE_START,DECOMPOSITION_CREATE_FINISH

  PUBLIC DECOMPOSITION_DESTROY

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_GET,DECOMPOSITION_ELEMENT_DOMAIN_SET

  PUBLIC DECOMPOSITION_MESH_COMPONENT_NUMBER_GET,DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  PUBLIC DECOMPOSITION_NUMBER_OF_DOMAINS_GET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET
  
  PUBLIC DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS,DecompositionTopology_DataPointCheckExists
  
  PUBLIC DecompositionTopology_DataProjectionCalculate

  PUBLIC DecompositionTopology_ElementDataPointLocalNumberGet

  PUBLIC DecompositionTopology_ElementDataPointUserNumberGet

  PUBLIC DecompositionTopology_NumberOfElementDataPointsGet
  
  PUBLIC DECOMPOSITION_TYPE_GET,DECOMPOSITION_TYPE_SET
  
  PUBLIC DECOMPOSITION_USER_NUMBER_FIND, DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION
  
  PUBLIC DECOMPOSITION_NODE_DOMAIN_GET

  PUBLIC DECOMPOSITION_CALCULATE_LINES_SET,DECOMPOSITION_CALCULATE_FACES_SET

  PUBLIC DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS

  PUBLIC DomainTopology_ElementBasisGet
  
  PUBLIC MeshTopologyElementCheckExists,MeshTopologyNodeCheckExists
  
  PUBLIC MESH_CREATE_START,MESH_CREATE_FINISH

  PUBLIC MESH_DESTROY
  
  PUBLIC MESH_NUMBER_OF_COMPONENTS_GET,MESH_NUMBER_OF_COMPONENTS_SET

  PUBLIC MESH_NUMBER_OF_ELEMENTS_GET,MESH_NUMBER_OF_ELEMENTS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_CREATE_START,MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH

  PUBLIC MESH_TOPOLOGY_ELEMENTS_DESTROY

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET,MeshElements_ElementNodeVersionSet

  PUBLIC MESH_TOPOLOGY_ELEMENTS_GET

  PUBLIC MeshElements_ElementUserNumberGet,MeshElements_ElementUserNumberSet
  
  PUBLIC MeshTopologyElementsUserNumbersAllSet

  PUBLIC MeshTopologyDataPointsCalculateProjection

  PUBLIC MeshTopologyNodeDerivativesGet

  PUBLIC MeshTopologyNodeNumberOfDerivativesGet

  PUBLIC MeshTopologyNodeNumberOfVersionsGet

  PUBLIC MeshTopologyNodesNumberOfNodesGet

  PUBLIC MeshTopologyNodesDestroy
  
  PUBLIC MeshTopologyNodesGet

  PUBLIC MESH_USER_NUMBER_FIND, MESH_USER_NUMBER_TO_MESH

  PUBLIC MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  PUBLIC MESH_EMBEDDING_CREATE,MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  PUBLIC MESH_EMBEDDING_SET_GAUSS_POINT_DATA

  PUBLIC MESHES_INITIALISE,MESHES_FINALISE

  ! from here mizawds contributions starts
  
  PUBLIC DECOMPOSITION_TPWGT_SET, DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET, &
   & DECOMPOSITION_UBVEC_SET, DECOMPOSITION_NUMBER_OF_CONSTRAINTS_SET,&
     & DECOMPOSITION_NODE_WEIGHT_SET, GET_SURROUNDING_NODES, &
       & GET_SURROUNDING_NODES_LINEAR_QUADRILATERAL_MESH, COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A,&
         & COUPLED_DECOMPOSITION_CREATE_START, COUPLED_DECOMPOSITION_ADD_COUPLED_MESH, COUPLED_DECOMPOSITION_ADD_INTERFACE, &
           & COUPLED_DECOMPOSITION_CREATE_FINISH, DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD, &
             & COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION, COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION  


CONTAINS

  !
  !================================================================================================================================
  !

  !>Finialises a decomposition adjacent element information and deallocates all memory
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
    IF(ALLOCATED(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises a decomposition adjacent element information.
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition on a given mesh. \see OPENCMISS::CMISSDecompositionCreateFinish
  SUBROUTINE DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(MESH_TYPE), POINTER :: MESH

    ENTERS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN

      IF(DECOMPOSITION%NODE_BASED_DECOMPOSITION .EQV. .FALSE.) THEN

        CALL DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)


      ELSE IF(DECOMPOSITION%NODE_BASED_DECOMPOSITION .EQV. .TRUE. & 
        & .AND. DECOMPOSITION%DECOMPOSITION_TYPE == DECOMPOSITION_CALCULATED_TYPE) THEN 

        IF( .NOT. ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(1)%PTR%NODES%NODES(1)%surroundingNodes)) & 
          & CALL GET_SURROUNDING_NODES(DECOMPOSITION%MESH, ERR, error, *999) !check to see if surroundingNodes data structure is  already calculated for the given mesh graph. 

        CALL DECOMPOSITION_NODE_BASED_DECOMPOSITION_PARAMETERS_SET(DECOMPOSITION,ERR,ERROR,*999) 

        CALL DECOMPOSITION_NODE_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)

      END IF

      !Initialise the topology information for this decomposition
      CALL DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the domain for this computational node
      CALL DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Calculate the decomposition topology

      CALL DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      DECOMPOSITION%DECOMPOSITION_FINISHED=.TRUE.
      ! 
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      MESH=>DECOMPOSITION%MESH
      IF(ASSOCIATED(MESH)) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Mesh = ",MESH%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of decompositions = ", &
          & MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS,ERR,ERROR,*999)
        DO decomposition_no=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Decomposition number = ",decomposition_no,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%GLOBAL_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%USER_NUMBER,ERR,ERROR,*999)
        ENDDO !decomposition_no
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a given mesh. \see OPENCMISS::CMISSDecompositionCreateStart
  SUBROUTINE DECOMPOSITION_CREATE_START(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition 
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to decompose
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On return a pointer to the created decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no, number_computational_nodes
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: NEW_DECOMPOSITION
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITION)
    NULLIFY(NEW_DECOMPOSITIONS)

    ENTERS("DECOMPOSITION_CREATE_START",ERR,ERROR,*999)

    NULLIFY(DECOMPOSITION)
    
    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
          IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
            CALL DECOMPOSITION_USER_NUMBER_FIND(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
            IF(ASSOCIATED(DECOMPOSITION)) THEN
              LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
                & " has already been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              !\todo Split this into an initialise and create start.
              ALLOCATE(NEW_DECOMPOSITION,STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decomposition.",ERR,ERROR,*999)
              ! extracting the number of computational nodes
              number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              !Set default decomposition properties
              NEW_DECOMPOSITION%GLOBAL_NUMBER=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1
              NEW_DECOMPOSITION%USER_NUMBER=USER_NUMBER
              NEW_DECOMPOSITION%DECOMPOSITION_FINISHED=.FALSE.
              NEW_DECOMPOSITION%CALCULATE_LINES=.TRUE. !Default
              NEW_DECOMPOSITION%CALCULATE_FACES=.FALSE. !Default
              NEW_DECOMPOSITION%DECOMPOSITIONS=>MESH%DECOMPOSITIONS
              NEW_DECOMPOSITION%MESH=>MESH
              !By default, the process of decompostion was done on the first mesh components. But the decomposition is the same for all mesh components, since the decomposition is element-based.
              NEW_DECOMPOSITION%MESH_COMPONENT_NUMBER=1
              !Default decomposition is all the mesh with one domain.
              NEW_DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
              NEW_DECOMPOSITION%NUMBER_OF_DOMAINS = number_computational_nodes  
              NEW_DECOMPOSITION%NODE_BASED_DECOMPOSITION = .FALSE.

              ALLOCATE(NEW_DECOMPOSITION%TPWGT(NEW_DECOMPOSITION%NUMBER_OF_DOMAINS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decomposition tpwgts.",ERR,ERROR,*999)

              ALLOCATE(NEW_DECOMPOSITION%UBVEC(1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decomposition ubvec.",ERR,ERROR,*999)

              NEW_DECOMPOSITION%TPWGT = REAL(1.0/(NEW_DECOMPOSITION%NUMBER_OF_DOMAINS),DP)
              NEW_DECOMPOSITION%NUMBER_OF_CONSTRAINTS=1_INTG
              NEW_DECOMPOSITION%UBVEC = 1.000001_RP 

              ALLOCATE(NEW_DECOMPOSITION%ELEMENT_DOMAIN(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decomposition element domain.",ERR,ERROR,*999)

              ALLOCATE(NEW_DECOMPOSITION%NODE_DOMAIN(0:MESH%TOPOLOGY(1)%PTR%NODES%numberOfnodes-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decomposition node domain.",ERR,ERROR,*999)

              NEW_DECOMPOSITION%ELEMENT_DOMAIN=0          
              !Nullify the domain
              NULLIFY(NEW_DECOMPOSITION%DOMAIN)
              !Nullify the topology
              NULLIFY(NEW_DECOMPOSITION%TOPOLOGY)
              !\todo change this to use move alloc.
              !Add new decomposition into list of decompositions on the mesh
              ALLOCATE(NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
              DO decomposition_no=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
                NEW_DECOMPOSITIONS(decomposition_no)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR
              ENDDO !decomposition_no
              NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1)%PTR=>NEW_DECOMPOSITION
              IF(ASSOCIATED(MESH%DECOMPOSITIONS%DECOMPOSITIONS)) DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
              MESH%DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
              MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1        
              DECOMPOSITION=>NEW_DECOMPOSITION
            ENDIF
          ELSE
            LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
              & " are not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITION)) THEN
      IF(ALLOCATED(NEW_DECOMPOSITION%ELEMENT_DOMAIN)) DEALLOCATE(NEW_DECOMPOSITION%ELEMENT_DOMAIN)
      DEALLOCATE(NEW_DECOMPOSITION)
    ENDIF
    IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    NULLIFY(DECOMPOSITION)
    ERRORSEXITS("DECOMPOSITION_CREATE_START",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a user number and deallocates all memory. \see OPENCMISS::CMISSDecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY_NUMBER(USER_NUMBER,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to destroy.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh containing the decomposition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITIONS)

    ENTERS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN

        !Find the decomposition identified by the user number
        FOUND=.FALSE.
        decomposition_position=0
        DO WHILE(decomposition_position<MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.FOUND)
          decomposition_position=decomposition_position+1
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          DECOMPOSITION=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR

          !Destroy all the decomposition components          
          IF(ALLOCATED(DECOMPOSITION%ELEMENT_DOMAIN)) DEALLOCATE(DECOMPOSITION%ELEMENT_DOMAIN)
          CALL DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
          CALL DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
          
          DEALLOCATE(DECOMPOSITION)

          !Remove the decomposition from the list of decompositions
          IF(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
            ALLOCATE(NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
            DO decomposition_idx=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
              IF(decomposition_idx<decomposition_position) THEN
                NEW_DECOMPOSITIONS(decomposition_idx)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ELSE IF(decomposition_idx>decomposition_position) THEN
                MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                  & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
                NEW_DECOMPOSITIONS(decomposition_idx-1)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ENDIF
            ENDDO !decomposition_idx
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
          ELSE
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    ERRORSEXITS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a pointer and deallocates all memory. \see OPENCMISS::CMISSDecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: DECOMPOSITIONS
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITIONS)

    ENTERS("DECOMPOSITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      DECOMPOSITIONS=>DECOMPOSITION%DECOMPOSITIONS
      IF(ASSOCIATED(DECOMPOSITIONS)) THEN
        decomposition_position=DECOMPOSITION%GLOBAL_NUMBER

        !Destroy all the decomposition components          
        IF(ALLOCATED(DECOMPOSITION%ELEMENT_DOMAIN)) DEALLOCATE(DECOMPOSITION%ELEMENT_DOMAIN)
        CALL DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
        CALL DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
        
        DEALLOCATE(DECOMPOSITION)
        
        !Remove the decomposition from the list of decompositions
        IF(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
          ALLOCATE(NEW_DECOMPOSITIONS(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
          DO decomposition_idx=1,DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
            IF(decomposition_idx<decomposition_position) THEN
              NEW_DECOMPOSITIONS(decomposition_idx)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
            ELSE IF(decomposition_idx>decomposition_position) THEN
              DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                & DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
              NEW_DECOMPOSITIONS(decomposition_idx-1)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
            ENDIF
          ENDDO !decomposition_idx
          DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
          DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
          DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
        ELSE
          DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
          DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        ENDIF
      ELSE
        CALL FlagError("Decomposition decompositions is not associated.",ERR,ERROR,*999)
       ENDIF
    ELSE
      CALL FlagError("Decompositions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    ERRORSEXITS("DECOMPOSITION_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition of a mesh. \see OPENCMISS::CMISSDecompositionElementDomainCalculate
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_elem_indicies,elem_index,elem_count,ne,nn,MY_COMPUTATIONAL_NODE_NUMBER,number_computational_nodes, &
      & no_computational_node,ELEMENT_START,ELEMENT_STOP,MY_ELEMENT_START,MY_ELEMENT_STOP,NUMBER_OF_ELEMENTS, &
      & MY_NUMBER_OF_ELEMENTS,MPI_IERROR,MAX_NUMBER_ELEMENTS_PER_NODE,component_idx,minNumberXi
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_COUNT(:),ELEMENT_PTR(:),ELEMENT_INDICIES(:),ELEMENT_DISTANCE(:),DISPLACEMENTS(:), &
      & RECEIVE_COUNTS(:)
    INTEGER(INTG) :: ELEMENT_WEIGHT(1),WEIGHT_FLAG,NUMBER_FLAG,NUMBER_OF_CONSTRAINTSS, &
      & NUMBER_OF_COMMON_NODES,PARMETIS_OPTIONS(0:2)
    !ParMETIS now has double for these
    !REAL(SP) :: UBVEC(1)
    !REAL(SP), ALLOCATABLE :: TPWGTS(:)
    REAL(DP) :: UBVEC(1)
    REAL(DP), ALLOCATABLE :: TPWGTS(:)
    REAL(DP) :: NUMBER_ELEMENTS_PER_NODE
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN

          component_idx=DECOMPOSITION%MESH_COMPONENT_NUMBER
          
          number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
          IF(ERR/=0) GOTO 999
          MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
          IF(ERR/=0) GOTO 999
          
          SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
          CASE(DECOMPOSITION_ALL_TYPE)
            !Do nothing. Decomposition checked below.
           CASE(DECOMPOSITION_CALCULATED_TYPE)
            !Calculate the general decomposition

            IF(DECOMPOSITION%NUMBER_OF_DOMAINS==1) THEN
              DECOMPOSITION%ELEMENT_DOMAIN=0
            ELSE
              number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              
              NUMBER_ELEMENTS_PER_NODE=REAL(MESH%NUMBER_OF_ELEMENTS,DP)/REAL(number_computational_nodes,DP)
              ELEMENT_START=1
              ELEMENT_STOP=0
              MAX_NUMBER_ELEMENTS_PER_NODE=-1
              ALLOCATE(RECEIVE_COUNTS(0:number_computational_nodes-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate recieve counts.",ERR,ERROR,*999)
              ALLOCATE(DISPLACEMENTS(0:number_computational_nodes-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate displacements.",ERR,ERROR,*999)
              ALLOCATE(ELEMENT_DISTANCE(0:number_computational_nodes),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element distance.",ERR,ERROR,*999)
              ELEMENT_DISTANCE(0)=0
              DO no_computational_node=0,number_computational_nodes-1
                ELEMENT_START=ELEMENT_STOP+1
                IF(no_computational_node==number_computational_nodes-1) THEN
                  ELEMENT_STOP=MESH%NUMBER_OF_ELEMENTS
                ELSE
                  ELEMENT_STOP=ELEMENT_START+NINT(NUMBER_ELEMENTS_PER_NODE,INTG)-1
                ENDIF
                IF((number_computational_nodes-1-no_computational_node)>(MESH%NUMBER_OF_ELEMENTS-ELEMENT_STOP)) &
                  & ELEMENT_STOP=MESH%NUMBER_OF_ELEMENTS-(number_computational_nodes-1-no_computational_node)
                IF(ELEMENT_START>MESH%NUMBER_OF_ELEMENTS) ELEMENT_START=MESH%NUMBER_OF_ELEMENTS
                IF(ELEMENT_STOP>MESH%NUMBER_OF_ELEMENTS) ELEMENT_STOP=MESH%NUMBER_OF_ELEMENTS
                DISPLACEMENTS(no_computational_node)=ELEMENT_START-1
                ELEMENT_DISTANCE(no_computational_node+1)=ELEMENT_STOP !C numbering
                NUMBER_OF_ELEMENTS=ELEMENT_STOP-ELEMENT_START+1
                RECEIVE_COUNTS(no_computational_node)=NUMBER_OF_ELEMENTS
                IF(NUMBER_OF_ELEMENTS>MAX_NUMBER_ELEMENTS_PER_NODE) MAX_NUMBER_ELEMENTS_PER_NODE=NUMBER_OF_ELEMENTS
                IF(no_computational_node==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                  MY_ELEMENT_START=ELEMENT_START
                  MY_ELEMENT_STOP=ELEMENT_STOP
                  MY_NUMBER_OF_ELEMENTS=ELEMENT_STOP-ELEMENT_START+1
                  number_elem_indicies=0
                  DO ne=MY_ELEMENT_START,MY_ELEMENT_STOP
                    BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                    number_elem_indicies=number_elem_indicies+BASIS%NUMBER_OF_NODES
                  ENDDO !ne
                ENDIF
              ENDDO !no_computational_node
              
              ALLOCATE(ELEMENT_PTR(0:MY_NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element pointer list.",ERR,ERROR,*999)
              ALLOCATE(ELEMENT_INDICIES(0:number_elem_indicies-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element indicies list.",ERR,ERROR,*999)
              ALLOCATE(TPWGTS(1:DECOMPOSITION%NUMBER_OF_DOMAINS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate tpwgts.",ERR,ERROR,*999)
              elem_index=0
              elem_count=0
              ELEMENT_PTR(0)=0
              minNumberXi=99999
              DO ne=MY_ELEMENT_START,MY_ELEMENT_STOP
                elem_count=elem_count+1
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                IF(BASIS%NUMBER_OF_XI<minNumberXi) minNumberXi=BASIS%NUMBER_OF_XI
                DO nn=1,BASIS%NUMBER_OF_NODES
                  ELEMENT_INDICIES(elem_index)=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)% &
                    & MESH_ELEMENT_NODES(nn)-1 !C numbering
                  elem_index=elem_index+1
                ENDDO !nn
                ELEMENT_PTR(elem_count)=elem_index !C numbering
              ENDDO !ne
              
              !Set up ParMETIS variables
              WEIGHT_FLAG=0 !No weights
              ELEMENT_WEIGHT(1)=1 !Isn't used due to weight flag
              NUMBER_FLAG=0 !C Numbering as there is a bug with Fortran numbering
              NUMBER_OF_CONSTRAINTSS=1
              IF(minNumberXi==1) THEN
                NUMBER_OF_COMMON_NODES=1
              ELSE
                NUMBER_OF_COMMON_NODES=2
              ENDIF
              !ParMETIS now has doule precision for these
              !TPWGTS=1.0_SP/REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,SP)
              !UBVEC=1.05_SP
              TPWGTS=1.0_DP/REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,DP)
              UBVEC=1.05_DP
              PARMETIS_OPTIONS(0)=1 !If zero, defaults are used, otherwise next two values are used
              PARMETIS_OPTIONS(1)=7 !Level of information to output
              PARMETIS_OPTIONS(2)=CMISS_RANDOM_SEEDS(1) !Seed for random number generator
              
              !Call ParMETIS to calculate the partitioning of the mesh graph.
              CALL PARMETIS_PARTMESHKWAY(ELEMENT_DISTANCE,ELEMENT_PTR,ELEMENT_INDICIES,ELEMENT_WEIGHT,WEIGHT_FLAG,NUMBER_FLAG, &
                & NUMBER_OF_CONSTRAINTSS,NUMBER_OF_COMMON_NODES,DECOMPOSITION%NUMBER_OF_DOMAINS,TPWGTS,UBVEC,PARMETIS_OPTIONS, &
                & DECOMPOSITION%NUMBER_OF_EDGES_CUT,DECOMPOSITION%ELEMENT_DOMAIN(DISPLACEMENTS(MY_COMPUTATIONAL_NODE_NUMBER)+1:), &
                & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,ERR,ERROR,*999)
              
              !Transfer all the element domain information to the other computational nodes so that each rank has all the info
              IF(number_computational_nodes>1) THEN
                !This should work on a single processor but doesn't for mpich2 under windows. Maybe a bug? Avoid for now.
                CALL MPI_ALLGATHERV(MPI_IN_PLACE,MAX_NUMBER_ELEMENTS_PER_NODE,MPI_INTEGER,DECOMPOSITION%ELEMENT_DOMAIN, &
                  & RECEIVE_COUNTS,DISPLACEMENTS,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
              ENDIF
              
              DEALLOCATE(DISPLACEMENTS)
              DEALLOCATE(RECEIVE_COUNTS)
              DEALLOCATE(ELEMENT_DISTANCE)
              DEALLOCATE(ELEMENT_PTR)
              DEALLOCATE(ELEMENT_INDICIES)
              DEALLOCATE(TPWGTS)

            ENDIF
            
          CASE(DECOMPOSITION_USER_DEFINED_TYPE)
            !Do nothing. Decomposition checked below.          
          CASE DEFAULT
            CALL FlagError("Invalid domain decomposition type.",ERR,ERROR,*999)            
          END SELECT

          !Check decomposition and check that each domain has an element in it.
          ALLOCATE(ELEMENT_COUNT(0:number_computational_nodes-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate element count.",ERR,ERROR,*999)
          ELEMENT_COUNT=0
          DO elem_index=1,MESH%NUMBER_OF_ELEMENTS
            no_computational_node=DECOMPOSITION%ELEMENT_DOMAIN(elem_index)
            IF(no_computational_node>=0.AND.no_computational_node<number_computational_nodes) THEN
              ELEMENT_COUNT(no_computational_node)=ELEMENT_COUNT(no_computational_node)+1
            ELSE
              LOCAL_ERROR="The computational node number of "//TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))// &
                & " for element number "//TRIM(NUMBER_TO_VSTRING(elem_index,"*",ERR,ERROR))// &
                & " is invalid. The computational node number must be between 0 and "// &
                & TRIM(NUMBER_TO_VSTRING(number_computational_nodes-1,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !elem_index
          DO no_computational_node=0,number_computational_nodes-1
            IF(ELEMENT_COUNT(no_computational_node)==0) THEN
              LOCAL_ERROR="Invalid decomposition. There are no elements in computational node "// &
                & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !no_computational_node
          DEALLOCATE(ELEMENT_COUNT)

          
        ELSE
          CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition for mesh number ",DECOMPOSITION%MESH%USER_NUMBER, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains = ", DECOMPOSITION%NUMBER_OF_DOMAINS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Element domains:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposition type = ", DECOMPOSITION%DECOMPOSITION_TYPE, &
        & ERR,ERROR,*999)
      IF(DECOMPOSITION%DECOMPOSITION_TYPE==DECOMPOSITION_CALCULATED_TYPE) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of edges cut = ",DECOMPOSITION%NUMBER_OF_EDGES_CUT, &
          & ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of elements = ",DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS, &
        & ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Domain = ",DECOMPOSITION%ELEMENT_DOMAIN(ne), &
          & ERR,ERROR,*999)
      ENDDO !ne
    ENDIF
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE")
    RETURN
999 IF(ALLOCATED(RECEIVE_COUNTS)) DEALLOCATE(RECEIVE_COUNTS)
    IF(ALLOCATED(DISPLACEMENTS)) DEALLOCATE(DISPLACEMENTS)
    IF(ALLOCATED(ELEMENT_DISTANCE)) DEALLOCATE(ELEMENT_DISTANCE)
    IF(ALLOCATED(ELEMENT_PTR)) DEALLOCATE(ELEMENT_PTR)
    IF(ALLOCATED(ELEMENT_INDICIES)) DEALLOCATE(ELEMENT_INDICIES)
    IF(ALLOCATED(TPWGTS)) DEALLOCATE(TPWGTS)
    ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Gets the domain for a given element in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OPENCMISS::CMISSDecompositionElementDomainGet
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET(DECOMPOSITION,USER_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables`
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS


    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR,*999)

    GLOBAL_ELEMENT_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
            IF(ASSOCIATED(MESH_ELEMENTS)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_ELEMENTS%ELEMENTS_TREE,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_ELEMENTS%ELEMENTS_TREE,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
                  DOMAIN_NUMBER=DECOMPOSITION%ELEMENT_DOMAIN(GLOBAL_ELEMENT_NUMBER)
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Decomposition mesh element corresponding to user number not found.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Decomposition mesh elements are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET
  
  !
  !================================================================================================================================
  !

  !>Sets the domain for a given element in a decomposition of a mesh. \todo move to user number, should be able to specify lists of elements. \see OPENCMISS::CMISSDecompositionElementDomainSet 
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET(DECOMPOSITION,GLOBAL_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_ELEMENT_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: DOMAIN_NUMBER !<The domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_computational_nodes
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
              number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              IF(DOMAIN_NUMBER>=0.AND.DOMAIN_NUMBER<number_computational_nodes) THEN
                DECOMPOSITION%ELEMENT_DOMAIN(GLOBAL_ELEMENT_NUMBER)=DOMAIN_NUMBER
              ELSE
                LOCAL_ERROR="Domain number "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid. The limits are 0 to "//TRIM(NUMBER_TO_VSTRING(number_computational_nodes,"*",ERR,ERROR))//"."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The limits are 1 to "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::CMISSDecompositionMeshComponentGet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the mesh component for
    INTEGER(INTG), INTENT(OUT) :: MESH_COMPONENT_NUMBER !<On return, the mesh component number to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET
  
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::CMISSDecompositionMeshComponentSet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS) THEN
            DECOMPOSITION%MESH_COMPONENT_NUMBER=MESH_COMPONENT_NUMBER
          ELSE
            LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
              & "is invalid. The component number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the number of domains for a decomposition. \see OPENCMISS::CMISSDecompositionNumberOfDomainsGet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the number of domains for.
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_DOMAINS !<On return, the number of domains to get.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        NUMBER_OF_DOMAINS=DECOMPOSITION%NUMBER_OF_DOMAINS
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition. \see OPENCMISS::CMISSDecompositionNumberOfDomainsSet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the number of domains for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !<The number of domains to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          IF(NUMBER_OF_DOMAINS==1) THEN
            DECOMPOSITION%NUMBER_OF_DOMAINS=1
          ELSE
            CALL FlagError("Can only have one domain for all decomposition type.",ERR,ERROR,*999)
          ENDIF
        CASE(DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE)
          IF(NUMBER_OF_DOMAINS>=1) THEN
            !wolfye???<=?
            IF(NUMBER_OF_DOMAINS<=DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS) THEN
              !Get the number of computational nodes
              NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              !!TODO: relax this later
              !IF(NUMBER_OF_DOMAINS==NUMBER_COMPUTATIONAL_NODES) THEN
                DECOMPOSITION%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS             
              !ELSE
              !  LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINSS,"*",ERR,ERROR))// &
              !    & ") is not equal to the number of computational nodes ("// &
              !    & TRIM(NUMBER_TO_VSTRING(NUMBER_COMPUTATIONAL_NODES,"*",ERR,ERROR))//")"
              !  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              !ENDIF
            ELSE
              LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINS,"*",ERR,ERROR))// &
                & ") must be <= the number of global elements ("// &
                & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//") in the mesh."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
             CALL FlagError("Number of domains must be >= 1.",ERR,ERROR,*999)
           ENDIF
         CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%DECOMPOSITION_TYPE,"*",ERR,ERROR))// &
            & " is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  !
  !================================================================================================================================
  !

  !>Calculates the topology for a decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
      !Calculate the elements topology
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      !Calculate the line topology
      IF(DECOMPOSITION%CALCULATE_LINES)THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      !Calculate the face topology
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      meshComponentNumber=DECOMPOSITION%MESH_COMPONENT_NUMBER
      IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsCalculate(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        END IF  
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DecompositionTopology_DataPointsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: localElement,globalElement,dataPointIdx,localData,meshComponentNumber
    INTEGER(INTG) :: INSERT_STATUS,MPI_IERROR,NUMBER_OF_COMPUTATIONAL_NODES,MY_COMPUTATIONAL_NODE_NUMBER,NUMBER_OF_GHOST_DATA, &
      & NUMBER_OF_LOCAL_DATA
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(MeshDataPointsType), POINTER :: meshData

    ENTERS("DecompositionTopology_DataPointsCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      decompositionData=>TOPOLOGY%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        decomposition=>decompositionData%DECOMPOSITION
        IF(ASSOCIATED(decomposition)) THEN
         decompositionElements=>TOPOLOGY%ELEMENTS
         IF(ASSOCIATED(decompositionElements)) THEN
           elementsMapping=>decomposition%DOMAIN(decomposition%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
           IF(ASSOCIATED(elementsMapping)) THEN
              meshComponentNumber=decomposition%MESH_COMPONENT_NUMBER
              meshData=>decomposition%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints
              IF(ASSOCIATED(meshData)) THEN
                NUMBER_OF_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
                IF(ERR/=0) GOTO 999
                MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                IF(ERR/=0) GOTO 999
                ALLOCATE(decompositionData%numberOfDomainLocal(0:NUMBER_OF_COMPUTATIONAL_NODES-1),STAT=ERR)
                ALLOCATE(decompositionData%numberOfDomainGhost(0:NUMBER_OF_COMPUTATIONAL_NODES-1),STAT=ERR)
                ALLOCATE(decompositionData%numberOfElementDataPoints(decompositionElements%NUMBER_OF_GLOBAL_ELEMENTS),STAT=ERR)
                ALLOCATE(decompositionData%elementDataPoint(decompositionElements%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate decomposition element data points.",ERR,ERROR,*999)
                CALL TREE_CREATE_START(decompositionData%dataPointsTree,ERR,ERROR,*999)
                CALL TREE_INSERT_TYPE_SET(decompositionData%dataPointsTree,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                CALL TREE_CREATE_FINISH(decompositionData%dataPointsTree,ERR,ERROR,*999)
                decompositionData%numberOfGlobalDataPoints=meshData%totalNumberOfProjectedData 
                DO globalElement=1,decompositionElements%NUMBER_OF_GLOBAL_ELEMENTS
                  decompositionData%numberOfElementDataPoints(globalElement)= &
                    & meshData%elementDataPoint(globalElement)%numberOfProjectedData
                ENDDO !globalElement
                localData=0;
                DO localElement=1,decompositionElements%TOTAL_NUMBER_OF_ELEMENTS
                  globalElement=decompositionElements%ELEMENTS(localElement)%GLOBAL_NUMBER
                  decompositionData%elementDataPoint(localElement)%numberOfProjectedData= &
                    & meshData%elementDataPoint(globalElement)%numberOfProjectedData
                  decompositionData%elementDataPoint(localElement)%globalElementNumber=globalElement
                  IF(localElement<elementsMapping%GHOST_START) THEN
                    decompositionData%numberOfDataPoints=decompositionData%numberOfDataPoints+ &
                      & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  END IF              
                  decompositionData%totalNumberOfDataPoints=decompositionData%totalNumberOfDataPoints+ &
                    & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  ALLOCATE(decompositionData%elementDataPoint(localElement)%dataIndices(decompositionData% &
                    & elementDataPoint(localElement)%numberOfProjectedData),STAT=ERR)
                  DO dataPointIdx=1,decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%userNumber
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%globalNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%globalNumber
                    localData=localData+1
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%localNumber=localData
                    CALL TREE_ITEM_INSERT(decompositionData%dataPointsTree,decompositionData% &
                      & elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber,localData, &
                      & INSERT_STATUS,ERR,ERROR,*999)
                  ENDDO !dataPointIdx
                ENDDO !localElement   
                !Calculate number of ghost data points on the current computational domain
                NUMBER_OF_LOCAL_DATA=decompositionData%numberOfDataPoints
                NUMBER_OF_GHOST_DATA=decompositionData%totalNumberOfDataPoints-decompositionData%numberOfDataPoints
                !Gather number of local data points on all computational nodes
                CALL MPI_ALLGATHER(NUMBER_OF_LOCAL_DATA,1,MPI_INTEGER,decompositionData% &
                  & numberOfDomainLocal,1,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                !Gather number of ghost data points on all computational nodes
                CALL MPI_ALLGATHER(NUMBER_OF_GHOST_DATA,1,MPI_INTEGER,decompositionData% &
                  & numberOfDomainGhost,1,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Mesh data points topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Element mapping  is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition elements topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition data points topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsCalculate")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsCalculate",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataPointsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology for a data projection (for data projections on fields).
  SUBROUTINE DecompositionTopology_DataProjectionCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DecompositionTopology_DataProjectionCalculate",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      CALL DecompositionTopology_DataPointsInitialise(decompositionTopology,err,error,*999)
      CALL DecompositionTopology_DataPointsCalculate(decompositionTopology,err,error,*999)
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN
999 ERRORS("DecompositionTopology_DataProjectionCalculate",err,error)
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataProjectionCalculate

  !
  !================================================================================================================================
  !

  !>Gets the local data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet(decompositionTopology,elementNumber,dataPointIndex, &
       & dataPointLocalNumber,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointLocalNumber !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        numberOfDataPoints = decompositionData%elementDataPoint(elementNumber)%numberOfProjectedData
        IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
          dataPointLocalNumber = decompositionData%elementDataPoint(elementNumber)%dataIndices(dataPointIndex)%localNumber
        ELSE
          localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
           & " out of range for element "//TRIM(NUMBER_TO_VSTRING(elementNumber,"*",ERR,ERROR))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the user (global) data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet(decompositionTopology,userElementNumber,dataPointIndex, &
       & dataPointUserNumber,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointUserNumber !<The data point user number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints,decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointUserNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints = decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
            IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
              dataPointUserNumber = decompositionData%elementDataPoint(decompositionLocalElementNumber)% &
                & dataIndices(dataPointIndex)%userNumber
            ELSE
              localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
               & " out of range for element "//TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",ERR,ERROR))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointUserNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of data points projected on an element
  SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet(decompositionTopology,userElementNumber, &
       & numberOfDataPoints,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(OUT) :: numberOfDataPoints !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_NumberOfElementDataPointsGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints = decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN
999 ERRORS("DecompositionTopology_NumberOfElementDataPointsGet",err,error)
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DecompositionTopology_DataPointCheckExists(decompositionTopology,userDataPointNumber,userDataPointExists, &
        & decompositionLocalDataPointNumber,ghostDataPoint,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to check the data point exists on
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: userDataPointExists !<On exit, is .TRUE. if the data point user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: decompositionLocalDataPointNumber !<On exit, if the data point exists the local number corresponding to the user data point number. If the data point does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostDataPoint !<On exit, is .TRUE. if the local data point (if it exists) is a ghost data point, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("DecompositionTopology_DataPointCheckExists",ERR,error,*999)

    userDataPointExists=.FALSE.
    decompositionLocalDataPointNumber=0
    ghostDataPoint=.FALSE.
    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        NULLIFY(treeNode)
        CALL TREE_SEARCH(decompositionData%dataPointsTree,userDataPointNumber,treeNode,err,error,*999)
        IF(ASSOCIATED(treeNode)) THEN
          CALL TREE_NODE_VALUE_GET(decompositionData%dataPointsTree,treeNode,decompositionLocalDataPointNumber,err,error,*999)
          userDataPointExists=.TRUE.
          ghostDataPoint=decompositionLocalDataPointNumber>decompositionData%numberOfDataPoints
        ENDIF
      ELSE
        CALL FlagError("Decomposition data point topology is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN
999 ERRORS("DecompositionTopology_DataPointCheckExists",err,error)
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointCheckExists
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(DECOMPOSITION_TOPOLOGY,USER_ELEMENT_NUMBER,ELEMENT_EXISTS, &
    & DECOMPOSITION_LOCAL_ELEMENT_NUMBER,GHOST_ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: DECOMPOSITION_TOPOLOGY !<A pointer to the decomposition topology to check the element exists on
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: ELEMENT_EXISTS !<On exit, is .TRUE. if the element user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DECOMPOSITION_LOCAL_ELEMENT_NUMBER !<On exit, if the element exists the local number corresponding to the user element number. If the element does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_ELEMENT !<On exit, is .TRUE. if the local element (if it exists) is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR,*999)

    ELEMENT_EXISTS=.FALSE.
    DECOMPOSITION_LOCAL_ELEMENT_NUMBER=0
    GHOST_ELEMENT=.FALSE.
    IF(ASSOCIATED(DECOMPOSITION_TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>DECOMPOSITION_TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,TREE_NODE,DECOMPOSITION_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
          ELEMENT_EXISTS=.TRUE.
          GHOST_ELEMENT=DECOMPOSITION_LOCAL_ELEMENT_NUMBER>DECOMPOSITION_ELEMENTS%NUMBER_OF_ELEMENTS
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS

  !
  !================================================================================================================================
  !

  !>Get the basis for an element in the domain identified by its local number
  SUBROUTINE DomainTopology_ElementBasisGet(domainTopology,userElementNumber, &
      & basis,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the element basis for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element user number to get the basis for
    TYPE(BASIS_TYPE), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    LOGICAL :: userElementExists,ghostElement
    INTEGER(INTG) :: localElementNumber

    ENTERS("DomainTopology_ElementBasisGet",err,error,*999)

    NULLIFY(basis)

    IF(ASSOCIATED(domainTopology)) THEN
      domainElements=>domainTopology%elements
      IF(ASSOCIATED(domainElements)) THEN
        decompositionTopology=>domainTopology%domain%decomposition%topology
        IF(ASSOCIATED(decompositionTopology)) THEN
          CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
            & userElementExists,localElementNumber,ghostElement,err,error,*999)
          IF(.NOT.userElementExists) THEN
            CALL FlagError("The specified user element number of "// &
              & TRIM(NumberToVstring(userElementNumber,"*",err,error))// &
              & " does not exist in the domain decomposition.",err,error,*999)
          END IF
          basis=>domainElements%elements(localElementNumber)%basis
        ELSE
          CALL FlagError("Decomposition topology is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Domain topology elements is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Domain topology is not associated.",err,error,*999)
    END IF

    EXITS("DomainTopology_ElementBasisGet")
    RETURN
    999 ERRORSEXITS("DomainTopology_ElementBasisGet",err,error)
    RETURN 1

    END SUBROUTINE DomainTopology_ElementBasisGet

  !
  !================================================================================================================================
  !

  !>Finalises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic !\todo add comment

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) THEN
      DO nic=LBOUND(ELEMENT%ADJACENT_ELEMENTS,1),UBOUND(ELEMENT%ADJACENT_ELEMENTS,1)
        CALL DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(ELEMENT%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    ENDIF
    IF(ALLOCATED(ELEMENT%ELEMENT_LINES)) DEALLOCATE(ELEMENT%ELEMENT_LINES)
    IF(ALLOCATED(ELEMENT%ELEMENT_FACES)) DEALLOCATE(ELEMENT%ELEMENT_FACES)
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%LOCAL_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    ELEMENT%BOUNDARY_ELEMENT=.FALSE.
  
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers adjacent to an element in a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the adjacent elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),NODE_POSITION_INDEX(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: candidate_idx,face_node_idx,node_idx,surrounding_el_idx,candidate_el,idx
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4),numberSurroundingElements
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:), surroundingElements(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST, surroundingElementsList
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    ENTERS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                  !Loop over the elements in the decomposition
                  DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                    !Create a list for every xi direction (plus and minus)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
                      CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                    ENDDO !nic
                    NUMBER_OF_NODES_XIC=1
                    NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)= &
                      & BASIS%NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)
                    !Place the current element in the surrounding list
                    CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER, &
                      & ERR,ERROR,*999)
                    SELECT CASE(BASIS%TYPE)
                    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
!!TODO: Calculate this and set it as part of the basis type
                      !Determine the collapsed "faces" if any
                      NODE_POSITION_INDEX=1
                      !Loop over the face normals of the element
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the face xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Reset the node_position_index in this xi direction
                        NODE_POSITION_INDEX(ni)=1
                        !Loop over the two faces with this normal
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni
                          FACE_COLLAPSED(xi_direction)=.FALSE.
                          DO j=1,2
                            xi_dir_check=FACE_XI(j)
                            IF(xi_dir_check<=BASIS%NUMBER_OF_XI) THEN
                              xi_dir_search=FACE_XI(3-j)
                              NODE_POSITION_INDEX(xi_dir_search)=1
                              XI_COLLAPSED=.TRUE.
                              DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                                !Get the first local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=1
                                nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                !Get the second local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=2
                                nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                IF(nn1/=0.AND.nn2/=0) THEN
                                  IF(DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn1)/= &
                                    & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn2)) XI_COLLAPSED=.FALSE.
                                ENDIF
                                NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                              ENDDO !xi_dir_search
                              IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                            ENDIF
                          ENDDO !j
                          NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                        ENDDO !direction_index
                      ENDDO !ni
                      !Loop over the xi directions and calculate the surrounding elements
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Loop over the two faces
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni                  
                          !Find nodes in the element on the appropriate face/line/point
                          NULLIFY(NODE_MATCH_LIST)
                          CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                          IF(direction_index==-1) THEN
                            NODE_POSITION_INDEX(ni)=1
                          ELSE
                            NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                          ENDIF
                          !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is
                          !also collapsed. This may indicate that we have a funny element in non-rc coordinates that goes around the
                          !central axis back to itself
                          IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                            !Do nothing - the match lists are already empty
                          ELSE
                            !Find the nodes to match and add them to the node match list
                            SELECT CASE(BASIS%NUMBER_OF_XI)
                            CASE(1)
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),1,1,1)
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            CASE(2)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),1,1)
                                IF(nn/=0) THEN
                                  np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                  CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                ENDIF
                              ENDDO !nn1
                            CASE(3)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2)),NUMBER_OF_NODES_XIC(FACE_XI(2))-1
                                  NODE_POSITION_INDEX(FACE_XI(2))=nn2
                                  nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                    & NODE_POSITION_INDEX(3),1)
                                  IF(nn/=0) THEN
                                    np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                    CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !nn2
                              ENDDO !nn1
                            CASE DEFAULT
                              LOCAL_ERROR="The number of xi directions in the basis of "// &
                                & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//" is invalid."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ENDIF
                          CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=0
                          IF(NUMBER_NODE_MATCHES>0) THEN
                            !NODE_MATCHES now contain the list of corner nodes in the current face with normal_xi of ni.
                            !Look at the surrounding elements of each of these nodes, if there is a repeated element that
                            !is not the current element ne, it's an adjacent element.
                            candidate_idx=0
                            NULLIFY(surroundingElementsList)
                            CALL LIST_CREATE_START(surroundingElementsList,ERR,ERROR,*999)
                            CALL LIST_DATA_TYPE_SET(surroundingElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                            CALL LIST_INITIAL_SIZE_SET(surroundingElementsList,2,ERR,ERROR,*999)
                            CALL LIST_CREATE_FINISH(surroundingElementsList,ERR,ERROR,*999)
                            DO face_node_idx=1,NUMBER_NODE_MATCHES
                              !Dump all the surrounding elements into an array, see if any are repeated
                              node_idx=NODE_MATCHES(face_node_idx)
                              DO surrounding_el_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                                candidate_el=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(surrounding_el_idx)
                                IF(candidate_el/=ne) THEN
                                  candidate_idx=candidate_idx+1
                                  CALL LIST_ITEM_ADD(surroundingElementsList,candidate_el,ERR,ERROR,*999)
                                ENDIF
                              ENDDO
                            ENDDO !face_node_idx
                            CALL LIST_DETACH_AND_DESTROY(surroundingElementsList,numberSurroundingElements,surroundingElements, &
                              & ERR,ERROR,*999)
                            DO idx=1,candidate_idx
                              ne1=surroundingElements(idx)
                              IF(COUNT(surroundingElements(1:numberSurroundingElements)==ne1)>=BASIS%NUMBER_OF_XI) THEN
                                !Found it, just exit
                                CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                                NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                                EXIT
                              ENDIF
                            ENDDO
                          ENDIF
                          IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                          IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
                        ENDDO !direction_index
                      ENDDO !ni
                    CASE(BASIS_SIMPLEX_TYPE)
                      !Loop over the xi coordinates and calculate the surrounding elements
                      DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                        !Find the other coordinates of the face/line/point
                        FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                        FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                        FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                        !Find nodes in the element on the appropriate face/line/point
                        NULLIFY(NODE_MATCH_LIST)
                        CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                        CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                        CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                        NODE_POSITION_INDEX(nic)=1 !Furtherest away from node with the nic'th coordinate
                        !Find the nodes to match and add them to the node match list
                        DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                          NODE_POSITION_INDEX(FACE_XIC(1))=nn1
                          DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                            NODE_POSITION_INDEX(FACE_XIC(2))=nn2
                            DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                              NODE_POSITION_INDEX(FACE_XIC(3))=nn3
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                & NODE_POSITION_INDEX(3),NODE_POSITION_INDEX(4))
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !nn3
                          ENDDO !nn2
                        ENDDO !nn1
                        CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                        IF(NUMBER_NODE_MATCHES>0) THEN
                          !Find list of elements surrounding those nodes
                          DO node_idx=1,NUMBER_NODE_MATCHES
                            np1=NODE_MATCHES(node_idx)
                            DO nep1=1,DOMAIN_NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                              ne1=DOMAIN_NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                              IF(ne1/=ne) THEN !Don't want the current element
                                ! grab the nodes list for current and this surrouding elements
                                ! current face : NODE_MATCHES
                                ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES 
                                ! if all of current face belongs to the candidate element, we will have found the neighbour
                                CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),DOMAIN_ELEMENTS%ELEMENTS(ne1)% &
                                  & ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                                IF(SUBSET) THEN
                                  CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                            ENDDO !nep1
                          ENDDO !node_idx
                        ENDIF
                        IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                      ENDDO !nic
                    CASE(BASIS_SERENDIPITY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_AUXILLIARY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_B_SPLINE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)                     
                    END SELECT
                    !Set the surrounding elements for this element
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
                      BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements.",ERR,ERROR,*999)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      CALL DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic), &
                        & ERR,ERROR,*999)
                      CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS,ERR,ERROR,*999)
                      ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS( &
                        DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate element adjacent elements.",ERR,ERROR,*999)
                      DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(1: &
                        & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)= &
                        ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)
                      IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
                    ENDDO !nic
                  ENDDO !ne           
                ELSE
                  CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not allocated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ",DECOMPOSITION_ELEMENTS% &
        & TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number : ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ERR,ERROR,*999)
          IF(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
              & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,8,8,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)% &
              & ADJACENT_ELEMENTS,'("        Adjacent elements :",8(X,I6))','(30x,8(X,I6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF

    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(surroundingElementsList)) CALL LIST_DESTROY(surroundingElementsList,DUMMY_ERR,DUMMY_ERROR,*997)
997 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*996)
    ENDDO !ni
996 ERRORS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR)
    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: global_element,INSERT_STATUS,local_element
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS

      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        DECOMPOSITION=>TOPOLOGY%DECOMPOSITION

        IF(ASSOCIATED(DECOMPOSITION)) THEN

          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR

          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY

            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS

              IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS

                IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                  DOMAIN_ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS

                  IF(ASSOCIATED(DOMAIN_ELEMENTS_MAPPING)) THEN
                    MESH=>DECOMPOSITION%MESH

                    IF(ASSOCIATED(MESH)) THEN
                      MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
                      IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
                        MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(MESH_ELEMENTS)) THEN
                          !Allocate the element topology arrays
                          ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate decomposition elements elements.",ERR,ERROR,*999)

                          DECOMPOSITION_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS
                          CALL TREE_CREATE_START(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                          CALL TREE_INSERT_TYPE_SET(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                          CALL TREE_CREATE_FINISH(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                          DO local_element=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                            CALL DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(local_element), &
                              & ERR,ERROR,*999)
                            global_element=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_element)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%USER_NUMBER=MESH_ELEMENTS%ELEMENTS(global_element)% &
                              & USER_NUMBER
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%LOCAL_NUMBER=local_element
                            CALL TREE_ITEM_INSERT(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,DECOMPOSITION_ELEMENTS% &
                              & ELEMENTS(local_element)%USER_NUMBER,local_element,INSERT_STATUS,ERR,ERROR,*999)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%GLOBAL_NUMBER=global_element
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%BOUNDARY_ELEMENT=MESH_ELEMENTS% &
                              & ELEMENTS(global_element)%BOUNDARY_ELEMENT
                          ENDDO !local_element
                          !Calculate the elements surrounding the elements in the decomposition topology
                          CALL DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given decomposition topology. \todo Pass in the decomposition elements pointer.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)) CALL TREE_DESTROY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Decomposition already has topology elements associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements.",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given decomposition. \todo Pass in a pointer to the decomposition topology
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      IF(DECOMPOSITION%CALCULATE_LINES) THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      DEALLOCATE(DECOMPOSITION%TOPOLOGY)
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
        CALL FlagError("Decomposition already has topology associated.",ERR,ERROR,*999)
      ELSE
        !Allocate decomposition topology
        ALLOCATE(DECOMPOSITION%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Decomposition topology could not be allocated.",ERR,ERROR,*999)
        DECOMPOSITION%TOPOLOGY%DECOMPOSITION=>DECOMPOSITION
        NULLIFY(DECOMPOSITION%TOPOLOGY%ELEMENTS)
        NULLIFY(DECOMPOSITION%TOPOLOGY%LINES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%FACES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%dataPoints)
        !Initialise the topology components
        CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        IF(DECOMPOSITION%CALCULATE_LINES) THEN !Default is currently true
          CALL DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF
        IF(DECOMPOSITION%CALCULATE_FACES) THEN !Default is currently false
          CALL DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        END IF     
        meshComponentNumber=DECOMPOSITION%MESH_COMPONENT_NUMBER
        IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsInitialise(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        END IF  
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(LINE%SURROUNDING_ELEMENTS)) DEALLOCATE(LINE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(LINE%ELEMENT_LINES)) DEALLOCATE(LINE%ELEMENT_LINES)
    LINE%ADJACENT_LINES=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a decomposition topology line.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    LINE%ADJACENT_LINES=0
    LINE%BOUNDARY_LINE=.FALSE.

    EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the lines in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,element_idx,surrounding_element_idx,basis_local_line_idx, &
      & surrounding_element_basis_local_line_idx,element_local_node_idx,basis_local_line_node_idx,derivative_idx,version_idx, &
      & local_line_idx,surrounding_element_local_line_idx,node_idx,local_node_idx,elem_idx,line_end_node_idx,basis_node_idx, &
      & NODES_IN_LINE(4),NUMBER_OF_LINES,MAX_NUMBER_OF_LINES,NEW_MAX_NUMBER_OF_LINES,LINE_NUMBER,COUNT
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_LINES(:)
    INTEGER(INTG), POINTER :: TEMP_LINES(:,:),NEW_TEMP_LINES(:,:)
    REAL(DP) :: APPROX_DIMENSION
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_LINE_TYPE), POINTER :: DECOMPOSITION_LINE,DECOMPOSITION_LINE2
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: DECOMPOSITION_LINES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE,DOMAIN_LINE2
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_LINES)
    NULLIFY(NEW_TEMP_LINES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_LINES=>TOPOLOGY%LINES
      IF(ASSOCIATED(DECOMPOSITION_LINES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish line
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Guestimate the number of lines
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      MAX_NUMBER_OF_LINES=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    CASE(2)
                      APPROX_DIMENSION=SQRT(REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP))
                      !This should give the maximum and will over estimate the number of lines for a "square mesh" by approx 33%
                      MAX_NUMBER_OF_LINES=NINT(3.0_DP*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE(3)
                      !This should give the maximum and will over estimate the number of lines for a "cube mesh" by approx 73%
                      APPROX_DIMENSION=REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP)**(1.0_DP/3.0_DP)
                      MAX_NUMBER_OF_LINES=NINT(11.0_DP*APPROX_DIMENSION*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain.",ERR,ERROR,*999)
                    END SELECT
                    DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES
                    IF(ASSOCIATED(DOMAIN_LINES)) THEN
                      ALLOCATE(TEMP_LINES(4,MAX_NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate temporary lines array.",ERR,ERROR,*999)
                      ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                      NODES_NUMBER_OF_LINES=0
                      NUMBER_OF_LINES=0
                      TEMP_LINES=0
                      !Loop over the elements in the topology
                      DO element_idx=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_LINES(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate element element lines.",ERR,ERROR,*999)
                        !Loop over the local lines of the element
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          !Calculate the topology node numbers that make up the line
                          NODES_IN_LINE=0
                          DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                            NODES_IN_LINE(basis_local_line_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                              & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                          ENDDO !basis_local_line_node_idx
                          !Try and find a previously created line that matches in the adjacent elements
                          FOUND=.FALSE.
                          node_idx=NODES_IN_LINE(1)
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            IF(surrounding_element_idx/=element_idx) THEN
                              IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_LINES)) THEN
                                BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                DO surrounding_element_basis_local_line_idx=1,BASIS2%NUMBER_OF_LOCAL_LINES
                                  local_line_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)% &
                                    & ELEMENT_LINES(surrounding_element_basis_local_line_idx)
                                  IF(ALL(NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx))== &
                                    & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx),local_line_idx))) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                  ENDIF
                                ENDDO !surrounding_element_basis_local_line_idx
                                IF(FOUND) EXIT
                              ENDIF
                            ENDIF
                          ENDDO !elem_idx
                          IF(FOUND) THEN
                            !Line has already been created
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=local_line_idx
                          ELSE
                            !Line has not been created
                            IF(NUMBER_OF_LINES==MAX_NUMBER_OF_LINES) THEN
                              !We are at maximum. Reallocate the LINES array to be 20% bigger and try again.
                              NEW_MAX_NUMBER_OF_LINES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_LINES,DP),INTG)
                              ALLOCATE(NEW_TEMP_LINES(4,NEW_MAX_NUMBER_OF_LINES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate new number of lines.",ERR,ERROR,*999)
                              NEW_TEMP_LINES(:,1:NUMBER_OF_LINES)=TEMP_LINES(:,1:NUMBER_OF_LINES)
                              NEW_TEMP_LINES(:,NUMBER_OF_LINES+1:NEW_MAX_NUMBER_OF_LINES)=0
                              DEALLOCATE(TEMP_LINES)
                              TEMP_LINES=>NEW_TEMP_LINES
                              NULLIFY(NEW_TEMP_LINES)
                              MAX_NUMBER_OF_LINES=NEW_MAX_NUMBER_OF_LINES
                            ENDIF
                            NUMBER_OF_LINES=NUMBER_OF_LINES+1
                            TEMP_LINES(:,NUMBER_OF_LINES)=NODES_IN_LINE
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=NUMBER_OF_LINES
                            DO basis_local_line_node_idx=1,SIZE(NODES_IN_LINE,1)
                              IF(NODES_IN_LINE(basis_local_line_node_idx)/=0) &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))= &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))+1
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      !Allocate the line arrays and set them from the LINES and NODE_LINES arrays
                      DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                        ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate node lines array.",ERR,ERROR,*999)
                        DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                      ENDDO !node_idx
                      DEALLOCATE(NODES_NUMBER_OF_LINES)
                      ALLOCATE(DECOMPOSITION_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology lines.",ERR,ERROR,*999)
                      DECOMPOSITION_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      ALLOCATE(DOMAIN_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate domain topology lines.",ERR,ERROR,*999)
                      DOMAIN_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                        CALL DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(DECOMPOSITION_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        DO basis_local_line_node_idx=1,SIZE(TEMP_LINES,1)
                          IF(TEMP_LINES(basis_local_line_node_idx,local_line_idx)/=0) THEN
                            node_idx=TEMP_LINES(basis_local_line_node_idx,local_line_idx)
                            DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES+1
                            DOMAIN_NODES%NODES(node_idx)%NODE_LINES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES)= &
                              & local_line_idx
                          ENDIF
                        ENDDO !basis_local_line_node_idx  
                      ENDDO !local_line_idx
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DOMAIN_LINE=>DOMAIN_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          IF(.NOT.ASSOCIATED(DOMAIN_LINE%BASIS)) THEN
                            DECOMPOSITION_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%ELEMENT_NUMBER=element_idx !Needs checking
                            DECOMPOSITION_LINE%XI_DIRECTION=BASIS%LOCAL_LINE_XI_DIRECTION(basis_local_line_idx)
                            DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                            ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line nodes in line.",ERR,ERROR,*999)
                            ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                              & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line derivatives in line.",ERR,ERROR,*999)
                            DOMAIN_LINE%NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx))= &
                              & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx),LINE_NUMBER)
                            DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                              !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                              version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                & basis_local_line_node_idx,basis_local_line_idx))
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                              IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                  & BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx), &
                                  & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                  & basis_local_line_node_idx,basis_local_line_idx),BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                  & basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
                              ENDIF
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      DEALLOCATE(TEMP_LINES)
                      !Calculate adjacent lines and the surrounding elements for each line
                      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                        BASIS=>DOMAIN_LINE%BASIS
                        IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                          DECOMPOSITION_LINE%BOUNDARY_LINE=.TRUE.
                          DOMAIN_LINE%BOUNDARY_LINE=.TRUE.
                        ENDIF
                        !Allocate the elements surrounding the line
                        ALLOCATE(DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line surrounding elements.",ERR,ERROR,*999)
                        ALLOCATE(DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line element lines.",ERR,ERROR,*999)
                        DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
                        DECOMPOSITION_LINE%ADJACENT_LINES=0
                        !Loop over the nodes at each end of the line
                        DO line_end_node_idx=0,1
                          FOUND=.FALSE.
                          node_idx=DOMAIN_LINE%NODES_IN_LINE(line_end_node_idx*(BASIS%NUMBER_OF_NODES-1)+1)
                          !Loop over the elements surrounding the node.
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                            DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                            !Loop over the local lines of the element
                            DO basis_local_line_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_LINES
                              surrounding_element_local_line_idx=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                              IF(surrounding_element_local_line_idx/=local_line_idx) THEN
                                DECOMPOSITION_LINE2=>DECOMPOSITION_LINES%LINES(surrounding_element_local_line_idx)
                                DOMAIN_LINE2=>DOMAIN_LINES%LINES(surrounding_element_local_line_idx)
                                IF(DECOMPOSITION_LINE2%XI_DIRECTION==DECOMPOSITION_LINE%XI_DIRECTION) THEN
                                  !Lines run in the same direction.
                                  BASIS2=>DOMAIN_LINE2%BASIS
                                  IF(line_end_node_idx==0) THEN
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES)
                                  ELSE
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(1)
                                  ENDIF
                                  IF(local_node_idx==node_idx) THEN
                                    !The node at the 'other' end of this line matches the node at the current end of the line.
                                    !Check it is not a coexistant line running the other way
                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
                                      COUNT=0
                                      DO basis_node_idx=1,BASIS%NUMBER_OF_NODES
                                        IF(DOMAIN_LINE2%NODES_IN_LINE(basis_node_idx)== &
                                          & DOMAIN_LINE%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES-basis_node_idx+1)) &
                                          & COUNT=COUNT+1
                                      ENDDO !basis_node_idx
                                      IF(COUNT<BASIS%NUMBER_OF_NODES) THEN
                                        FOUND=.TRUE.
                                        EXIT
                                      ENDIF
                                    ELSE
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDIF
                            ENDDO !basis_local_line_idx
                            IF(FOUND) EXIT
                          ENDDO !element_idx
                          IF(FOUND) DECOMPOSITION_LINE%ADJACENT_LINES(line_end_node_idx)=surrounding_element_local_line_idx
                        ENDDO !line_end_node_idx
                      ENDDO !local_line_idx
                      !Set the surrounding elements
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=element_idx
                          DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=basis_local_line_idx
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                    ELSE
                      CALL FlagError("Domain topology lines is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain lines
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES                      
                          IF(ASSOCIATED(DOMAIN_LINES)) THEN
                            ALLOCATE(DOMAIN_LINES%LINES(DECOMPOSITION_LINES%NUMBER_OF_LINES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain lines lines.",ERR,ERROR,*999)
                            DOMAIN_LINES%NUMBER_OF_LINES=DECOMPOSITION_LINES%NUMBER_OF_LINES
                            ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                            NODES_NUMBER_OF_LINES=0
                            !Loop over the lines in the topology
                            DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                              DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                element_idx=DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(1)
                                basis_local_line_idx=DECOMPOSITION_LINE%ELEMENT_LINES(1)
                                CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                                DOMAIN_LINE%NUMBER=local_line_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_LINE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                                DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                                ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in line.",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in line.",ERR,ERROR,*999)
                                DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx, &
                                    & basis_local_line_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                    & basis_local_line_node_idx,basis_local_line_idx))
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                                  IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx),element_local_node_idx)
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx),BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx))
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
                                  ENDIF
                                  NODES_NUMBER_OF_LINES(node_idx)=NODES_NUMBER_OF_LINES(node_idx)+1
                                ENDDO !basis_local_line_node_idx
                              ELSE
                                CALL FlagError("Line is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !local_line_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node lines.",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_LINES)
                            DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              BASIS=>DOMAIN_LINE%BASIS
                              DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_LINES=DOMAIN_NODE%NUMBER_OF_NODE_LINES+1
                                DOMAIN_NODE%NODE_LINES(DOMAIN_NODE%NUMBER_OF_NODE_LINES)=local_line_idx
                              ENDDO !basis_local_line_node_idx
                            ENDDO !local_line_idx
                          ELSE
                            CALL FlagError("Domain lines is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
            END IF                       
          ELSE
            CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology lines is not associated.",ERR,ERROR,*999)

      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology lines:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of lines = ",DECOMPOSITION_LINES%NUMBER_OF_LINES,ERR,ERROR,*999)
      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Line number = ",DECOMPOSITION_LINE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",DECOMPOSITION_LINE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%ELEMENT_LINES,'("      Element lines        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_LINE%ADJACENT_LINES, &
          & '("      Adjacent lines       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary line = ",DECOMPOSITION_LINE%BOUNDARY_LINE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_LINE=>DOMAIN%TOPOLOGY%LINES%LINES(local_line_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_LINE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_LINE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in lines = ",DOMAIN_LINE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_LINE%NODES_IN_LINE, &
            & '("        Nodes in line        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_line_node_idx=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_line_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(1,:,basis_local_line_node_idx),'("            Derivatives in line  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(2,:,basis_local_line_node_idx), &
              & '("            Derivatives Versions in line  :",4(X,I8))','(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_line_node_idx
        ENDDO !component_idx
      ENDDO !local_line_idx
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_LINES)) DEALLOCATE(TEMP_LINES)
    IF(ASSOCIATED(NEW_TEMP_LINES)) DEALLOCATE(NEW_TEMP_LINES)
    IF(ALLOCATED(NODES_NUMBER_OF_LINES)) DEALLOCATE(NODES_NUMBER_OF_LINES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given decomposition topology. \todo Pass in the topology lines
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DECOMPOSITION_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines.",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_DataPointsInitialise(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DecompositionTopology_DataPointsInitialise",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Decomposition already has topology data points associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points.",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%numberOfDataPoints=0
        TOPOLOGY%dataPoints%totalNumberOfDataPoints=0
        TOPOLOGY%dataPoints%numberOfGlobalDataPoints=0
        NULLIFY(TOPOLOGY%dataPoints%dataPointsTree)
        TOPOLOGY%dataPoints%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION      
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsInitialise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsInitialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataPointsInitialise
  
  !
  !================================================================================================================================
  !
 
  !>Finalises a face in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(FACE%SURROUNDING_ELEMENTS)) DEALLOCATE(FACE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(FACE%ELEMENT_FACES)) DEALLOCATE(FACE%ELEMENT_FACES)
!    FACE%ADJACENT_FACES=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a decomposition topology face.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!    FACE%ADJACENT_FACES=0
    FACE%BOUNDARY_FACE=.FALSE.
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,ne,surrounding_element_idx,basis_local_face_idx,surrounding_element_basis_local_face_idx, &
      & element_local_node_idx,basis_local_face_node_idx,basis_local_face_derivative_idx,derivative_idx,version_idx,face_idx, &
      & node_idx,elem_idx,NODES_IN_FACE(16),NUMBER_OF_FACES,MAX_NUMBER_OF_FACES,NEW_MAX_NUMBER_OF_FACES,FACE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_FACES(:)
    INTEGER(INTG), POINTER :: TEMP_FACES(:,:),NEW_TEMP_FACES(:,:)
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMPOSITION_FACE!,DECOMPOSITION_FACE2
    TYPE(DECOMPOSITION_FACES_TYPE), POINTER :: DECOMPOSITION_FACES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE!,DOMAIN_FACE2
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_FACES)
    NULLIFY(NEW_TEMP_FACES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_FACES=>TOPOLOGY%FACES
      IF(ASSOCIATED(DECOMPOSITION_FACES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish face
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Estimate the number of faces
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      ! Faces not calculated in 1D 
                    CASE(2)
                      ! Faces not calculated in 2D
                    CASE(3)
                      !This should give the maximum and will over estimate the number of faces for a "cube mesh" by approx 33%
                      MAX_NUMBER_OF_FACES= &
                        & NINT(((REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP)*5.0_DP)+1.0_DP)*(4.0_DP/3.0_DP),INTG)

                      DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                      IF(ASSOCIATED(DOMAIN_FACES)) THEN
                        ALLOCATE(TEMP_FACES(16,MAX_NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate temporary faces array",ERR,ERROR,*999)
                        ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                        NODES_NUMBER_OF_FACES=0
                        NUMBER_OF_FACES=0
                        TEMP_FACES=0
                        !Loop over the elements in the topology and fill temp_faces with node numbers for each element
                        DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_FACES(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate element faces of element",ERR,ERROR,*999)
                          !Loop over the local faces of the element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            !Calculate the topology node numbers that make up the face
                            NODES_IN_FACE=0
                            !Check whether face has already been read out
                            DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                              !Read out node numbers of local face from ELEMENT_NODES
                              NODES_IN_FACE(basis_local_face_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                                & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                            ENDDO !basis_local_face_node_idx
                            !Try and find a previously created face that matches in the adjacent elements
                            FOUND=.FALSE.
                            node_idx=NODES_IN_FACE(1)
                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                              surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                              IF(surrounding_element_idx/=ne) THEN
                                IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES)) THEN
                                  BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                  DO surrounding_element_basis_local_face_idx=1,BASIS2%NUMBER_OF_LOCAL_FACES
                                    face_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES( &
                                      & surrounding_element_basis_local_face_idx)
                                    IF(ALL(NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx))== &
                                      & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx),face_idx))) THEN
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !surrounding_element_basis_local_face_idx
                                  IF(FOUND) EXIT
                                ENDIF
                              ENDIF
                            ENDDO !elem_idx
                            IF(FOUND) THEN
                              !Face has already been created
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=face_idx
                            ELSE
                              !Face has not been created
                              IF(NUMBER_OF_FACES==MAX_NUMBER_OF_FACES) THEN
                                !We are at maximum. Reallocate the FACES array to be 20% bigger and try again.
                                NEW_MAX_NUMBER_OF_FACES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_FACES,DP),INTG)
                                !\todo: Change 16 to a variable and above for NODES_IN_FACE
                                ALLOCATE(NEW_TEMP_FACES(16,NEW_MAX_NUMBER_OF_FACES),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate new number of faces",ERR,ERROR,*999)
                                NEW_TEMP_FACES(:,1:NUMBER_OF_FACES)=TEMP_FACES(:,1:NUMBER_OF_FACES)
                                NEW_TEMP_FACES(:,NUMBER_OF_FACES+1:NEW_MAX_NUMBER_OF_FACES)=0
                                DEALLOCATE(TEMP_FACES)
                                TEMP_FACES=>NEW_TEMP_FACES
                                NULLIFY(NEW_TEMP_FACES)
                                MAX_NUMBER_OF_FACES=NEW_MAX_NUMBER_OF_FACES
                              ENDIF
                              NUMBER_OF_FACES=NUMBER_OF_FACES+1
                              TEMP_FACES(:,NUMBER_OF_FACES)=NODES_IN_FACE(:)
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=NUMBER_OF_FACES
                              DO basis_local_face_node_idx=1,SIZE(NODES_IN_FACE,1)
                                IF(NODES_IN_FACE(basis_local_face_node_idx)/=0) &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))= &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))+1
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        !Allocate the face arrays and set them from the FACES and NODE_FACES arrays
                        DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                          ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate node faces array",ERR,ERROR,*999)
                          DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                        ENDDO !node_idx
                        DEALLOCATE(NODES_NUMBER_OF_FACES)
                        ALLOCATE(DECOMPOSITION_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology faces",ERR,ERROR,*999)
                        DECOMPOSITION_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        ALLOCATE(DOMAIN_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate domain topology faces",ERR,ERROR,*999)
                        DOMAIN_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                          CALL DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(DECOMPOSITION_FACES%FACES(face_idx),ERR,ERROR,*999)
                          CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                          DO basis_local_face_node_idx=1,SIZE(TEMP_FACES,1)
                            IF(TEMP_FACES(basis_local_face_node_idx,face_idx)/=0) THEN
                              node_idx=TEMP_FACES(basis_local_face_node_idx,face_idx)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES+1
                              DOMAIN_NODES%NODES(node_idx)%NODE_FACES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES)=face_idx
                            ENDIF
                          ENDDO !basis_local_face_node_idx  
                        ENDDO !face_idx

                        !Set nodes in face and derivatives of nodes in face for domain faces
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          !Loop over local faces of element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DOMAIN_FACE=>DOMAIN_FACES%FACES(FACE_NUMBER)
                            DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                            IF(.NOT.ASSOCIATED(DOMAIN_FACE%BASIS)) THEN
                              DECOMPOSITION_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%ELEMENT_NUMBER=ne !! Needs checking
!                              DECOMPOSITION_FACE%ELEMENT_NUMBER=DECOMPOSITION_ELEMENT%NUMBER
!                              DOMAIN_FACE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                              DECOMPOSITION_FACE%XI_DIRECTION=BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx)
                              DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(ABS(DECOMPOSITION_FACE%XI_DIRECTION))%PTR
                              ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                & STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face nodes in face",ERR,ERROR,*999)
                              ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face derivatives in face",ERR,ERROR,*999)
                              DOMAIN_FACE%DERIVATIVES_IN_FACE=0
                              !Set nodes in face based upon face number
                              DOMAIN_FACE%NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx))= &
                                & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx),FACE_NUMBER)
                              !Set derivatives of nodes in domain face from derivatives of nodes in element
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx, &
                                  & basis_local_face_idx)
                                !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_FACE( &
                                  & basis_local_face_node_idx,basis_local_face_idx))
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                  DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(1,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(2,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=version_idx
                                  ENDDO !basis_local_face_derivative_idx
                                ENDIF
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        DEALLOCATE(TEMP_FACES)
                        !\todo Note: Adjacency will be left out of faces calculation for the time being
                        !Calculate adjacent faces and the surrounding elements for each face
                        DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                          DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                          DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                          BASIS=>DOMAIN_FACE%BASIS
                          IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                            DECOMPOSITION_FACE%BOUNDARY_FACE=.TRUE.
                            DOMAIN_FACE%BOUNDARY_FACE=.TRUE.
                          ENDIF
                          !Allocate the elements surrounding the face
                          ALLOCATE(DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face surrounding elements",ERR,ERROR,*999)

                          ALLOCATE(DECOMPOSITION_FACE%ELEMENT_FACES(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face element faces",ERR,ERROR,*999)
!                          DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!                          DECOMPOSITION_FACE%ADJACENT_FACES=0

                           !Loop over the nodes at each end of the face
!                          DO node_idx1=0,1
!                           DO node_idx2=0,1
!                            FOUND=.FALSE.
!                            node_idx=DOMAIN_FACE%NODES_IN_FACE((node_idx2*BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION*(BASIS%NUMBER_OF_FACES-1))&
!                                                                             &+(node_idx1*(BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
                             !Loop over the elements surrounding the node.
!                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
!                              ne=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
!                              DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
!                              DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                               !Loop over the local faces of the element
!                              DO basis_local_face_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_FACES
!                                nf2=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
!                                IF(nf2/=face_idx) THEN
!                                  DECOMPOSITION_FACE2=>DECOMPOSITION_FACES%FACES(nf2)
!                                  DOMAIN_FACE2=>DOMAIN_FACES%FACES(nf2)
                                   !Check whether XI of face have same direction
!                                  IF((OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),2,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),2,1)).OR.&
!                                     &(OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),3,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),3,1))) THEN
                                     !Loop over nodes in face of surrounding element
!                                    BASIS2=>DOMAIN_FACE2%BASIS
!                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
!                                      NODE_COUNT=0
!                                      DO node_idx3=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                        DO node_idx4=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                          np2=DOMAIN_FACE2%NODES_IN_FACE((node_idx4*(BASIS2%NUMBER_OF_FACES-1))&
!                                                                      &+(node_idx3*(BASIS2%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
!                                          IF(np2==node_idx) NODE_COUNT=NODE_COUNT+1
!                                        ENDDO !node_idx4
!                                      ENDDO !node_idx3
!                                      IF(NODE_COUNT<BASIS%NUMBER_OF_NODES) THEN
!                                        FOUND=.TRUE.
!                                        EXIT
!                                      ENDIF
!                                    ENDIF
!                                  ENDIF
!                                ENDIF
!                              ENDDO !basis_local_face_idx
!                                IF(FOUND) EXIT
!                            ENDDO !elem_idx
!                            IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx2)=nf2
!                           ENDDO !node_idx2
!                           IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx1)=nf2
!                          ENDDO !node_idx1
                        ENDDO !face_idx

                        !Set the surrounding elements
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DO face_idx=1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS
                              DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(face_idx)=ne
                              DECOMPOSITION_FACE%ELEMENT_FACES(face_idx)=basis_local_face_idx
                            ENDDO
                          ENDDO !basis_local_face_idx
                        ENDDO !ne
                      ELSE
                        CALL FlagError("Domain topology faces is not associated",ERR,ERROR,*999)
                      ENDIF
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain",ERR,ERROR,*999)
                    END SELECT
                 ELSE
                    CALL FlagError("Domain topology elements is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain faces
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                          IF(ASSOCIATED(DOMAIN_FACES)) THEN
                            ALLOCATE(DOMAIN_FACES%FACES(DECOMPOSITION_FACES%NUMBER_OF_FACES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain faces faces",ERR,ERROR,*999)
                            DOMAIN_FACES%NUMBER_OF_FACES=DECOMPOSITION_FACES%NUMBER_OF_FACES
                            ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                            NODES_NUMBER_OF_FACES=0
                            !Loop over the faces in the topology
                            DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                              DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                ne=DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(1)
                                basis_local_face_idx=DECOMPOSITION_FACE%ELEMENT_FACES(1)
                                CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                                DOMAIN_FACE%NUMBER=face_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(ABS(DECOMPOSITION_FACE%XI_DIRECTION))%PTR
                                ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in face",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in face",ERR,ERROR,*999)
                                !Set derivatives of nodes in domain face from derivatives of nodes in element
                                DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx, &
                                    & basis_local_face_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_FACE( &
                                    & basis_local_face_node_idx,basis_local_face_idx))
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                  IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                      derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%DERIVATIVES_IN_FACE(1,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=derivative_idx
                                      version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%DERIVATIVES_IN_FACE(2,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=version_idx
                                    ENDDO !basis_local_face_derivative_idx
                                  ENDIF
                                  NODES_NUMBER_OF_FACES(node_idx)=NODES_NUMBER_OF_FACES(node_idx)+1
                                ENDDO !basis_local_face_node_idx
                              ELSE
                                CALL FlagError("Face is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !face_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node faces",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_FACES)
                            DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              BASIS=>DOMAIN_FACE%BASIS
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_FACES=DOMAIN_NODE%NUMBER_OF_NODE_FACES+1
                                !Set the face numbers a node is on
                                DOMAIN_NODE%NODE_FACES(DOMAIN_NODE%NUMBER_OF_NODE_FACES)=face_idx
                              ENDDO !basis_local_face_node_idx
                            ENDDO !face_idx
                          ELSE
                            CALL FlagError("Domain faces is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
            END IF                       
          ELSE
            CALL FlagError("Topology decomposition is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology faces is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology faces:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of faces = ",DECOMPOSITION_FACES%NUMBER_OF_FACES,ERR,ERROR,*999)
      DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
        DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
        DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Face number = ",DECOMPOSITION_FACE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction (Normal to Face) = &
                                                                         &",DECOMPOSITION_FACE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%ELEMENT_FACES,'("      Element faces        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
!        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_FACE%ADJACENT_FACES, &
!          & '("      Adjacent faces       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary face = ",DECOMPOSITION_FACE%BOUNDARY_FACE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_FACE=>DOMAIN%TOPOLOGY%FACES%FACES(face_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_FACE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_FACE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in faces = ",DOMAIN_FACE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_FACE%NODES_IN_FACE, &
            & '("        Nodes in face        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_face_node_idx=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_face_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(1,:,basis_local_face_node_idx),'("            Derivatives in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(2,:,basis_local_face_node_idx),'("            Derivatives Versions in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_face_node_idx
        ENDDO !component_idx
      ENDDO !face_idx
    ENDIF

    EXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_FACES)) DEALLOCATE(TEMP_FACES)
    IF(ASSOCIATED(NEW_TEMP_FACES)) DEALLOCATE(NEW_TEMP_FACES)
    IF(ALLOCATED(NODES_NUMBER_OF_FACES)) DEALLOCATE(NODES_NUMBER_OF_FACES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DECOMPOSITION_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the decomposition type for a decomposition. \see OPENCMISS::CMISSDecompositionTypeGet
  SUBROUTINE DECOMPOSITION_TYPE_GET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the type for
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the decomposition type for the specified decomposition \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        TYPE=DECOMPOSITION%DECOMPOSITION_TYPE
      ELSE
        CALL FlagError("Decomposition has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_GET",ERR,ERROR)
    RETURN
  END SUBROUTINE DECOMPOSITION_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the decomposition type for a decomposition.  \see OPENCMISS::CMISSDecompositionTypeSet
  SUBROUTINE DECOMPOSITION_TYPE_SET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The decomposition type to set \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          !heye: three types for decomposition--decompostion_all_type means no decomposition
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
        CASE(DECOMPOSITION_CALCULATED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_CALCULATED_TYPE
        CASE(DECOMPOSITION_USER_DEFINED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_USER_DEFINED_TYPE
        CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether lines should be calculated in the the decomposition. \see OPENCMISS::CMISSDecompositionCalculateLinesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET(DECOMPOSITION,CALCULATE_LINES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_LINES_FLAG !<The boolean flag to determine whether the lines should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_LINES=CALCULATE_LINES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_LINES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether faces should be calculated in the the decomposition. \see OPENCMISS::CMISSDecompositionCalculateFacesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET(DECOMPOSITION,CALCULATE_FACES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_FACES_FLAG !<The boolean flag to determine whether the faces should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_FACES=CALCULATE_FACES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_FACES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in DECOMPOSITION a pointer to the decomposition identified by USER_NUMBER in the given MESH. If no decomposition with that USER_NUMBER exists DECOMPOSITION is left nullified.
  SUBROUTINE DECOMPOSITION_USER_NUMBER_FIND(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to find
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh containing the decomposition to find
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On return a pointer to the decomposition with the specified user number. If no decomposition with that user number exists then DECOMPOSITION is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_USER_NUMBER_FIND",ERR,ERROR,*999)

    NULLIFY(DECOMPOSITION)
    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        decomposition_idx=1
        DO WHILE(decomposition_idx<=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.ASSOCIATED(DECOMPOSITION))
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            DECOMPOSITION=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
          ELSE
            decomposition_idx=decomposition_idx+1
          ENDIF
        ENDDO
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_USER_NUMBER_FIND")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_USER_NUMBER_FIND",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises the domain decompositions for a given mesh. \todo Pass in a pointer to the decompositions.
  SUBROUTINE DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise the decomposition for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        DO WHILE(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>0)
          CALL DECOMPOSITION_DESTROY(MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR,ERR,ERROR,*999)
        ENDDO !no_decomposition
       DEALLOCATE(MESH%DECOMPOSITIONS)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain decompositions for a given mesh.
  SUBROUTINE DECOMPOSITIONS_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise the decompositions for

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        CALL FlagError("Mesh already has decompositions associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(MESH%DECOMPOSITIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Mesh decompositions could not be allocated.",ERR,ERROR,*999)
        MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        NULLIFY(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
        MESH%DECOMPOSITIONS%MESH=>MESH
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the domain for a given decomposition and deallocates all memory. \todo Pass in a pointer to the domain
  SUBROUTINE DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the domain for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    
    ENTERS("DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          DO component_idx=1,DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS
            IF(ALLOCATED(DECOMPOSITION%DOMAIN(component_idx)%PTR%NODE_DOMAIN))  &
              & DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR%NODE_DOMAIN)
            CALL DOMAIN_MAPPINGS_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)        
            CALL DOMAIN_TOPOLOGY_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR)
          ENDDO !component_idx
          DEALLOCATE(DECOMPOSITION%DOMAIN)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain for a given decomposition.
  SUBROUTINE DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the domain for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx

    ENTERS("DOMAIN_INITIALISE",ERR,ERROR,*999)
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          CALL FlagError("Decomposition already has a domain associated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Decomposition domain could not be allocated.",ERR,ERROR,*999)
          DO component_idx=1,DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS !Mesh component
            ALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Decomposition domain component could not be allocated.",ERR,ERROR,*999)

            DECOMPOSITION%DOMAIN(component_idx)%PTR%DECOMPOSITION=>DECOMPOSITION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH=>DECOMPOSITION%MESH
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH_COMPONENT_NUMBER=component_idx
            DECOMPOSITION%DOMAIN(component_idx)%PTR%REGION=>DECOMPOSITION%MESH%REGION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_DIMENSIONS=DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS

            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_ELEMENTS=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_FACES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_LINES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_NODES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_MESH_DOFS=0
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%MAPPINGS)
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%TOPOLOGY)

            CALL DOMAIN_MAPPINGS_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            CALL DOMAIN_TOPOLOGY_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
          ENDDO !component_idx
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_INITIALISE

  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs mapping in the given domain mappings. \todo Pass in the domain mappings dofs
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%DOFS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_DOFS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises the dofs mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL FlagError("Domain dofs mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings dofs.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%DOFS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the local/global element mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the element mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,adjacent_element,domain_no,domain_idx,ne,nn,np,NUMBER_OF_DOMAINS, &
      & NUMBER_OF_ADJACENT_ELEMENTS,MY_COMPUTATIONAL_NODE_NUMBER,component_idx
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:),DOMAINS(:),LOCAL_ELEMENT_NUMBERS(:)
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS_LIST(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
          ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
          IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            IF(ASSOCIATED(DOMAIN%MESH)) THEN
              MESH=>DOMAIN%MESH
              component_idx=DOMAIN%MESH_COMPONENT_NUMBER
              MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999        

              !Calculate the local and global numbers and set up the mappings
              ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element mapping global to local map.",ERR,ERROR,*999)
              ELEMENTS_MAPPING%NUMBER_OF_GLOBAL=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS
              !Loop over the global elements and calculate local numbers
              ALLOCATE(LOCAL_ELEMENT_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate local element numbers.",ERR,ERROR,*999)
              LOCAL_ELEMENT_NUMBERS=0
              ALLOCATE(ADJACENT_ELEMENTS_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements list.",ERR,ERROR,*999)
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                NULLIFY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR)
                CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,MAX(INT(MESH%NUMBER_OF_ELEMENTS/2),1), &
                  & ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
              ENDDO !domain_idx
            
              DO ne=1,MESH%NUMBER_OF_ELEMENTS
                !Calculate the local numbers
                domain_no=DECOMPOSITION%ELEMENT_DOMAIN(ne)
                LOCAL_ELEMENT_NUMBERS(domain_no)=LOCAL_ELEMENT_NUMBERS(domain_no)+1
                !Calculate the adjacent elements to the computational domains and the adjacent domain numbers themselves
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                NULLIFY(ADJACENT_DOMAINS_LIST)
                CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)
                DO nn=1,BASIS%NUMBER_OF_NODES
                  np=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                  DO no_adjacent_element=1,MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%numberOfSurroundingElements
                    adjacent_element=MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%surroundingElements(no_adjacent_element)
                    IF(DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element)/=domain_no) THEN
                      CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(domain_no)%PTR,adjacent_element,ERR,ERROR,*999)
                      CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element),ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_adjacent_element
                ENDDO !nn
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                DEALLOCATE(DOMAINS)
                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne),ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map domain number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local type.",ERR,ERROR,*999)
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS=1
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(1)=LOCAL_ELEMENT_NUMBERS(domain_no)
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(1)=DECOMPOSITION%ELEMENT_DOMAIN(ne) 
                IF(NUMBER_OF_DOMAINS==1) THEN
                  !Element is an internal element
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                ELSE
                  !Element is on the boundary of computational domains
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                ENDIF
              ENDDO !ne
              
              !Compute ghost element mappings
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,NUMBER_OF_ADJACENT_ELEMENTS, &
                  & ADJACENT_ELEMENTS,ERR,ERROR,*999)
                DO no_adjacent_element=1,NUMBER_OF_ADJACENT_ELEMENTS
                  adjacent_element=ADJACENT_ELEMENTS(no_adjacent_element)
                  LOCAL_ELEMENT_NUMBERS(domain_idx)=LOCAL_ELEMENT_NUMBERS(domain_idx)+1
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS= &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS+1
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=LOCAL_ELEMENT_NUMBERS(domain_idx)
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%DOMAIN_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=domain_idx
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_TYPE( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)= &
                    & DOMAIN_LOCAL_GHOST
                ENDDO !no_adjacent_element
                IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
              ENDDO !domain_idx
              
              DEALLOCATE(ADJACENT_ELEMENTS_LIST)
              DEALLOCATE(LOCAL_ELEMENT_NUMBERS)
              
              !Calculate element local to global maps from global to local map
              CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(ELEMENTS_MAPPING,ERR,ERROR,*999)
                            
            ELSE
              CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ne=1,MESH%NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER, &
          & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ne
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ne=1,ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global element = ", &
          & ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ne),ERR,ERROR,*999)
      ENDDO !ne
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal elements = ", &
          & ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%INTERNAL_START:ELEMENTS_MAPPING%INTERNAL_FINISH), &
          & '("    Internal elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary elements = ", &
          & ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%BOUNDARY_START:ELEMENTS_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost elements = ", &
          & ELEMENTs_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%GHOST_START:ELEMENTS_MAPPING%GHOST_FINISH), &
          & '("    Ghost elements   :",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR( &
        & ELEMENTS_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,ELEMENTS_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)    
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE
  
!================================================================================================================================
  !>Calculates the local/global node/DOFs mappings for a domain decomposition.
!================================================================================================================================
 SUBROUTINE DOMAIN_MAPPINGS_DOFS_NODE_CALCULATE_NEW(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the node mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: adjacent_node, adjacent_node_idx, adjacent_element,  & 
        & component_idx, derivative_idx, DUMMY_ERR, domain_idx, domain_no, MY_COMPUTATIONAL_NODE_NUMBER, ne, nn, np, ny, &
          & NUMBER_OF_DOMAINS, NUMBER_OF_ADJACENT_ELEMENTS, NUMBER_OF_ADJACENT_NODES, &
            & version_idx, node_idx, no_adjacent_node, no_adjacent_element
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_NOdes(:),DOMAINS(:),LOCAL_NODE_NUMBERS(:),LOCAL_DOF_NUMBERS(:)
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_NODES_LIST(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS_MAPPING,NODE_MAPPING
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("DOMAIN_MAPPINGS_DOFS_NODE_CALCULATE_NEW",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) THEN
          NODE_MAPPING=>DOMAIN%MAPPINGS%NODES
          IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) THEN
            DOFS_MAPPING=>DOMAIN%MAPPINGS%DOFS
            IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
              DECOMPOSITION=>DOMAIN%DECOMPOSITION
                IF(ASSOCIATED(DOMAIN%MESH)) THEN

                  MESH=>DOMAIN%MESH
                  component_idx=DOMAIN%MESH_COMPONENT_NUMBER
                  MESH_TOPOLOGY=>MESH%TOPOLOGY(component_idx)%PTR
                  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999 
                  MESH%NUMBER_OF_NODES = MESH%TOPOLOGY(1)%PTR%NODES%numberOfnodes  
                  IF(ALLOCATED(DOMAIN%NODE_DOMAIN)) DEALLOCATE(DOMAIN%NODE_DOMAIN)
                  ALLOCATE(DOMAIN%NODE_DOMAIN(MESH%TOPOLOGY(1)%PTR%NODES%numberOfnodes)) 
                  DOMAIN%NODE_DOMAIN = DECOMPOSITION%NODE_DOMAIN
                  !Calculate the local and global numbers and set up the mappings
                  ALLOCATE(NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate Node mapping global to local map.",ERR,ERROR,*999)
                  NODE_MAPPING%NUMBER_OF_GLOBAL=MESH%TOPOLOGY(component_idx)%PTR%NODES%numberOfnodes
                  !Loop over the global elements and calculate local numbers
                  ALLOCATE(LOCAL_NODE_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local node numbers.",ERR,ERROR,*999)
                  LOCAL_NODE_NUMBERS=0
                  ALLOCATE(ADJACENT_NODES_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate adjacent node list.",ERR,ERROR,*999)
                  ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%dofs%numberOfDofs),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate dofs mapping global to local map.",ERR,ERROR,*999)
                  DOFS_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%DOFS%numberOfDofs
                  ALLOCATE(LOCAL_DOF_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local dof numbers.",ERR,ERROR,*999)
                  LOCAL_DOF_NUMBERS=0
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    NULLIFY(ADJACENT_NODES_LIST(domain_idx)%PTR)
                    CALL LIST_CREATE_START(ADJACENT_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ADJACENT_NODES_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ADJACENT_NODES_LIST(domain_idx)%PTR,MAX(INT(MESH%NUMBER_OF_NODES/2),1), &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ADJACENT_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !domain_idx
                  DO ne=1,SIZE(MESH_TOPOLOGY%NODES%NODES)
                   !Calculate the local numbers
                   domain_no=DECOMPOSITION%NODE_DOMAIN(ne-1) !look here
                   LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                   LOCAL_DOF_NUMBERS(domain_no)= LOCAL_DOF_NUMBERS(domain_no)+1
                   !Calculate the adjacent elements to the computational domains and the adjacent domain numbers themselves
                   NULLIFY(ADJACENT_DOMAINS_LIST)
                   CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                   CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                   CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                   CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                   CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)

                   DO no_adjacent_element=1,MESH_TOPOLOGY%NODES%NODES(ne)%numberOfSurroundingElements

                     adjacent_element=MESH_TOPOLOGY%NODES%NODES(ne)%surroundingElements(no_adjacent_element)

                     DO no_adjacent_node=1, SIZE(MESH_TOPOLOGY%ELEMENTS%ELEMENTS(adjacent_element)%MESH_ELEMENT_NODES)

                       adjacent_node = &
                         & MESH_TOPOLOGY%ELEMENTS%ELEMENTS(adjacent_element)%MESH_ELEMENT_NODES(no_adjacent_node)
  
                       domain_no=DECOMPOSITION%NODE_DOMAIN(ne-1)
          
                       IF(DECOMPOSITION%NODE_DOMAIN(adjacent_node-1)/=domain_no) THEN
                         CALL LIST_ITEM_ADD(ADJACENT_NODES_LIST(domain_no)%PTR,adjacent_node,ERR,ERROR,*999)
                         CALL LIST_ITEM_ADD& 
                           & (ADJACENT_DOMAINS_LIST,DECOMPOSITION%NODE_DOMAIN(adjacent_node-1),ERR,ERROR,*999)
                       END IF

                     ENDDO !no_adjacent_node
                   ENDDO !no_adjacent_element

                   CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                   CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                   DEALLOCATE(DOMAINS)
                   CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne),ERR,ERROR,*999)
                   ALLOCATE(NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)

                   IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local number.",ERR,ERROR,*999)
                   ALLOCATE(NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                   IF(ERR/=0) CALL FlagError("Could not allocate element global to local map domain number.",ERR,ERROR,*999)
                   ALLOCATE(NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
                   IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local type.",ERR,ERROR,*999)
                   DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%numberOfDerivatives
                     DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%numberOfVersions
                       ny=MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                       ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                       IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local number.",ERR,ERROR,*999)
                       ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                       IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map domain number.",ERR,ERROR,*999)
                       ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
                       IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local type.",ERR,ERROR,*999)
                     ENDDO !version_idx
                   ENDDO !derivative_idx

                   NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS=1
                   NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                   NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(1)=DECOMPOSITION%NODE_DOMAIN(ne-1)

                   DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%numberOfDerivatives
                     DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%numberOfVersions
                       ny=MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                       DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                       DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                       DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=DECOMPOSITION%NODE_DOMAIN(ne-1)
                     ENDDO !version_idx
                   ENDDO !derivative_idx

                   IF(NUMBER_OF_DOMAINS==1) THEN
                     !Node/Dof is an internal element
                     NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                     DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%numberOfDerivatives
                       DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%numberOfVersions
                         ny=MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                         DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                       ENDDO !version_idx
                     ENDDO !derivative_idx
                   ELSE
                     !Node/Dof is on the boundary of computational domains
                     NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                     DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%numberOfDerivatives
                       DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%numberOfVersions
                         ny=MESH_TOPOLOGY%NODES%NODES(ne)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                         DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                       ENDDO !version_idx
                     ENDDO !derivative_idx
                   ENDIF
                 ENDDO !ne

                 DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1

                   CALL LIST_REMOVE_DUPLICATES(ADJACENT_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)

                   CALL LIST_DETACH_AND_DESTROY(ADJACENT_NODES_LIST(domain_idx)%PTR,NUMBER_OF_ADJACENT_NODES, &
                     & ADJACENT_NODES,ERR,ERROR,*999)

                   DO no_adjacent_node=1,NUMBER_OF_ADJACENT_NODES

                     adjacent_node=ADJACENT_NODES(no_adjacent_node)

                     LOCAL_NODE_NUMBERS(domain_idx)=LOCAL_NODE_NUMBERS(domain_idx)+1

                     NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%NUMBER_OF_DOMAINS= &
                       & NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%NUMBER_OF_DOMAINS+1

                     NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%LOCAL_NUMBER( &
                       & NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%NUMBER_OF_DOMAINS)=LOCAL_NODE_NUMBERS(domain_idx)

                     NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%DOMAIN_NUMBER( &
                       & NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%NUMBER_OF_DOMAINS)=domain_idx

                     NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%LOCAL_TYPE( &
                       & NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%NUMBER_OF_DOMAINS)= &
                         & DOMAIN_LOCAL_GHOST

                     DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(adjacent_node)%numberOfDerivatives

                       DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(adjacent_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                      
                         ny=MESH_TOPOLOGY%NODES%NODES(adjacent_node)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)

                         LOCAL_DOF_NUMBERS(domain_idx)=LOCAL_DOF_NUMBERS(domain_idx)+1

                         DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS= &
                           & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS+1

                         DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER( &
                           & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                            & LOCAL_DOF_NUMBERS(domain_idx)

                         DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER( &
                           & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)=domain_idx

                         DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE( &
                           & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                             & DOMAIN_LOCAL_GHOST

                       ENDDO !version_idx
                     ENDDO !derivative_idx

                   ENDDO !no_adjacent_node
                   IF(ALLOCATED(ADJACENT_NODES)) DEALLOCATE(ADJACENT_NODES)
                 ENDDO !domain_idx

                 DEALLOCATE(ADJACENT_NODES_LIST)
                 DEALLOCATE(LOCAL_NODE_NUMBERS)
                 !Calculate element local to global maps from global to local map
                 CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(NODE_MAPPING,ERR,ERROR,*999)
                 CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOFS_MAPPING,ERR,ERROR,*999)
         
               ELSE
                 CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
               ENDIF
             ELSE
               CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
             ENDIF
           ELSE
             CALL FlagError("Domain mappings dofs is not associated.",ERR,ERROR,*999)
           ENDIF
         ELSE
           CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
         ENDIF
       ELSE
         CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
       ENDIF
     ELSE
       CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
     ENDIF
   
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node decomposition :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Domain = ",DOMAIN%NODE_DOMAIN(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
            & NUMBER_OF_DOMAINS,8,8,NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODE_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !node_idx     
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO node_idx=1,NODE_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global node = ", &
          & NODE_MAPPING%LOCAL_TO_GLOBAL_MAP(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal nodes = ", &
          & NODE_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & NODE_MAPPING%DOMAIN_LIST(NODE_MAPPING%INTERNAL_START:NODE_MAPPING%INTERNAL_FINISH), &
          & '("    Internal nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary nodes = ", &
          & NODE_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & NODE_MAPPING%DOMAIN_LIST(NODE_MAPPING%BOUNDARY_START:NODE_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost nodes = ", &
          & NODE_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%NUMBER_OF_GHOST,8,8, &
          & NODE_MAPPING%DOMAIN_LIST(NODE_MAPPING%GHOST_START:NODE_MAPPING%GHOST_FINISH), &
          & '("    Ghost nodes   :",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & NODE_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & NODE_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%ADJACENT_DOMAINS_PTR( &
        & NODE_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,NODE_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,NODE_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & NODE_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & NODE_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,NODE_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & NODE_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODE_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,NODE_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Dofs mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ny=1,MESH_TOPOLOGY%DOFS%numberOfDofs
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
            & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ny   
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ny=1,DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global dof = ", &
          & DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(ny),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(.false.) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%INTERNAL_START:DOFS_MAPPING%INTERNAL_FINISH), &
          & '("    Internal dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%BOUNDARY_START:DOFS_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%GHOST_START:DOFS_MAPPING%GHOST_FINISH), &
          & '("    Ghost dofs   :",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & DOFS_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS_PTR( &
        & DOFS_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,DOFS_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx

    ENDIF

    EXITS("DOMAIN_MAPPINGS_NODES_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ADJACENT_NODES)) DEALLOCATE(ADJACENT_NODES)    
    IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("DOMAIN_MAPPINGS_DOFS_NODE_CALCULATE_NEW",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_NODE_CALCULATE_NEW

!=====================================================================================================================

  !>Finalises the mappings in the given domain. \todo pass in the domain mappings
  SUBROUTINE DOMAIN_MAPPINGS_FINALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%MAPPINGS)
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the element mapping in the given domain mapping. \todo pass in the domain mappings elements
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the elements for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%ELEMENTS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL FlagError("Domain elements mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings elements.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%ELEMENTS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the mappings for a domain decomposition. \todo finalise on error.
  SUBROUTINE DOMAIN_MAPPINGS_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to initialise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        CALL FlagError("Domain already has mappings associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN%MAPPINGS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings.",ERR,ERROR,*999)
        DOMAIN%MAPPINGS%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%MAPPINGS%ELEMENTS)
        NULLIFY(DOMAIN%MAPPINGS%NODES)
        NULLIFY(DOMAIN%MAPPINGS%DOFS)
        !Calculate the node and element mappings

        CALL DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)

        IF(DOMAIN%DECOMPOSITION%NODE_BASED_DECOMPOSITION .EQV. .FALSE.) THEN ! For element-wise decomposition.

           CALL DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*999)
           CALL DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*999)

        ELSE                                                                ! For vertex-wise decomposition.

           CALL DOMAIN_MAPPINGS_DOFS_NODE_CALCULATE_NEW(DOMAIN,ERR,ERROR,*999)
           CALL DOMAIN_MAPPINGS_ELEMENTS_CALCULATE_NEW(DOMAIN,ERR,ERROR,*999)
        END IF

      END IF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_INITIALISE
  
!==============================================================================================================================
  !>Calculates the local/global element mappings for a domain decomposition.
!==============================================================================================================================
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE_NEW(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: adjacent_node, adjacent_element, component_idx, & 
      & derivative_idx, decomposition_idx, domain_idx, domain_idx2, domain_no, DUMMY_ERR, & 
        & element_idx, ghost_element,MY_COMPUTATIONAL_NODE_NUMBER, MAX_NUMBER_DOMAINS, ne, no_adjacent_node, &
          & no_ghost_element, no_adjacent_element, no_computational_node, NUMBER_OF_DOMAINS, & 
            & NUMBER_OF_ELEMENTS_PER_DOMAIN, NUMBER_OF_GHOST_ELEMENTS, ny, number_computational_nodes, version_idx
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_COUNT(:),LOCAL_ELEMENT_NUMBERS(:),LOCAL_DOF_NUMBERS(:), & 
      & NUMBER_INTERNAL_ELEMENTS(:), NUMBER_BOUNDARY_ELEMENTS(:)
    INTEGER(INTG), ALLOCATABLE :: ALL_DOMAINS(:), DOMAINS(:), GHOST_ELEMENTS(:)
    LOGICAL :: BOUNDARY_DOMAIN
    TYPE(LIST_TYPE), POINTER ::   ADJACENT_DOMAINS_LIST, ALL_ADJACENT_DOMAINS_LIST 
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_ELEMENTS_LIST(:)
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE_NEW",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
          ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
           IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) THEN
             NODES_MAPPING=>DOMAIN%MAPPINGS%NODES
              IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
                DECOMPOSITION=>DOMAIN%DECOMPOSITION
                IF(ASSOCIATED(DOMAIN%MESH)) THEN

                  MESH=>DOMAIN%MESH
                  component_idx=DOMAIN%MESH_COMPONENT_NUMBER
                  MESH_TOPOLOGY=>MESH%TOPOLOGY(component_idx)%PTR
                  number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999

                  !Calculate the local and global numbers and set up the mappings
                  ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate element mapping global to local map.",ERR,ERROR,*999)
                  ELEMENTS_MAPPING%NUMBER_OF_GLOBAL=MESH%NUMBER_OF_ELEMENTS
                  ALLOCATE(LOCAL_ELEMENT_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local element numbers.",ERR,ERROR,*999)
                  LOCAL_ELEMENT_NUMBERS=0
                  ALLOCATE(GHOST_ELEMENTS_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ghost element list.",ERR,ERROR,*999)
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    NULLIFY(GHOST_ELEMENTS_LIST(domain_idx)%PTR)
                    CALL LIST_CREATE_START(GHOST_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(GHOST_ELEMENTS_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(GHOST_ELEMENTS_LIST(domain_idx)%PTR,CEILING(MESH%NUMBER_OF_ELEMENTS/2.0), &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(GHOST_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !domain_idx

                  ALLOCATE(NUMBER_INTERNAL_ELEMENTS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of internal element.",ERR,ERROR,*999)
                  NUMBER_INTERNAL_ELEMENTS=0
                  ALLOCATE(NUMBER_BOUNDARY_ELEMENTS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of boundary element.",ERR,ERROR,*999)
                  NUMBER_BOUNDARY_ELEMENTS=0    

                  !For the first pass just determine the internal and boundary elements
                  DO element_idx=1,MESH%NUMBER_OF_ELEMENTS

                    NULLIFY(ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    NULLIFY(ALL_ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ALL_ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ALL_ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)

                    DO no_adjacent_node=1, SIZE(MESH_TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%MESH_ELEMENT_NODES,1)

                     adjacent_node=MESH_TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%MESH_ELEMENT_NODES(no_adjacent_node)

                     DO no_adjacent_element=1,MESH_TOPOLOGY%NODES%NODES(adjacent_node)%numberOfSurroundingElements

                       domain_no=DECOMPOSITION%NODE_DOMAIN(adjacent_node-1)

                       CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)
                       DO domain_idx=1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)%NUMBER_OF_DOMAINS
                         CALL LIST_ITEM_ADD(ALL_ADJACENT_DOMAINS_LIST,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_node)% &
                           & DOMAIN_NUMBER(domain_idx),ERR,ERROR,*999)
                       ENDDO !domain_idx

                     END DO !no_adjacent_element
                    END DO  !no_adjacent_node

                    CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                    CALL LIST_REMOVE_DUPLICATES(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ALL_ADJACENT_DOMAINS_LIST,MAX_NUMBER_DOMAINS,ALL_DOMAINS,ERR,ERROR,*999)

                    IF(NUMBER_OF_DOMAINS/=MAX_NUMBER_DOMAINS) THEN !Ghost element
                      DO domain_idx=1,MAX_NUMBER_DOMAINS
                        domain_no=ALL_DOMAINS(domain_idx)
                        BOUNDARY_DOMAIN=.FALSE.
                        DO domain_idx2=1,NUMBER_OF_DOMAINS
                          IF(domain_no==DOMAINS(domain_idx2)) THEN
                            BOUNDARY_DOMAIN=.TRUE.
                            EXIT
                          ENDIF
                        ENDDO !domain_idx2
                        IF(.NOT.BOUNDARY_DOMAIN) then 
                          CALL LIST_ITEM_ADD(GHOST_ELEMENTS_LIST(domain_no)%PTR,element_idx,ERR,ERROR,*999)
                        end if
                      ENDDO !domain_idx
                    ENDIF
                    ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(element_idx)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local number.",ERR,ERROR,*999)
                    ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(element_idx)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate element global to local map domain number.",ERR,ERROR,*999)
                    ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(element_idx)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local type.",ERR,ERROR,*999)
                    IF(NUMBER_OF_DOMAINS==1) THEN
                      !Node is an internal node
                      domain_no=DOMAINS(1)
                      NUMBER_INTERNAL_ELEMENTS(domain_no)=NUMBER_INTERNAL_ELEMENTS(domain_no)+1
                      !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%NUMBER_OF_DOMAINS=1
                      !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(DOMAINS(1))
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(element_idx)%LOCAL_NUMBER(1)=-1
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(element_idx)%DOMAIN_NUMBER(1)=DOMAINS(1) 
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(element_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                    ELSE

                      !Node is on the boundary of computational domains
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=DOMAINS(domain_idx)
                        !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                        !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(domain_idx)=LOCAL_NODE_NUMBERS(domain_no)
                        ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%LOCAL_NUMBER(domain_idx)=-1
                        ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%DOMAIN_NUMBER(domain_idx)=domain_no
                        ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                      ENDDO !domain_idx
                    ENDIF
                    DEALLOCATE(DOMAINS) 
                    DEALLOCATE(ALL_DOMAINS)
                  ENDDO !node_idx
                  !For the second pass assign boundary nodes to one domain on the boundary and set local node numbers.
                  NUMBER_OF_ELEMENTS_PER_DOMAIN=FLOOR(REAL(MESH%NUMBER_OF_ELEMENTS,DP)/ &
                    & REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,DP))
                  IF(ALLOCATED(DECOMPOSITION%ELEMENT_DOMAIN)) DEALLOCATE(DECOMPOSITION%ELEMENT_DOMAIN)
                  ALLOCATE(DECOMPOSITION%ELEMENT_DOMAIN(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)



                  IF(ERR/=0) CALL FlagError("Could not allocate ELEMENT domain",ERR,ERROR,*999)
                  DECOMPOSITION%ELEMENT_DOMAIN=-1
                  DO element_idx=1,MESH%NUMBER_OF_ELEMENTS
                    IF(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%NUMBER_OF_DOMAINS==1) THEN !Internal node
                      domain_no=ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%DOMAIN_NUMBER(1)
                      DECOMPOSITION%ELEMENT_DOMAIN(ELEMENT_idx)=domain_no
                      LOCAL_ELEMENT_NUMBERS(domain_no)=LOCAL_ELEMENT_NUMBERS(domain_no)+1
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%LOCAL_NUMBER(1)=LOCAL_ELEMENT_NUMBERS(domain_no)
                    ELSE !Boundary node
                      NUMBER_OF_DOMAINS=ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%NUMBER_OF_DOMAINS
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%DOMAIN_NUMBER(domain_idx)
                        IF( DECOMPOSITION%ELEMENT_DOMAIN(ELEMENT_idx)<0) THEN
                          IF((NUMBER_INTERNAL_ELEMENTS(domain_no)+ & 
                            & NUMBER_BOUNDARY_ELEMENTS(domain_no)<NUMBER_OF_ELEMENTS_PER_DOMAIN).OR. &
                              & (domain_idx==ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%NUMBER_OF_DOMAINS)) THEN
                            !Allocate the node to this domain
                            DECOMPOSITION%ELEMENT_DOMAIN(ELEMENT_idx)=domain_no
                            NUMBER_BOUNDARY_ELEMENTS(domain_no)=NUMBER_BOUNDARY_ELEMENTS(domain_no)+1
                            LOCAL_ELEMENT_NUMBERS(domain_no)=LOCAL_ELEMENT_NUMBERS(domain_no)+1
                            !Reset the boundary information to be in the first domain index. The remaining domain indicies will
                            !be overwritten when the ghost nodes are calculated below. 
                            ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%NUMBER_OF_DOMAINS=1
                            ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%LOCAL_NUMBER(1)=LOCAL_ELEMENT_NUMBERS(domain_no)
                            ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%DOMAIN_NUMBER(1)=domain_no
                            ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ELEMENT_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                          ELSE
                            !The node as already been assigned to a domain so it must be a ghost node in this domain
                            CALL LIST_ITEM_ADD(GHOST_ELEMENTS_LIST(domain_no)%PTR,ELEMENT_idx,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          !The node as already been assigned to a domain so it must be a ghost node in this domain
                          CALL LIST_ITEM_ADD(GHOST_ELEMENTS_LIST(domain_no)%PTR,ELEMENT_idx,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !domain_idx
                    ENDIF
                  ENDDO !node_idx
                  DEALLOCATE(NUMBER_INTERNAL_ELEMENTS)

                 !Calculate ghost node and dof mappings
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    CALL LIST_REMOVE_DUPLICATES(GHOST_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(&
                      & GHOST_ELEMENTS_LIST(domain_idx)%PTR,NUMBER_OF_GHOST_ELEMENTS,GHOST_ELEMENTS,ERR,ERROR,*999)

                    DO no_ghost_element=1,NUMBER_OF_GHOST_ELEMENTS
                      ghost_element=GHOST_ELEMENTS(no_ghost_element)
                      LOCAL_ELEMENT_NUMBERS(domain_idx)=LOCAL_ELEMENT_NUMBERS(domain_idx)+1
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%NUMBER_OF_DOMAINS= &
                        & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%NUMBER_OF_DOMAINS+1

                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%LOCAL_NUMBER( &
                        & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%NUMBER_OF_DOMAINS)= &
                        & LOCAL_ELEMENT_NUMBERS(domain_idx)
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%DOMAIN_NUMBER( &
                        & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%NUMBER_OF_DOMAINS)=domain_idx
                      ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%LOCAL_TYPE( &
                        & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_element)%NUMBER_OF_DOMAINS)= &
                        & DOMAIN_LOCAL_GHOST
                    ENDDO !no_ghost_node
                    DEALLOCATE(GHOST_ELEMENTS)
                  ENDDO !domain_idx

                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(ELEMENTS_MAPPING,ERR,ERROR,*999)



                  IF(DIAGNOSTICS1) THEN
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element mappings :",ERR,ERROR,*999)
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
                    DO ne=1,MESH%NUMBER_OF_ELEMENTS
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global element = ",ne,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
                        & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
                        & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER, &
                          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
                        & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER, &
                          & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
                        & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE, &
                          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
                    ENDDO !ne
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
                    DO ne=1,ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local element = ",ne,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global element = ", &
                        & ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ne),ERR,ERROR,*999)
                    ENDDO !ne
                    IF(DIAGNOSTICS2) THEN
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal elements :",ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal elements = ", &
                        & ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
                        & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%INTERNAL_START:ELEMENTS_MAPPING%INTERNAL_FINISH), &
                          & '("    Internal elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary elements :",ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary elements = ", &
                        & ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
                        & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%BOUNDARY_START:ELEMENTS_MAPPING%BOUNDARY_FINISH), &
                          & '("    Boundary elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost elements :",ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost elements = ", &
                        & ELEMENTs_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_GHOST,8,8, &
                        & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%GHOST_START:ELEMENTS_MAPPING%GHOST_FINISH), &
                          & '("    Ghost elements   :",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
                    ENDIF
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
                      & ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
                      & ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR, & 
                        & '("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR( &
                      & ELEMENTS_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,ELEMENTS_MAPPING%ADJACENT_DOMAINS_LIST, &
                        & '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
                    DO domain_idx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
                        & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
                        & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
                        & NUMBER_OF_SEND_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
                        & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
                        & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
                        & NUMBER_OF_RECEIVE_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
                         & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
                    ENDDO !domain_idx
                  ENDIF

                ELSE
                  CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
              ENDIF
           ELSE
             CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
           ENDIF
        ELSE
          CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ALL_DOMAINS)) DEALLOCATE(ALL_DOMAINS)
    IF(ALLOCATED(GHOST_ELEMENTS)) DEALLOCATE(GHOST_ELEMENTS)
    IF(ALLOCATED(NUMBER_INTERNAL_ELEMENTS)) DEALLOCATE(NUMBER_INTERNAL_ELEMENTS)
    IF(ALLOCATED(NUMBER_BOUNDARY_ELEMENTS)) DEALLOCATE(NUMBER_BOUNDARY_ELEMENTS)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE_NEW",ERR,ERROR)
    RETURN 1

  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE_NEW

!==============================================================================================================================
  !>Calculates the local/global node and dof mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,no_computational_node,no_ghost_node,adjacent_element,ghost_node, &
      & NUMBER_OF_NODES_PER_DOMAIN,domain_idx,domain_idx2,domain_no,node_idx,derivative_idx,version_idx,ny,NUMBER_OF_DOMAINS, &
      & MAX_NUMBER_DOMAINS,NUMBER_OF_GHOST_NODES,MY_COMPUTATIONAL_NODE_NUMBER,number_computational_nodes,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_NODE_NUMBERS(:),LOCAL_DOF_NUMBERS(:),NODE_COUNT(:),NUMBER_INTERNAL_NODES(:), &
      & NUMBER_BOUNDARY_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: DOMAINS(:),ALL_DOMAINS(:),GHOST_NODES(:)
    LOGICAL :: BOUNDARY_DOMAIN
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST,ALL_ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_NODES_LIST(:)
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) THEN
          NODES_MAPPING=>DOMAIN%MAPPINGS%NODES
          IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) THEN
            DOFS_MAPPING=>DOMAIN%MAPPINGS%DOFS
            IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
              ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
              IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
                DECOMPOSITION=>DOMAIN%DECOMPOSITION
                IF(ASSOCIATED(DOMAIN%MESH)) THEN
                  MESH=>DOMAIN%MESH
                  component_idx=DOMAIN%MESH_COMPONENT_NUMBER
                  MESH_TOPOLOGY=>MESH%TOPOLOGY(component_idx)%PTR
                  
                  number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  
                  !Calculate the local and global numbers and set up the mappings
                  ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%NODES%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node mapping global to local map.",ERR,ERROR,*999)
                  NODES_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%NODES%numberOfNodes
                  ALLOCATE(LOCAL_NODE_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local node numbers.",ERR,ERROR,*999)
                  LOCAL_NODE_NUMBERS=0
                  ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%dofs%numberOfDofs),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate dofs mapping global to local map.",ERR,ERROR,*999)
                  DOFS_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%DOFS%numberOfDofs
                  ALLOCATE(LOCAL_DOF_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local dof numbers.",ERR,ERROR,*999)
                  LOCAL_DOF_NUMBERS=0
                  ALLOCATE(GHOST_NODES_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ghost nodes list.",ERR,ERROR,*999)
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    NULLIFY(GHOST_NODES_LIST(domain_idx)%PTR)
                    CALL LIST_CREATE_START(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(GHOST_NODES_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(GHOST_NODES_LIST(domain_idx)%PTR,INT(MESH_TOPOLOGY%NODES%numberOfNodes/2), &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !domain_idx
                  ALLOCATE(NUMBER_INTERNAL_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of internal nodes.",ERR,ERROR,*999)
                  NUMBER_INTERNAL_NODES=0
                  ALLOCATE(NUMBER_BOUNDARY_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of boundary nodes.",ERR,ERROR,*999)
                  NUMBER_BOUNDARY_NODES=0

                  !For the first pass just determine the internal and boundary nodes
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    NULLIFY(ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    NULLIFY(ALL_ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ALL_ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ALL_ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    DO no_adjacent_element=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfSurroundingElements
                      adjacent_element=MESH_TOPOLOGY%NODES%NODES(node_idx)%surroundingElements(no_adjacent_element)
                      domain_no=DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element)
                      CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)
                      DO domain_idx=1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS
                        CALL LIST_ITEM_ADD(ALL_ADJACENT_DOMAINS_LIST,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)% &
                          & DOMAIN_NUMBER(domain_idx),ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDDO !no_adjacent_element
                    CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                    CALL LIST_REMOVE_DUPLICATES(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ALL_ADJACENT_DOMAINS_LIST,MAX_NUMBER_DOMAINS,ALL_DOMAINS,ERR,ERROR,*999)
                    IF(NUMBER_OF_DOMAINS/=MAX_NUMBER_DOMAINS) THEN !Ghost node
                      DO domain_idx=1,MAX_NUMBER_DOMAINS
                        domain_no=ALL_DOMAINS(domain_idx)
                        BOUNDARY_DOMAIN=.FALSE.
                        DO domain_idx2=1,NUMBER_OF_DOMAINS
                          IF(domain_no==DOMAINS(domain_idx2)) THEN
                            BOUNDARY_DOMAIN=.TRUE.
                            EXIT
                          ENDIF
                        ENDDO !domain_idx2
                        IF(.NOT.BOUNDARY_DOMAIN) CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDIF
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map local number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map domain number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map local type.",ERR,ERROR,*999)
                    DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                      DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                        ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map domain number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local type.",ERR,ERROR,*999)
                      ENDDO !version_idx
                    ENDDO !derivative_idx
                    IF(NUMBER_OF_DOMAINS==1) THEN
                      !Node is an internal node
                      domain_no=DOMAINS(1)
                      NUMBER_INTERNAL_NODES(domain_no)=NUMBER_INTERNAL_NODES(domain_no)+1
                      !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=1
                      !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(DOMAINS(1))
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=-1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)=DOMAINS(1) 
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=-1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ELSE
                      !Node is on the boundary of computational domains
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=DOMAINS(domain_idx)
                        !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                        !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(domain_idx)=LOCAL_NODE_NUMBERS(domain_no)
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(domain_idx)=-1
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(domain_idx)=domain_no
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                        DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                          DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                            ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)=-1
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)=domain_no
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                          ENDDO !version_idx
                        ENDDO !derivative_idx
                      ENDDO !domain_idx
                    ENDIF
                    DEALLOCATE(DOMAINS) 
                    DEALLOCATE(ALL_DOMAINS)
                  ENDDO !node_idx

                  !For the second pass assign boundary nodes to one domain on the boundary and set local node numbers.
                  NUMBER_OF_NODES_PER_DOMAIN=FLOOR(REAL(MESH_TOPOLOGY%NODES%numberOfNodes,DP)/ &
                    & REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,DP))
                  IF(ALLOCATED(DOMAIN%NODE_DOMAIN)) DEALLOCATE(DOMAIN%NODE_DOMAIN) !mirzawd
                  ALLOCATE(DOMAIN%NODE_DOMAIN(MESH_TOPOLOGY%NODES%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node domain",ERR,ERROR,*999)
                  DOMAIN%NODE_DOMAIN=-1
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    IF(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS==1) THEN !Internal node
                      domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)
                      DOMAIN%NODE_DOMAIN(node_idx)=domain_no
                      LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ELSE !Boundary node
                      NUMBER_OF_DOMAINS=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(domain_idx)
                        IF(DOMAIN%NODE_DOMAIN(node_idx)<0) THEN
                          IF((NUMBER_INTERNAL_NODES(domain_no)+NUMBER_BOUNDARY_NODES(domain_no)<NUMBER_OF_NODES_PER_DOMAIN).OR. &
                            & (domain_idx==NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS)) THEN
                            !Allocate the node to this domain
                            DOMAIN%NODE_DOMAIN(node_idx)=domain_no
                            NUMBER_BOUNDARY_NODES(domain_no)=NUMBER_BOUNDARY_NODES(domain_no)+1
                            LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                            !Reset the boundary information to be in the first domain index. The remaining domain indicies will
                            !be overwritten when the ghost nodes are calculated below. 
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=1
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)=domain_no
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                            DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                              DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                                ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                                LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                              ENDDO !version_idx
                            ENDDO !derivative_idx
                          ELSE
                            !The node as already been assigned to a domain so it must be a ghost node in this domain
                            CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          !The node as already been assigned to a domain so it must be a ghost node in this domain
                          CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !domain_idx
                    ENDIF
                  ENDDO !node_idx
                  DEALLOCATE(NUMBER_INTERNAL_NODES)
                  
                  !Calculate ghost node and dof mappings
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    CALL LIST_REMOVE_DUPLICATES(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(GHOST_NODES_LIST(domain_idx)%PTR,NUMBER_OF_GHOST_NODES,GHOST_NODES,ERR,ERROR,*999)
                    DO no_ghost_node=1,NUMBER_OF_GHOST_NODES
                      ghost_node=GHOST_NODES(no_ghost_node)
                      LOCAL_NODE_NUMBERS(domain_idx)=LOCAL_NODE_NUMBERS(domain_idx)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS= &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%LOCAL_NUMBER( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)= &
                        & LOCAL_NODE_NUMBERS(domain_idx)
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%DOMAIN_NUMBER( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)=domain_idx
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%LOCAL_TYPE( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)= &
                        & DOMAIN_LOCAL_GHOST
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(ghost_node)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          LOCAL_DOF_NUMBERS(domain_idx)=LOCAL_DOF_NUMBERS(domain_idx)+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS= &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                            & LOCAL_DOF_NUMBERS(domain_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)=domain_idx
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                            & DOMAIN_LOCAL_GHOST
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ENDDO !no_ghost_node
                    DEALLOCATE(GHOST_NODES)
                  ENDDO !domain_idx
                  
                  !Check decomposition and check that each domain has a node in it.
                  ALLOCATE(NODE_COUNT(0:number_computational_nodes-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node count.",ERR,ERROR,*999)
                  NODE_COUNT=0
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    no_computational_node=DOMAIN%NODE_DOMAIN(node_idx)
                    IF(no_computational_node>=0.AND.no_computational_node<number_computational_nodes) THEN
                      NODE_COUNT(no_computational_node)=NODE_COUNT(no_computational_node)+1
                    ELSE
                      LOCAL_ERROR="The computational node number of "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))// &
                        & " for node number "//TRIM(NUMBER_TO_VSTRING(node_idx,"*",ERR,ERROR))// &
                        & " is invalid. The computational node number must be between 0 and "// &
                        & TRIM(NUMBER_TO_VSTRING(number_computational_nodes-1,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !node_idx
                  DO no_computational_node=0,number_computational_nodes-1
                    IF(NODE_COUNT(no_computational_node)==0) THEN
                      LOCAL_ERROR="Invalid decomposition. There are no nodes in computational node "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_computational_node
                  DEALLOCATE(NODE_COUNT)
          
                  DEALLOCATE(GHOST_NODES_LIST)
                  DEALLOCATE(LOCAL_NODE_NUMBERS)
                  
                  !Calculate node and dof local to global maps from global to local map
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(NODES_MAPPING,ERR,ERROR,*999)
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOFS_MAPPING,ERR,ERROR,*999)
                  
                ELSE
                  CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Domain mappings dofs is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain mappings nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node decomposition :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Domain = ",DOMAIN%NODE_DOMAIN(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
            & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !node_idx     
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO node_idx=1,NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global node = ", &
          & NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal nodes = ", &
          & NODES_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%INTERNAL_START:NODES_MAPPING%INTERNAL_FINISH), &
          & '("    Internal nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary nodes = ", &
          & NODES_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%BOUNDARY_START:NODES_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost nodes = ", &
          & NODES_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_GHOST,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%GHOST_START:NODES_MAPPING%GHOST_FINISH), &
          & '("    Ghost nodes   :",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & NODES_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS_PTR( &
        & NODES_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,NODES_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Dofs mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ny=1,MESH_TOPOLOGY%DOFS%numberOfDofs
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
            & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ny   
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ny=1,DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global dof = ", &
          & DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(ny),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%INTERNAL_START:DOFS_MAPPING%INTERNAL_FINISH), &
          & '("    Internal dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%BOUNDARY_START:DOFS_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%GHOST_START:DOFS_MAPPING%GHOST_FINISH), &
          & '("    Ghost dofs   :",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & DOFS_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS_PTR( &
        & DOFS_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,DOFS_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ALL_DOMAINS)) DEALLOCATE(ALL_DOMAINS)
    IF(ALLOCATED(GHOST_NODES)) DEALLOCATE(GHOST_NODES)
    IF(ALLOCATED(NUMBER_INTERNAL_NODES)) DEALLOCATE(NUMBER_INTERNAL_NODES)
    IF(ALLOCATED(NUMBER_BOUNDARY_NODES)) DEALLOCATE(NUMBER_BOUNDARY_NODES)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*997)
997 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE
 !============================================================================
  !>Finalises the node mapping in the given domain mappings. \todo pass in the nodes mapping
  SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mapping to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%NODES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the node mapping in the given domain mapping. \todo finalise on error
  SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL FlagError("Domain nodes mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings nodes.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%NODES,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the domain topology.
  SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne,np
    TYPE(BASIS_TYPE), POINTER :: BASIS

    ENTERS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Find maximum number of element parameters for all elements in the domain toplogy
      TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=-1
      DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        IF(ASSOCIATED(BASIS)) THEN
          IF(BASIS%NUMBER_OF_ELEMENT_PARAMETERS>TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS) &
            & TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
        ELSE
          CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDDO !ne
      !Find maximum number of derivatives for all nodes in the domain toplogy
      TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=-1
      DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
        IF(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES>TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES) &
            & TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
      ENDDO !np
      !Calculate the elements surrounding the nodes in the domain topology
      CALL DomainTopology_NodesSurroundingElementsCalculate(TOPOLOGY,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_CALCULATE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Initialises the local domain topology from the mesh topology.
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to initialise the domain topology from the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_element,global_element,local_node,global_node,version_idx,derivative_idx,node_idx,dof_idx, &
      & component_idx
    INTEGER(INTG) :: ne,nn,nkk,INSERT_STATUS
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(MeshNodesType), POINTER :: MESH_NODES
    TYPE(DOMAIN_DOFS_TYPE), POINTER :: DOMAIN_DOFS

    ENTERS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
          IF(ASSOCIATED(DOMAIN%MESH)) THEN
            MESH=>DOMAIN%MESH
            component_idx=DOMAIN%MESH_COMPONENT_NUMBER
            IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
              MESH_ELEMENTS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS
              DOMAIN_ELEMENTS=>DOMAIN%TOPOLOGY%ELEMENTS
              MESH_NODES=>MESH%TOPOLOGY(component_idx)%PTR%NODES
              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
              DOMAIN_DOFS=>DOMAIN%TOPOLOGY%DOFS
              ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(DOMAIN%MAPPINGS%ELEMENTS%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain elements elements.",ERR,ERROR,*999)
              DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_LOCAL
              DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_GLOBAL
              ALLOCATE(DOMAIN_NODES%NODES(DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain nodes nodes.",ERR,ERROR,*999)
              DOMAIN_NODES%NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%NUMBER_OF_LOCAL
              DOMAIN_NODES%TOTAL_NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_NODES%NUMBER_OF_GLOBAL_NODES=DOMAIN%MAPPINGS%NODES%NUMBER_OF_GLOBAL
              ALLOCATE(DOMAIN_DOFS%DOF_INDEX(3,DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain dofs dof index.",ERR,ERROR,*999)
              DOMAIN_DOFS%NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%NUMBER_OF_LOCAL
              DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_DOFS%NUMBER_OF_GLOBAL_DOFS=DOMAIN%MAPPINGS%DOFS%NUMBER_OF_GLOBAL
              !Loop over the domain nodes and calculate the parameters from the mesh nodes
              CALL TREE_CREATE_START(DOMAIN_NODES%NODES_TREE,ERR,ERROR,*999)
              CALL TREE_INSERT_TYPE_SET(DOMAIN_NODES%NODES_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
              CALL TREE_CREATE_FINISH(DOMAIN_NODES%NODES_TREE,ERR,ERROR,*999)
              dof_idx=0
              DO local_node=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                CALL DOMAIN_TOPOLOGY_NODE_INITIALISE(DOMAIN_NODES%NODES(local_node),ERR,ERROR,*999)
                global_node=DOMAIN%MAPPINGS%NODES%LOCAL_TO_GLOBAL_MAP(local_node)
                DOMAIN_NODES%NODES(local_node)%LOCAL_NUMBER=local_node
                DOMAIN_NODES%NODES(local_node)%MESH_NUMBER=global_node
                DOMAIN_NODES%NODES(local_node)%GLOBAL_NUMBER=MESH_NODES%NODES(global_node)%globalNumber
                DOMAIN_NODES%NODES(local_node)%USER_NUMBER=MESH_NODES%NODES(global_node)%userNumber
                CALL TREE_ITEM_INSERT(DOMAIN_NODES%NODES_TREE,DOMAIN_NODES%NODES(local_node)%USER_NUMBER,local_node, &
                  & INSERT_STATUS,ERR,ERROR,*999)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_SURROUNDING_ELEMENTS=0
                NULLIFY(DOMAIN_NODES%NODES(local_node)%SURROUNDING_ELEMENTS)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES=MESH_NODES%NODES(global_node)%numberOfDerivatives
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(MESH_NODES%NODES(global_node)%numberOfDerivatives),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain node derivatives.",ERR,ERROR,*999)
                DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                  CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE( &
                    & DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%globalDerivativeIndex
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX= &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%partialDerivativeIndex
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%numberOfVersions= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                  ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%userVersionNumbers( &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node derivative version numbers.",ERR,ERROR,*999)
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%userVersionNumbers(1: &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions)= &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%userVersionNumbers(1: &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions)
                  ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%DOF_INDEX( &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node dervative versions dof index.",ERR,ERROR,*999)
                  DO version_idx=1,DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                    dof_idx=dof_idx+1
                    DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)=dof_idx
                    DOMAIN_DOFS%DOF_INDEX(1,dof_idx)=version_idx
                    DOMAIN_DOFS%DOF_INDEX(2,dof_idx)=derivative_idx
                    DOMAIN_DOFS%DOF_INDEX(3,dof_idx)=local_node
                  ENDDO !version_idx
                ENDDO !derivative_idx
                DOMAIN_NODES%NODES(local_node)%BOUNDARY_NODE=MESH_NODES%NODES(global_node)%boundaryNode
              ENDDO !local_node
              !Loop over the domain elements and renumber from the mesh elements
              DO local_element=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                CALL DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(DOMAIN_ELEMENTS%ELEMENTS(local_element),ERR,ERROR,*999)
                global_element=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(local_element)               
                BASIS=>MESH_ELEMENTS%ELEMENTS(global_element)%BASIS
                DOMAIN_ELEMENTS%ELEMENTS(local_element)%NUMBER=local_element
                DOMAIN_ELEMENTS%ELEMENTS(local_element)%BASIS=>BASIS
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element nodes.",ERR,ERROR,*999)
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                  & BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element derivatives.",ERR,ERROR,*999)
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementVersions(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                  & BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element versions.",ERR,ERROR,*999)
                DO nn=1,BASIS%NUMBER_OF_NODES
                  global_node=MESH_ELEMENTS%ELEMENTS(global_element)%MESH_ELEMENT_NODES(nn)
                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(global_node)%LOCAL_NUMBER(1)
                  DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_NODES(nn)=local_node
                  DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                    !Find equivalent node derivative by matching partial derivative index
                    !/todo Take a look at this - is it needed any more?
                    FOUND=.FALSE.
                    DO nkk=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                      IF(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(nkk)%PARTIAL_DERIVATIVE_INDEX == &
                        & BASIS%PARTIAL_DERIVATIVE_INDEX(derivative_idx,nn)) THEN
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !nkk
                    IF(FOUND) THEN
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(derivative_idx,nn)=nkk
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementVersions(derivative_idx,nn) = & 
                        & MESH_ELEMENTS%ELEMENTS(global_element)%USER_ELEMENT_NODE_VERSIONS(derivative_idx,nn)
                    ELSE
                      CALL FlagError("Could not find equivalent node derivative",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !derivative_idx
                ENDDO !nn
              ENDDO !local_element                       
            ELSE
              CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh is not associated",ERR,ERROR,*999)

          ENDIF
        ELSE
          CALL FlagError("Domain mapping is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Initialised domain topology :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain nodes = ",DOMAIN_NODES%TOTAL_NUMBER_OF_NODES, &
        & ERR,ERROR,*999)
      DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",DOMAIN_NODES%NODES(node_idx)%LOCAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",DOMAIN_NODES%NODES(node_idx)%MESH_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",DOMAIN_NODES%NODES(node_idx)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",DOMAIN_NODES%NODES(node_idx)%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
          & DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES,ERR,ERROR,*999)
        DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node local derivative number = ",derivative_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Global derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Partial derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions,4,4, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX, &
            & '("        Degree-of-freedom index(version_idx)  :",4(X,I9))','(36X,4(X,I9))',ERR,ERROR,*999)
        ENDDO !derivative_idx
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary node = ", &
          & DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE,ERR,ERROR,*999)
     ENDDO !node_idx
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total Number of domain dofs = ",DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS, &
        & ERR,ERROR,*999)
      DO dof_idx=1,DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dof number = ",dof_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3, &
          & DOMAIN_DOFS%DOF_INDEX(:,dof_idx),'("    Degree-of-freedom index :",3(X,I9))','(29X,3(X,I9))', &
          & ERR,ERROR,*999)
      ENDDO !dof_idx
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain elements = ", &
        & DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",DOMAIN_ELEMENTS%ELEMENTS(ne)%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis user number = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local nodes = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES,'("    Element nodes(nn) :",8(X,I9))','(23X,8(X,I9))', &
          & ERR,ERROR,*999)
        DO nn=1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Local node number : ",nn,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(:,nn), &
            & '("        Element derivatives :",8(X,I2))','(29X,8(X,I2))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%elementVersions(:,nn), &
            & '("        Element versions    :",8(X,I2))','(29X,8(X,I2))',ERR,ERROR,*999)
        ENDDO !nn
      ENDDO !ne
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH
  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs in the given domain topology. \todo pass in the dofs topolgy
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        IF(ALLOCATED(TOPOLOGY%DOFS%DOF_INDEX)) DEALLOCATE(TOPOLOGY%DOFS%DOF_INDEX)
        DEALLOCATE(TOPOLOGY%DOFS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the dofs data structures for a domain topology. \todo finalise on exit
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        CALL FlagError("Decomposition already has topology dofs associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology dofs",ERR,ERROR,*999)
        TOPOLOGY%DOFS%NUMBER_OF_DOFS=0
        TOPOLOGY%DOFS%TOTAL_NUMBER_OF_DOFS=0
        TOPOLOGY%DOFS%NUMBER_OF_GLOBAL_DOFS=0
        TOPOLOGY%DOFS%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ELEMENT_NODES)) DEALLOCATE(ELEMENT%ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%ELEMENT_DERIVATIVES)) DEALLOCATE(ELEMENT%ELEMENT_DERIVATIVES)
    IF(ALLOCATED(ELEMENT%elementVersions)) DEALLOCATE(ELEMENT%elementVersions)
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%NUMBER=0
    NULLIFY(ELEMENT%BASIS)
  
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_INITALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_INITALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given domain topology. \todo pass in the domain elements
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DOMAIN_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Decomposition already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=0
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE
  
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given domain. \todo pass in domain topology
  SUBROUTINE DOMAIN_TOPOLOGY_FINALISE(DOMAIN,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_TOPOLOGY_NODES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_DOFS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_LINES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_FACES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%TOPOLOGY)
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given domain. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !A pointer to the domain to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        CALL FlagError("Domain already has topology associated",ERR,ERROR,*999)
      ELSE
        !Allocate domain topology
        ALLOCATE(DOMAIN%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Domain topology could not be allocated",ERR,ERROR,*999)
        DOMAIN%TOPOLOGY%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%TOPOLOGY%ELEMENTS)
        NULLIFY(DOMAIN%TOPOLOGY%NODES)
        NULLIFY(DOMAIN%TOPOLOGY%DOFS)
        NULLIFY(DOMAIN%TOPOLOGY%LINES)
        NULLIFY(DOMAIN%TOPOLOGY%FACES)
        !Initialise the topology components
        CALL DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_NODES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_DOFS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_LINES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_FACES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        !Initialise the domain topology from the domain mappings and the mesh it came from
        CALL DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*999)
        !Calculate the topological information.
        CALL DOMAIN_TOPOLOGY_CALCULATE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    IF(ALLOCATED(LINE%NODES_IN_LINE)) DEALLOCATE(LINE%NODES_IN_LINE)
    IF(ALLOCATED(LINE%DERIVATIVES_IN_LINE)) DEALLOCATE(LINE%DERIVATIVES_IN_LINE)
 
    EXITS("DOMAIN_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a domain topology line.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    LINE%BOUNDARY_LINE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given domain topology. \todo pass in domain lines
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DOMAIN_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !
  !>Finalises a face in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    IF(ALLOCATED(FACE%NODES_IN_FACE)) DEALLOCATE(FACE%NODES_IN_FACE)
    IF(ALLOCATED(FACE%DERIVATIVES_IN_FACE)) DEALLOCATE(FACE%DERIVATIVES_IN_FACE)
 
    EXITS("DOMAIN_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a domain topology face.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    FACE%BOUNDARY_FACE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given domain topology. \todo pass in domain faces
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DOMAIN_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a domain nodes topolgoy. 
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(DOMAIN_TOPOLOGY,USER_NODE_NUMBER,NODE_EXISTS,DOMAIN_LOCAL_NODE_NUMBER, &
    & GHOST_NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY !<A pointer to the domain topology to check the node exists on
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: NODE_EXISTS !<On exit, is .TRUE. if the node user number exists in the domain nodes topolgoy (even if it is a ghost node), .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_LOCAL_NODE_NUMBER !<On exit, if the node exists the local number corresponding to the user node number. If the node does not exist then global number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_NODE !<On exit, is .TRUE. if the local node (if it exists) is a ghost node, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    ENTERS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    DOMAIN_LOCAL_NODE_NUMBER=0
    GHOST_NODE=.FALSE.
    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
      IF(ASSOCIATED(DOMAIN_NODES)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(DOMAIN_NODES%NODES_TREE,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(DOMAIN_NODES%NODES_TREE,TREE_NODE,DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
          NODE_EXISTS=.TRUE.
          GHOST_NODE=DOMAIN_LOCAL_NODE_NUMBER>DOMAIN_NODES%NUMBER_OF_NODES
        ENDIF
      ELSE
        CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node derivative and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE_DERIVATIVE%userVersionNumbers)) DEALLOCATE(NODE_DERIVATIVE%userVersionNumbers)
    IF(ALLOCATED(NODE_DERIVATIVE%DOF_INDEX)) DEALLOCATE(NODE_DERIVATIVE%DOF_INDEX)

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR,*999)

    NODE_DERIVATIVE%numberOfVersions=0
    NODE_DERIVATIVE%GLOBAL_DERIVATIVE_INDEX=0
    NODE_DERIVATIVE%PARTIAL_DERIVATIVE_INDEX=0

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx

    ENTERS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE%DERIVATIVES)) THEN
      DO derivative_idx=1,NODE%NUMBER_OF_DERIVATIVES
        CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
      ENDDO !derivative_idx
      DEALLOCATE(NODE%DERIVATIVES)
    ENDIF
    IF(ASSOCIATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(NODE%NODE_LINES)) DEALLOCATE(NODE%NODE_LINES)
 
    EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%LOCAL_NUMBER=0
    NODE%MESH_NUMBER=0
    NODE%GLOBAL_NUMBER=0
    NODE%USER_NUMBER=0
    NODE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    NODE%NUMBER_OF_NODE_LINES=0
    NODE%BOUNDARY_NODE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the nodees in the given domain topology. \todo pass in domain nodes
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np

    ENTERS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
          CALL DOMAIN_TOPOLOGY_NODE_FINALISE(TOPOLOGY%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) DEALLOCATE(TOPOLOGY%NODES%NODES)
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES_TREE)) CALL TREE_DESTROY(TOPOLOGY%NODES%NODES_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%NODES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        CALL FlagError("Decomposition already has topology nodes associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology nodes",ERR,ERROR,*999)
        TOPOLOGY%NODES%NUMBER_OF_NODES=0
        TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES=0
        TOPOLOGY%NODES%NUMBER_OF_GLOBAL_NODES=0
        TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=0
        TOPOLOGY%NODES%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%NODES%NODES)
        NULLIFY(TOPOLOGY%NODES%NODES_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a domain.
  SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate the elements surrounding the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_no,insert_position,ne,nn,np,surrounding_elem_no
    INTEGER(INTG), POINTER :: NEW_SURROUNDING_ELEMENTS(:)
    LOGICAL :: FOUND_ELEMENT
    TYPE(BASIS_TYPE), POINTER :: BASIS

    NULLIFY(NEW_SURROUNDING_ELEMENTS)

    ENTERS("DomainTopology_NodesSurroundingElementsCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) THEN
            DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
              TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS=0
              IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)) &
                & DEALLOCATE(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)
            ENDDO !np
            DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
              BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
              DO nn=1,BASIS%NUMBER_OF_NODES
                np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                FOUND_ELEMENT=.FALSE.
                element_no=1
                insert_position=1
                DO WHILE(element_no<=TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS.AND..NOT.FOUND_ELEMENT)
                  surrounding_elem_no=TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(element_no)
                  IF(surrounding_elem_no==ne) THEN
                    FOUND_ELEMENT=.TRUE.
                  ENDIF
                  element_no=element_no+1
                  IF(ne>=surrounding_elem_no) THEN
                    insert_position=element_no
                  ENDIF
                ENDDO
                IF(.NOT.FOUND_ELEMENT) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(NEW_SURROUNDING_ELEMENTS(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate new surrounding elements",ERR,ERROR,*999)
                  IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)) THEN
                    NEW_SURROUNDING_ELEMENTS(1:insert_position-1)=TOPOLOGY%NODES%NODES(np)% &
                      & SURROUNDING_ELEMENTS(1:insert_position-1)
                    NEW_SURROUNDING_ELEMENTS(insert_position)=ne
                    NEW_SURROUNDING_ELEMENTS(insert_position+1:TOPOLOGY%NODES%NODES(np)% &
                      & NUMBER_OF_SURROUNDING_ELEMENTS+1)=TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(insert_position: &
                      & TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS)
                    DEALLOCATE(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)
                  ELSE
                    NEW_SURROUNDING_ELEMENTS(1)=ne
                  ENDIF
                  TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS=>NEW_SURROUNDING_ELEMENTS
                  TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS= &
                    & TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1
                ENDIF
              ENDDO !nn
            ENDDO !ne
          ELSE
            CALL FlagError("Domain topology nodes nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain topology nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN
999 IF(ASSOCIATED(NEW_SURROUNDING_ELEMENTS)) DEALLOCATE(NEW_SURROUNDING_ELEMENTS)
    ERRORS("DomainTopology_NodesSurroundingElementsCalculate",ERR,ERROR)
    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN 1   
  END SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finialises the mesh adjacent elements information and deallocates all memory
  SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ADJACENT_ELEMENT_TYPE) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
    IF(ALLOCATED(MESH_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(MESH_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)
       
    EXITS("MESH_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises the mesh adjacent elements information.
  SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ADJACENT_ELEMENT_TYPE) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
       
    EXITS("MESH_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a mesh. \see OPENCMISS::CMISSMeshCreateFinish
  SUBROUTINE MESH_CREATE_FINISH(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    LOGICAL :: FINISHED
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
        !Check that the mesh component elements have been finished
        FINISHED=.TRUE.
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
            IF(.NOT.MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS_FINISHED) THEN
              LOCAL_ERROR="The elements for mesh component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " have not been finished"
              FINISHED=.FALSE.
              EXIT
            ENDIF
          ELSE
            LOCAL_ERROR="The elements for mesh topology component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " are not associated"
            FINISHED=.FALSE.
            EXIT
          ENDIF
        ENDDO !component_idx
        IF(.NOT.FINISHED) CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        MESH%MESH_FINISHED=.TRUE.
        !Calulcate the mesh topology
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL MeshTopologyCalculate(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
        ENDDO !component_idx
      ELSE
        CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh user number       = ",MESH%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",MESH%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ",MESH%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_CREATE_FINISH",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh 
  SUBROUTINE MESH_CREATE_START_GENERIC(MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<The pointer to the meshes
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,mesh_idx
    TYPE(MESH_TYPE), POINTER :: NEW_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(NEW_MESH)
    NULLIFY(NEW_MESHES)

    ENTERS("MESH_CREATE_START_GENERIC",ERR,ERROR,*997)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*997)
      ELSE
        CALL MESH_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Set default mesh values
        NEW_MESH%USER_NUMBER=USER_NUMBER
        NEW_MESH%GLOBAL_NUMBER=MESHES%NUMBER_OF_MESHES+1
        NEW_MESH%MESHES=>MESHES
        NEW_MESH%NUMBER_OF_DIMENSIONS=NUMBER_OF_DIMENSIONS
        NEW_MESH%NUMBER_OF_COMPONENTS=1
        NEW_MESH%SURROUNDING_ELEMENTS_CALCULATE=.true. !default true
        !Initialise mesh topology and decompositions
        CALL MeshTopologyInitialise(NEW_MESH,ERR,ERROR,*999)
        CALL DECOMPOSITIONS_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Add new mesh into list of meshes 
        ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES+1),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
        DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
          NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
        ENDDO !mesh_idx
        NEW_MESHES(MESHES%NUMBER_OF_MESHES+1)%PTR=>NEW_MESH
        IF(ASSOCIATED(MESHES%MESHES)) DEALLOCATE(MESHES%MESHES)
        MESHES%MESHES=>NEW_MESHES
        MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES+1
        MESH=>NEW_MESH
      ENDIF
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*997)
    ENDIF
      
    EXITS("MESH_CREATE_START_GENERIC")
    RETURN
999 CALL MESH_FINALISE(NEW_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    NULLIFY(MESH)    
997 ERRORSEXITS("MESH_CREATE_START_GENERIC",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_GENERIC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in an interface. \see OPENCMISS::CMISSMeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_INTERFACE(USER_NUMBER,INTERFACE,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(INTERFACE%MESHES)) THEN
          CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on interface number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(INTERFACE%INTERFACES)) THEN
              PARENT_REGION=>INTERFACE%INTERFACES%PARENT_REGION
              IF(ASSOCIATED(PARENT_REGION)) THEN
                IF(ASSOCIATED(PARENT_REGION%COORDINATE_SYSTEM)) THEN                  
                  IF(NUMBER_OF_DIMENSIONS>0) THEN
                    IF(NUMBER_OF_DIMENSIONS<=PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                      CALL MESH_CREATE_START_GENERIC(INTERFACE%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                      MESH%INTERFACE=>INTERFACE
                    ELSE
                      LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & ") must be <= number of parent region dimensions ("// &
                        & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Parent region coordinate system is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interfaces parent region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface interfaces is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_START_INTERFACE")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_INTERFACE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in the region identified by REGION. \see OPENCMISS::CMISSMeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_REGION(USER_NUMBER,REGION,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(REGION%MESHES)) THEN
          CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
              IF(NUMBER_OF_DIMENSIONS>0) THEN
                IF(NUMBER_OF_DIMENSIONS<=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                  CALL MESH_CREATE_START_GENERIC(REGION%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                  MESH%REGION=>REGION
                ELSE
                  LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                    & ") must be <= number of region dimensions ("// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The coordinate system on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
                & " are not associated."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_CREATE_START_REGION")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_REGION",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_REGION

  !
  !================================================================================================================================
  !

  !>Destroys the mesh identified by a user number on the given region and deallocates all memory. \see OPENCMISS::CMISSMeshDestroy
  SUBROUTINE MESH_DESTROY_NUMBER(USER_NUMBER,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to destroy
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN

!!TODO: have a mesh_destory_ptr and mesh_destroy_number
        
        !Find the problem identified by the user number
        FOUND=.FALSE.
        mesh_position=0
        DO WHILE(mesh_position<REGION%MESHES%NUMBER_OF_MESHES.AND..NOT.FOUND)
          mesh_position=mesh_position+1
          IF(REGION%MESHES%MESHES(mesh_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          MESH=>REGION%MESHES%MESHES(mesh_position)%PTR

          CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

          !Remove the mesh from the list of meshes
          IF(REGION%MESHES%NUMBER_OF_MESHES>1) THEN
            ALLOCATE(NEW_MESHES(REGION%MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
            DO mesh_idx=1,REGION%MESHES%NUMBER_OF_MESHES
              IF(mesh_idx<mesh_position) THEN
                NEW_MESHES(mesh_idx)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ELSE IF(mesh_idx>mesh_position) THEN
                REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
                NEW_MESHES(mesh_idx-1)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ENDIF
            ENDDO !mesh_idx
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%MESHES=>NEW_MESHES
            REGION%MESHES%NUMBER_OF_MESHES=REGION%MESHES%NUMBER_OF_MESHES-1
          ELSE
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%NUMBER_OF_MESHES=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys the mesh and deallocates all memory. \see OPENCMISS::CMISSMeshDestroy
  SUBROUTINE MESH_DESTROY(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    TYPE(MESHES_TYPE), POINTER :: MESHES
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      MESHES=>MESH%MESHES
      IF(ASSOCIATED(MESHES)) THEN
        mesh_position=MESH%GLOBAL_NUMBER
          
        CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

        !Remove the mesh from the list of meshes
        IF(MESHES%NUMBER_OF_MESHES>1) THEN
          ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new meshes.",ERR,ERROR,*999)
          DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
            IF(mesh_idx<mesh_position) THEN
              NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ELSE IF(mesh_idx>mesh_position) THEN
              MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
              NEW_MESHES(mesh_idx-1)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ENDIF
          ENDDO !mesh_idx
          DEALLOCATE(MESHES%MESHES)
          MESHES%MESHES=>NEW_MESHES
          MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES-1
        ELSE
          DEALLOCATE(MESHES%MESHES)
          MESHES%NUMBER_OF_MESHES=0
        ENDIF
      ELSE
        CALL FlagError("The mesh meshes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a mesh and deallocates all memory.
  SUBROUTINE MESH_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL MeshTopologyFinalise(MESH,ERR,ERROR,*999)
      CALL DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*999)
!      IF(ASSOCIATED(MESH%INTF)) CALL INTERFACE_MESH_FINALISE(MESH,ERR,ERROR,*999)  ! <<??>>
      DEALLOCATE(MESH)
    ENDIF
 
    EXITS("MESH_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a region nodes pointer corresponding to the mesh global nodes accounting for interfaces.
  SUBROUTINE MeshGlobalNodesGet(mesh,nodes,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the global nodes for
    TYPE(NODES_TYPE), POINTER :: nodes !<On return, the nodes pointer corresponding to the global nodes for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshGlobalNodesGet",err,error,*999)
    
    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(nodes)) THEN
        CALL FlagError("Nodes is already associated.",err,error,*999)
      ELSE
        NULLIFY(nodes)
        region=>mesh%region
        IF(ASSOCIATED(region)) THEN
          nodes=>region%nodes
        ELSE
          INTERFACE=>mesh%INTERFACE
          IF(ASSOCIATED(interface)) THEN
            nodes=>interface%nodes
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have an associated region or interface."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
        IF(.NOT.ASSOCIATED(nodes)) THEN
          IF(ASSOCIATED(region)) THEN
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have any nodes associated with the mesh region."
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have any nodes associated with the mesh interface."
          ENDIF
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
        
    EXITS("MeshGlobalNodesGet")
    RETURN
999 ERRORSEXITS("MeshGlobalNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshGlobalNodesGet
  
  !
  !================================================================================================================================
  !

  !>Initialises a mesh.
  SUBROUTINE MESH_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(MESH,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate new mesh.",ERR,ERROR,*999)
      MESH%USER_NUMBER=0
      MESH%GLOBAL_NUMBER=0
      MESH%MESH_FINISHED=.FALSE.
      NULLIFY(MESH%MESHES)
      NULLIFY(MESH%REGION)
      NULLIFY(MESH%INTERFACE)
      NULLIFY(MESH%GENERATED_MESH)
      MESH%NUMBER_OF_DIMENSIONS=0
      MESH%NUMBER_OF_COMPONENTS=0
      MESH%MESH_EMBEDDED=.FALSE.
      NULLIFY(MESH%EMBEDDING_MESH)
      MESH%NUMBER_OF_EMBEDDED_MESHES=0
      NULLIFY(MESH%EMBEDDED_MESHES)
      MESH%NUMBER_OF_ELEMENTS=0
      NULLIFY(MESH%TOPOLOGY)
      NULLIFY(MESH%DECOMPOSITIONS)
    ENDIF
    
    EXITS("MESH_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the number of mesh components for a mesh identified by a pointer. \see OPENCMISS::CMISSMeshNumberOfComponentsGet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of components for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COMPONENTS !<On return, the number of components in the specified mesh.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_COMPONENTS=MESH%NUMBER_OF_COMPONENTS
      ELSE
        CALL FlagError("Mesh has not finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of mesh components for a mesh. \see OPENCMISS::CMISSMeshNumberOfComponentsSet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of components for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MeshComponentTopologyPtrType), POINTER :: NEW_TOPOLOGY(:)

    NULLIFY(NEW_TOPOLOGY)
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_COMPONENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FlagError("Mesh has been finished",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_COMPONENTS/=MESH%NUMBER_OF_COMPONENTS) THEN
            ALLOCATE(NEW_TOPOLOGY(NUMBER_OF_COMPONENTS),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new topology",ERR,ERROR,*999)
            IF(NUMBER_OF_COMPONENTS<MESH%NUMBER_OF_COMPONENTS) THEN
              DO component_idx=1,NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
            ELSE !NUMBER_OF_COMPONENTS>MESH%NUMBER_OF_COMPONENTS
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
!!TODO \todo sort out mesh_topology initialise/finalise so that they allocate and deal with this below then call that routine
              DO component_idx=MESH%NUMBER_OF_COMPONENTS+1,NUMBER_OF_COMPONENTS
                ALLOCATE(NEW_TOPOLOGY(component_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate new topology component",ERR,ERROR,*999)
                NEW_TOPOLOGY(component_idx)%PTR%mesh=>mesh
                NEW_TOPOLOGY(component_idx)%PTR%meshComponentNumber=component_idx
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%elements)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%nodes)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%dofs)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%dataPoints)
                !Initialise the topology components
                CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MeshTopologyNodesInitialise(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MeshTopologyDofsInitialise(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
              ENDDO !component_idx
            ENDIF
            IF(ASSOCIATED(MESH%TOPOLOGY)) DEALLOCATE(MESH%TOPOLOGY)
            MESH%TOPOLOGY=>NEW_TOPOLOGY
            MESH%NUMBER_OF_COMPONENTS=NUMBER_OF_COMPONENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of mesh components ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
          & ") is illegal. You must have >0 mesh components"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF     
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_SET")
    RETURN
!!TODO: tidy up memory deallocation on error
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET

  !
  !================================================================================================================================
  !
  
  !>Gets the number of elements for a mesh identified by a pointer. \see OPENCMISS::CMISSMeshNumberOfElementsGet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_ELEMENTS !<On return, the number of elements in the specified mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS
      ELSE
        CALL FlagError("Mesh has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of elements for a mesh. \see OPENCMISS::CMISSMeshNumberOfElementsSet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of elements for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS !<The number of elements to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_ELEMENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_ELEMENTS/=MESH%NUMBER_OF_ELEMENTS) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
                  IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
                    IF(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS>0) THEN
!!TODO: Reallocate the elements and copy information. 
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  CALL FlagError("Mesh topology component pointer is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
            ENDIF
            MESH%NUMBER_OF_ELEMENTS=NUMBER_OF_ELEMENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of elements ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ELEMENTS,"*",ERR,ERROR))// &
          & ") is invalid. You must have > 0 elements."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF     
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET

  !
  !================================================================================================================================
  !

  !>Returns the region for a mesh accounting for regions and interfaces
  SUBROUTINE MeshRegionGet(mesh,region,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the region for
    TYPE(REGION_TYPE), POINTER :: region !<On return, the meshes region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: parentRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshRegionGet",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(region)) THEN
        CALL FlagError("Region is already associated.",err,error,*999)
      ELSE
        NULLIFY(region)
        NULLIFY(interface)
        region=>mesh%region
        IF(.NOT.ASSOCIATED(region)) THEN
          interface=>mesh%interface
          IF(ASSOCIATED(interface)) THEN
            parentRegion=>interface%PARENT_REGION
            IF(ASSOCIATED(parentRegion)) THEN
              region=>parentRegion
            ELSE
              localError="The parent region not associated for mesh number "// &
                & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//" of interface number "// &
                & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The region or interface is not associated for mesh number "// &
              & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF

    EXITS("MeshRegionGet")
    RETURN
999 ERRORSEXITS("MeshRegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshRegionGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the surrounding elements calculate flag. \see OPENCMISS::CMISSMeshSurroundingElementsCalculateSet
  SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET(MESH,SURROUNDING_ELEMENTS_CALCULATE_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the surrounding elements calculate flag for
    LOGICAL, INTENT(IN) :: SURROUNDING_ELEMENTS_CALCULATE_FLAG !<The surrounding elements calculate flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
      ELSE
        MESH%SURROUNDING_ELEMENTS_CALCULATE=SURROUNDING_ELEMENTS_CALCULATE_FLAG
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET")
    RETURN
999 ERRORSEXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology.
  SUBROUTINE MeshTopologyCalculate(topology,err,error,*)    

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      !Calculate the nodes used in the mesh
      CALL MeshTopologyNodesCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the nodes in a mesh
      CALL MeshTopologySurroundingElementsCalculate(topology,err,error,*999)
      !Calculate the number of derivatives at each node in a mesh
      CALL MeshTopologyNodesDerivativesCalculate(topology,err,error,*999)
      !Calculate the number of versions for each derivative at each node in a mesh
      CALL MeshTopologyNodesVersionCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the elements in the mesh
      CALL MeshTopology_ElementsAdjacentElementsCalculate(topology,err,error,*999)
      !Calculate the boundary nodes and elements in the mesh
      CALL MeshTopologyBoundaryCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the elements in the mesh
      CALL MeshTopologyDofsCalculate(topology,err,error,*999)
    ELSE
      CALL FlagError("Topology is not associated",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyCalculate")
    RETURN
999 ERRORSEXITS("MeshTopologyCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyCalculate
  
  !
  !===============================================================================================================================
  !

  !>Calculates the boundary nodes and elements for a mesh topology. 
  SUBROUTINE MeshTopologyBoundaryCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the boundary for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,localNodeIdx,matchIndex,nodeIdx,xiCoordIdx,xiDirection
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopologyBoundaryCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      nodes=>topology%nodes
      IF(ASSOCIATED(nodes)) THEN
        elements=>topology%elements
        IF(ASSOCIATED(elements)) THEN
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            basis=>elements%elements(elementIdx)%basis
            SELECT CASE(basis%type)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              DO xiCoordIdx=-basis%NUMBER_OF_XI_COORDINATES,basis%NUMBER_OF_XI_COORDINATES
                IF(xiCoordIdx/=0) THEN
                  IF(elements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)%NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                    elements%elements(elementIdx)%BOUNDARY_ELEMENT=.TRUE.
                    IF(xiCoordIdx<0) THEN
                      xiDirection=-xiCoordIdx
                      matchIndex=1
                    ELSE
                      xiDirection=xiCoordIdx
                      matchIndex=BASIS%NUMBER_OF_NODES_XIC(xiCoordIdx)
                    END IF                   
                    DO localNodeIdx=1,BASIS%NUMBER_OF_NODES
                      IF(basis%NODE_POSITION_INDEX(localNodeIdx,XIDIRECTION)==matchIndex) THEN
                        nodeIdx=elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)
                        nodes%nodes(nodeIdx)%boundaryNode=.TRUE.
                      ENDIF
                    ENDDO !nn
                  ENDIF
                ENDIF
              ENDDO !xiCoordIdx            
            CASE(BASIS_SIMPLEX_TYPE)
              elements%elements(elementIdx)%BOUNDARY_ELEMENT=.FALSE.
              DO xiCoordIdx=1,basis%NUMBER_OF_XI_COORDINATES
                elements%elements(elementIdx)%BOUNDARY_ELEMENT=elements%elements(elementIdx)%BOUNDARY_ELEMENT.OR. &
                  & elements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)%NUMBER_OF_ADJACENT_ELEMENTS==0
                IF(elements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)%NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                  DO localNodeIdx=1,basis%NUMBER_OF_NODES
                    IF(basis%NODE_POSITION_INDEX(localNodeIdx,xiCoordIdx)==1) THEN
                      nodeIdx=elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)
                      nodes%nodes(nodeIdx)%boundaryNode=.TRUE.
                    ENDIF
                  ENDDO !localNodeIdx                  
                ENDIF
              ENDDO !xiCoordIdx
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !elementIdx
        ELSE
          CALL FlagError("Topology elements is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology nodes is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary elements:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",elements%NUMBER_OF_ELEMENTS,err,error,*999)
      DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Element : ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary element = ",elements%elements(elementIdx)%BOUNDARY_ELEMENT, &
          & err,error,*999)        
      ENDDO !elementIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary nodes:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Node : ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary node = ",nodes%nodes(nodeIdx)%boundaryNode,err,error,*999)        
      ENDDO !elementIdx            
    ENDIF
 
    EXITS("MeshTopologyBoundaryCalculate")
    RETURN
999 ERRORSEXITS("MeshTopologyBoundaryCalculate",err,error)
    RETURN 1
  END SUBROUTINE MeshTopologyBoundaryCalculate

  !
  !===============================================================================================================================
  !

  !>Calculates the degrees-of-freedom for a mesh topology. 
  SUBROUTINE MeshTopologyDofsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,nodeIdx,numberOfDofs,versionIdx
    TYPE(MeshDofsType), POINTER :: dofs
    TYPE(MeshNodesType), POINTER :: nodes

    ENTERS("MeshTopologyDofsCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      nodes=>topology%nodes
      IF(ASSOCIATED(nodes)) THEN
        dofs=>topology%dofs
        IF(ASSOCIATED(dofs)) THEN
          numberOfDofs=0
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex( &
                & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate mesh topology node derivative version dof index.",err,error,*999)
              DO versionIdx=1,nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                numberOfDofs=numberOfDofs+1
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)=numberOfDofs
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
          dofs%numberOfDofs=numberOfDofs
        ELSE
          CALL FlagError("Topology dofs is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology nodes is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF
 
    EXITS("MeshTopologyDofsCalculate")
    RETURN
999 ERRORSEXITS("MeshTopologyDofsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyDofsCalculate

  !
  !===============================================================================================================================
  !

  !>Finalises the dof data structures for a mesh topology and deallocates any memory. \todo pass in dofs
  SUBROUTINE MeshTopologyDofsFinalise(dofs,err,error,*)

    !Argument variables
    TYPE(MeshDofsType), POINTER :: dofs !<A pointer to the mesh topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyDofsFinalise",err,error,*999)

    IF(ASSOCIATED(dofs)) THEN
      DEALLOCATE(dofs)
    ENDIF
 
    EXITS("MeshTopologyDofsFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyDofsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyDofsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the dofs in a given mesh topology.
  SUBROUTINE MeshTopologyDofsInitialise(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("MeshTopologyDofsInitialise",err,error,*998)

    IF(ASSOCIATED(topology)) THEN
      IF(ASSOCIATED(topology%dofs)) THEN
        CALL FlagError("Mesh already has topology dofs associated",err,error,*998)
      ELSE
        ALLOCATE(topology%dofs,STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate topology dofs",err,error,*999)
        topology%dofs%numberOfDofs=0
        topology%dofs%meshComponentTopology=>topology
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",err,error,*998)
    ENDIF
    
    EXITS("MeshTopologyDofsInitialise")
    RETURN
999 CALL MeshTopologyDofsFinalise(topology%dofs,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopologyDofsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyDofsInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating elements for a specified mesh component in a mesh topology. \see OPENCMISS::CMISSMeshElementsCreateFinish
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology

    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Mesh elements have already been finished.",ERR,ERROR,*999)
      ELSE        
        ELEMENTS%ELEMENTS_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      meshComponentTopology=>elements%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(MESH)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of global elements = ",MESH%NUMBER_OF_ELEMENTS, &
            & ERR,ERROR,*999)
          DO ne=1,MESH%NUMBER_OF_ELEMENTS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element = ",ne,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER, &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ",ELEMENTS%ELEMENTS(ne)%USER_NUMBER, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(ELEMENTS%ELEMENTS(ne)%BASIS)) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis number         = ",ELEMENTS%ELEMENTS(ne)%BASIS% &
                & USER_NUMBER,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES)) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)% BASIS%NUMBER_OF_NODES,8,8, &
                & ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES,'("    User element nodes   =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("User element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES)) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
                & ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES,'("    Global element nodes =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("Global element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !ne
        ELSE
          CALL FlagError("Mesh component topology mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh elements mesh component topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating elements in the mesh component identified by MESH and component_idx. The elements will be created with a default basis of BASIS. ELEMENTS is the returned pointer to the MESH_ELEMENTS data structure. \see OPENCMISS::CMISSMeshElementsCreateStart
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,MESH_COMPONENT_NUMBER,BASIS,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to start creating the elements on
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the default basis to use
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the created mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,INSERT_STATUS,ne
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(MESH)) THEN     
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FlagError("Elements is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
              IF(ASSOCIATED(ELEMENTS%ELEMENTS)) THEN
                CALL FlagError("Mesh topology already has elements associated",ERR,ERROR,*998)
              ELSE
                IF(ASSOCIATED(BASIS)) THEN
                  MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%meshComponentNumber=MESH_COMPONENT_NUMBER
                  ALLOCATE(ELEMENTS%ELEMENTS(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate individual elements",ERR,ERROR,*999)
                  ELEMENTS%NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS !Psuedo inheritance of the number of elements
                  CALL TREE_CREATE_START(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  CALL TREE_INSERT_TYPE_SET(ELEMENTS%ELEMENTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                  CALL TREE_CREATE_FINISH(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  ELEMENTS%ELEMENTS_FINISHED=.FALSE.
                  !Set up the default values and allocate element structures
                  DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                    CALL MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER=ne
                    ELEMENTS%ELEMENTS(ne)%USER_NUMBER=ne
                    CALL TREE_ITEM_INSERT(ELEMENTS%ELEMENTS_TREE,ne,ne,INSERT_STATUS,ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%BASIS=>BASIS
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate user element nodes",ERR,ERROR,*999)
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES=1
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES=1
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                      & BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes versions",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS = 1
                  ENDDO !ne
                ELSE
                  CALL FlagError("Basis is not associated",ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Mesh topology elements is not associated",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Mesh topology is not associated",ERR,ERROR,*998)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START")
    RETURN
999 CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,DUMMY_ERR,DUMMY_ERROR,*998)
998 NULLIFY(ELEMENTS)
    ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys the elements in a mesh topology. \todo as this is a user routine it should take a mesh pointer like create start and finish? Split this into destroy and finalise?
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to destroy 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic
    
    ENTERS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODE_VERSIONS)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODE_VERSIONS)
    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODES)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%GLOBAL_ELEMENT_NODES)) DEALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%MESH_ELEMENT_NODES)) DEALLOCATE(ELEMENT%MESH_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) THEN
      DO nic=LBOUND(ELEMENT%ADJACENT_ELEMENTS,1),UBOUND(ELEMENT%ADJACENT_ELEMENTS,1)
        CALL MESH_ADJACENT_ELEMENT_FINALISE(ELEMENT%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh elements for a given mesh component.
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_GET(MESH,MESH_COMPONENT_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the elements for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to get the elements for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("MESH_TOPOLOGY_ELEMENTS_GET",ERR,ERROR,*998)
    
    IF(ASSOCIATED(MESH)) THEN
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FlagError("Elements is already associated.",ERR,ERROR,*998)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
            ELSE
              CALL FlagError("Mesh topology elements is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_GET")
    RETURN
999 NULLIFY(ELEMENTS)
998 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_GET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    NULLIFY(ELEMENT%BASIS)
    ELEMENT%BOUNDARY_ELEMENT=.FALSE.
    
    EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: Take user number
  
  !>Gets the basis for a mesh element identified by a given global number. \todo should take user number \see OPENCMISS::CMISSMeshElementsBasisGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to get the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to get the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<On return, a pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
          BASIS=>ELEMENT%BASIS
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the basis for a mesh element identified by a given global number. \todo should take user number
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: NEW_USER_ELEMENT_NODES(:),NEW_GLOBAL_ELEMENT_NODES(:),NEW_USER_ELEMENT_NODE_VERSIONS(:,:)
    INTEGER(INTG) :: OVERLAPPING_NUMBER_NODES,OVERLAPPING_NUMBER_DERIVATIVES
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(ASSOCIATED(BASIS)) THEN
            ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
            IF(ELEMENT%BASIS%NUMBER_OF_NODES/=BASIS%NUMBER_OF_NODES.OR. &
                & ELEMENT%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES/=BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES) THEN
              !Allocate new user and global element nodes
              ALLOCATE(NEW_USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                & BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element node versions",ERR,ERROR,*999)

              OVERLAPPING_NUMBER_NODES=MIN(BASIS%NUMBER_OF_NODES,ELEMENT%BASIS%NUMBER_OF_NODES)
              OVERLAPPING_NUMBER_DERIVATIVES=MIN(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,ELEMENT%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES)

              !Set default values
              NEW_USER_ELEMENT_NODE_VERSIONS=1
              NEW_USER_ELEMENT_NODES(OVERLAPPING_NUMBER_NODES+1:)=0
              NEW_GLOBAL_ELEMENT_NODES(OVERLAPPING_NUMBER_NODES+1:)=0
              !Copy previous values
              NEW_USER_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%USER_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)
              NEW_GLOBAL_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%GLOBAL_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)
              NEW_USER_ELEMENT_NODE_VERSIONS(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)= &
                & ELEMENT%USER_ELEMENT_NODE_VERSIONS(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)

              !Replace arrays with new ones
              CALL MOVE_ALLOC(NEW_USER_ELEMENT_NODE_VERSIONS,ELEMENT%USER_ELEMENT_NODE_VERSIONS)
              CALL MOVE_ALLOC(NEW_USER_ELEMENT_NODES,ELEMENT%USER_ELEMENT_NODES)
              CALL MOVE_ALLOC(NEW_GLOBAL_ELEMENT_NODES,ELEMENT%GLOBAL_ELEMENT_NODES)
            END IF           
            ELEMENT%BASIS=>BASIS
          ELSE
            CALL FlagError("Basis is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the adjacent element number for a mesh element identified by a global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNo
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET(GLOBAL_NUMBER,ELEMENTS,ADJACENT_ELEMENT_XI,ADJACENT_ELEMENT_NUMBER, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to get the adjacent element for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements of a mesh component from which to get the adjacent element from.
    INTEGER(INTG), INTENT(IN) :: ADJACENT_ELEMENT_XI !< The xi coordinate direction to get the adjacent element Note that -xiCoordinateDirection gives the adjacent element before the element in the xiCoordinateDirection'th direction and +xiCoordinateDirection gives the adjacent element after the element in the xiCoordinateDirection'th direction. The xiCoordinateDirection=0 index will give the information on the current element.
    INTEGER(INTG), INTENT(OUT) :: ADJACENT_ELEMENT_NUMBER !<On return, the adjacent element number in the specified xi coordinate direction. Return 0 if the specified element has no adjacent elements in the specified xi coordinate direction.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(ADJACENT_ELEMENT_XI>=-ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI .AND. &
            & ADJACENT_ELEMENT_XI<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI) THEN
            IF(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%ADJACENT_ELEMENTS(ADJACENT_ELEMENT_XI)%NUMBER_OF_ADJACENT_ELEMENTS > 0) THEN !\todo Currently returns only the first adjacent element for now as the python binding require the output array size of the adjacent element to be known a-prior. Add routine to first output number of adjacent elements and then loop over all adjacent elements
              ADJACENT_ELEMENT_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%ADJACENT_ELEMENTS(ADJACENT_ELEMENT_XI)%ADJACENT_ELEMENTS(1)
            ELSE !Return 0 indicating the specified element has no adjacent elements in the specified xi coordinate direction.
              ADJACENT_ELEMENT_NUMBER=0
            ENDIF
          ELSE
            LOCAL_ERROR="The specified adjacent element xi is invalid. The supplied xi is "// &
            & TRIM(NUMBER_TO_VSTRING(ADJACENT_ELEMENT_XI,"*",ERR,ERROR))//" and needs to be >=-"// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//" and <="// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET",ERR,ERROR)
    RETURN 1

  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET

  !
  !================================================================================================================================
  !

  !>Gets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(OUT) :: USER_ELEMENT_NODES(:) !<On return, USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)>=SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1)) THEN
            USER_ELEMENT_NODES=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES
          ELSE
            LOCAL_ERROR="The size of USER_ELEMENT_NODES is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(USER_ELEMENT_NODES,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1),"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesSet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODES(:) !<USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nn,NUMBER_OF_BAD_NODES,GLOBAL_NODE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:),BAD_NODES(:)
    LOGICAL :: ELEMENT_NODES_OK,NODE_EXISTS
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)==ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>MESH%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>MESH%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%PARENT_REGION
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(ASSOCIATED(NODES)) THEN
                  ELEMENT_NODES_OK=.TRUE.
                  ALLOCATE(GLOBAL_ELEMENT_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate global element nodes.",ERR,ERROR,*999)
                  ALLOCATE(BAD_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate bad nodes.",ERR,ERROR,*999)
                  NUMBER_OF_BAD_NODES=0
                  DO nn=1,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES
                    CALL NODE_CHECK_EXISTS(NODES,USER_ELEMENT_NODES(nn),NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                    IF(NODE_EXISTS) THEN
                      GLOBAL_ELEMENT_NODES(nn)=GLOBAL_NODE_NUMBER
                    ELSE
                      NUMBER_OF_BAD_NODES=NUMBER_OF_BAD_NODES+1
                      BAD_NODES(NUMBER_OF_BAD_NODES)=USER_ELEMENT_NODES(nn)
                      ELEMENT_NODES_OK=.FALSE.
                    ENDIF
                  ENDDO !nn
                  IF(ELEMENT_NODES_OK) THEN
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES=USER_ELEMENT_NODES
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES=GLOBAL_ELEMENT_NODES
                  ELSE
                    IF(NUMBER_OF_BAD_NODES==1) THEN
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in region "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
                          & " of parent region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))
                      DO nn=2,NUMBER_OF_BAD_NODES-1
                        LOCAL_ERROR=LOCAL_ERROR//","//TRIM(NUMBER_TO_VSTRING(BAD_NODES(nn),"*",ERR,ERROR))
                      ENDDO !nn
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" of parent region number "// &
                          &  TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ENDIF
                    ENDIF
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  IF(ASSOCIATED(REGION)) THEN                   
                    CALL FlagError("The elements mesh region does not have any associated nodes.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("The elements mesh interface does not have any associated nodes.",ERR,ERROR,*999) 
                  ENDIF
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of element nodes does not match number of basis nodes for this element.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets an element node's version for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesSet
  SUBROUTINE MeshElements_ElementNodeVersionSet(GLOBAL_NUMBER,ELEMENTS,VERSION_NUMBER,DERIVATIVE_NUMBER, &
      & USER_ELEMENT_NODE_INDEX,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The version number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODE_INDEX !< The node index of the specified element node to set a version for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementNodeVersionSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(USER_ELEMENT_NODE_INDEX>=1.AND.USER_ELEMENT_NODE_INDEX<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>mesh%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>mesh%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%PARENT_REGION
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(DERIVATIVE_NUMBER>=1.AND.DERIVATIVE_NUMBER<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS% &
                  & NUMBER_OF_DERIVATIVES(USER_ELEMENT_NODE_INDEX)) THEN !Check if the specified derivative exists
                  IF(VERSION_NUMBER>=1) THEN !Check if the specified version is greater than 1
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODE_VERSIONS(DERIVATIVE_NUMBER,USER_ELEMENT_NODE_INDEX) & 
                      & = VERSION_NUMBER
                    !\todo : There is redunancy in USER_ELEMENT_NODE_VERSIONS since it was allocated in MESH_TOPOLOGY_ELEMENTS_CREATE_START based on MAXIMUM_NUMBER_OF_DERIVATIVES for that elements basis:ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
                  ELSE
                    LOCAL_ERROR="The specified node version number of "//TRIM(NUMBER_TO_VSTRING(VERSION_NUMBER,"*", & 
                      & ERR,ERROR))//" is invalid. The element node index should be greater than 1."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified node derivative number of "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*", & 
                    & ERR,ERROR))//" is invalid. The element node derivative index should be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_DERIVATIVES( &
                    & USER_ELEMENT_NODE_INDEX),"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",err,error,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified element node index of "//TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NODE_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The element node index should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementNodeVersionSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementNodeVersionSet",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementNodeVersionSet

  !
  !================================================================================================================================
  !

 !>Calculates the element numbers surrounding an element in a mesh topology.
  SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the elements adjacent to elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,node_idx,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),NODE_POSITION_INDEX(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4)
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    ENTERS("MeshTopology_ElementsAdjacentElementsCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          !Loop over the global elements in the mesh
          DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
            !%%%% first we initialize lists that are required to find the adjacent elements list
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
            DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
              CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
            ENDDO !ni
            NUMBER_OF_NODES_XIC=1
            NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)=BASIS%NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)
            !Place the current element in the surrounding list
            CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER,ERR,ERROR,*999)
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              !Determine the collapsed "faces" if any
              NODE_POSITION_INDEX=1
              !Loop over the face normals of the element
              DO ni=1,BASIS%NUMBER_OF_XI
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Reset the node_position_index in this xi direction
                NODE_POSITION_INDEX(ni)=1
                !Loop over the two faces with this normal
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni
                  FACE_COLLAPSED(xi_direction)=.FALSE.
                  DO j=1,2
                    xi_dir_check=FACE_XI(j)
                    IF(xi_dir_check<=BASIS%NUMBER_OF_XI) THEN
                      xi_dir_search=FACE_XI(3-j)
                      NODE_POSITION_INDEX(xi_dir_search)=1
                      XI_COLLAPSED=.TRUE.
                      DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                        !Get the first local node along the xi check direction
                        NODE_POSITION_INDEX(xi_dir_check)=1
                        nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        !Get the second local node along the xi check direction
                        NODE_POSITION_INDEX(xi_dir_check)=2
                        nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        IF(nn1/=0.AND.nn2/=0) THEN
                          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn1)/= &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn2)) XI_COLLAPSED=.TRUE.
                        ENDIF
                        NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                      ENDDO !xi_dir_search
                      IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                    ENDIF
                  ENDDO !j
                  NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                ENDDO !direction_index
              ENDDO !ni
              !Loop over the xi directions and calculate the surrounding elements
              DO ni=1,BASIS%NUMBER_OF_XI
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Loop over the two faces
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni
                  !Find nodes in the element on the appropriate face/line/point
                  NULLIFY(NODE_MATCH_LIST)
                  CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)

                  CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                  IF(direction_index==-1) THEN
                    NODE_POSITION_INDEX(ni)=1
                  ELSE
                    NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                  ENDIF
                  !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is also
                  !collpased. This may indicate that we have a funny element in non-rc coordinates that goes around the central
                  !axis back to itself
                  IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                    !Do nothing - the match lists are already empty
                  ELSE
                    !Find the nodes to match and add them to the node match list
                    DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1))
                      NODE_POSITION_INDEX(FACE_XI(1))=nn1
                      DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2))
                        NODE_POSITION_INDEX(FACE_XI(2))=nn2
                        nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        IF(nn/=0) THEN
                          np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                          CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !nn2
                    ENDDO !nn1
                  ENDIF
                  CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                  NUMBER_SURROUNDING=0
                  IF(NUMBER_NODE_MATCHES>0) THEN
                    !Find list of elements surrounding those nodes
                    np1=NODE_MATCHES(1)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%numberOfSurroundingElements
                      ne1=TOPOLOGY%NODES%NODES(np1)%surroundingElements(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES ! should this be GLOBAL_ELEMENT_NODES?
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & MESH_ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDIF
                  IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                ENDDO !direction_index
              ENDDO !ni
            CASE(BASIS_SIMPLEX_TYPE)
              !Loop over the xi coordinates and calculate the surrounding elements
              DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                !Find the other coordinates of the face/line/point
                FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                !Find nodes in the element on the appropriate face/line/point
                NULLIFY(NODE_MATCH_LIST)
                CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                NODE_POSITION_INDEX(nic)=1 !Furtherest away from node with the nic'th coordinate
                !Find the nodes to match and add them to the node match list
                DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                  NODE_POSITION_INDEX(FACE_XIC(1))=nn1
                  DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                    NODE_POSITION_INDEX(FACE_XIC(2))=nn2
                    DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                      NODE_POSITION_INDEX(FACE_XIC(3))=nn3
                      nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3), &
                        NODE_POSITION_INDEX(4))
                      IF(nn/=0) THEN
                        np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                        CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !nn3
                  ENDDO !nn2
                ENDDO !nn1
                CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                IF(NUMBER_NODE_MATCHES>0) THEN
                  !Find list of elements surrounding those nodes
                  DO node_idx=1,NUMBER_NODE_MATCHES
                    np1=NODE_MATCHES(node_idx)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%numberOfSurroundingElements
                      ne1=TOPOLOGY%NODES%NODES(np1)%surroundingElements(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES 
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & MESH_ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDDO !node_idx
                ENDIF
                IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
              ENDDO !nic
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            !Set the surrounding elements for this element
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
              & BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements.",ERR,ERROR,*999)
            DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              CALL MESH_ADJACENT_ELEMENT_INITIALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
              CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS,ERR,ERROR,*999)
              ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element adjacent elements.",ERR,ERROR,*999)
              TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS) = ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS% &
                & ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)
              IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
            ENDDO !nic
          ENDDO !ne           
        ELSE
          CALL FlagError("Mesh topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not allocated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global element number : ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ERR,ERROR,*999)
          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
              & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,8,8,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)% &
              & ADJACENT_ELEMENTS,'("        Adjacent elements :",8(X,I8))','(30x,8(X,I8))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF
    
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
    ENDDO !ni
997 ERRORS("MeshTopology_ElementsAdjacentElementsCalculate",ERR,ERROR)
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements data structures for a mesh topology and deallocates any memory. \todo pass in elements
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
        CALL MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
      ENDDO !ne
      DEALLOCATE(ELEMENTS%ELEMENTS)
      IF(ASSOCIATED(ELEMENTS%ELEMENTS_TREE)) CALL TREE_DESTROY(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
      DEALLOCATE(ELEMENTS)
    ENDIF
 
    EXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Mesh already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%meshComponentTopology=>TOPOLOGY
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Mesh already has topology data points associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%totalNumberOfProjectedData=0
        TOPOLOGY%dataPoints%meshComponentTopology=>topology
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: ditto.
  
  !>Gets the user number for a global element identified by a given global number. \todo Check that the user number doesn't already exist. \see OPENCMISS::CMISSMeshElementsUserNumberGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the elements to get the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET

  !
  !================================================================================================================================
  !

  !>Returns the user number for a global element identified by a given global number. \see OPENCMISS::CMISSMeshElementsUserNumberGet
  SUBROUTINE MeshElements_ElementUserNumberGet(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberGet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN          
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Elements have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberGet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberGet",ERR,ERROR)
    RETURN 1

   
  END SUBROUTINE MeshElements_ElementUserNumberGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a global element identified by a given global number. \see OPENCMISS::CMISSMeshElementsUserNumberSet
  SUBROUTINE MeshElements_ElementUserNumberSet(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to set.
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the element to set
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,INSERT_STATUS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          NULLIFY(TREE_NODE)
          CALL TREE_SEARCH(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
          IF(ASSOCIATED(TREE_NODE)) THEN
            CALL TREE_NODE_VALUE_GET(ELEMENTS%ELEMENTS_TREE,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
            LOCAL_ERROR="Element user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " is already used by global element number "// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            CALL TREE_ITEM_DELETE(ELEMENTS%ELEMENTS_TREE,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
            CALL TREE_ITEM_INSERT(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
            ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberSet",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MeshElements_ElementUserNumberSet
  
  !
  !================================================================================================================================
  !

  !>Changes/sets the user numbers for all elements.
  SUBROUTINE MeshTopologyElementsUserNumbersAllSet(elements,userNumbers,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: elements !<A pointer to the elements to set all the user numbers for 
    INTEGER(INTG), INTENT(IN) :: userNumbers(:) !<The user numbers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,insertStatus
    TYPE(TREE_TYPE), POINTER :: newElementsTree
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newElementsTree)

    ENTERS("MeshTopologyElementsUserNumbersAllSet",err,error,*999)

    IF(ASSOCIATED(elements)) THEN
      IF(elements%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",err,error,*999)
      ELSE
        IF(elements%NUMBER_OF_ELEMENTS==SIZE(userNumbers,1)) THEN
          !Check the users numbers to ensure that there are no duplicates          
          CALL TREE_CREATE_START(newElementsTree,err,error,*999)
          CALL TREE_INSERT_TYPE_SET(newElementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
          CALL TREE_CREATE_FINISH(newElementsTree,err,error,*999)
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            CALL TREE_ITEM_INSERT(newElementsTree,userNumbers(elementIdx),elementIdx,insertStatus,err,error,*999)
            IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) THEN
              localError="The specified user number of "//TRIM(NumberToVstring(userNumbers(elementIdx),"*",err,error))// &
                & " for global element number "//TRIM(NUMBER_TO_VSTRING(elementIdx,"*",err,error))// &
                & " is a duplicate. The user element numbers must be unique."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !elementIdx
          CALL TREE_DESTROY(elements%ELEMENTS_TREE,err,error,*999)
          elements%ELEMENTS_TREE=>newElementsTree
          NULLIFY(newElementsTree)
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            elements%ELEMENTS(elementIdx)%GLOBAL_NUMBER=elementIdx
            elements%ELEMENTS(elementIdx)%USER_NUMBER=userNumbers(elementIdx)
          ENDDO !elementIdx
        ELSE
          localError="The number of specified element user numbers ("// &
            TRIM(NumberToVstring(SIZE(userNumbers,1),"*",err,error))// &
            ") does not match number of elements ("// &
            TRIM(NumberToVstring(elements%NUMBER_OF_ELEMENTS,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyElementsUserNumbersAllSet")
    RETURN
999 IF(ASSOCIATED(newElementsTree)) CALL TREE_DESTROY(newElementsTree,err,error,*998)
998 ERRORSEXITS("MeshTopologyElementsUserNumbersAllSet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopologyElementsUserNumbersAllSet
  
  !
  !================================================================================================================================
  !

  !>Calculates the data points in the given mesh topology.
  SUBROUTINE MeshTopologyDataPointsCalculateProjection(mesh,dataProjection,err,error,*)
  
    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh topology to calcualte the data projection for
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection !<A pointer to the data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DATA_POINTS_TYPE), POINTER :: dataPoints !<A pointer to the data points
    TYPE(MeshDataPointsType), POINTER :: dataPointsTopology
    TYPE(DATA_PROJECTION_RESULT_TYPE), POINTER :: dataProjectionResult
    TYPE(MeshElementsType), POINTER :: elements
    INTEGER(INTG) :: dataPointIdx,elementIdx,countIdx,projectionNumber,globalCountIdx,elementNumber

    ENTERS("MeshTopologyDataPointsCalculateProjection",ERR,ERROR,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(dataProjection%DATA_PROJECTION_FINISHED) THEN 
        dataPoints=>dataProjection%DATA_POINTS
        !Default the first mesh component topology to contain data points ! \TODO: need to be changed once the data points topology is moved under meshTopologyType.
        dataPointsTopology=>mesh%TOPOLOGY(1)%PTR%dataPoints
        !Extract the global number of the data projection 
        projectionNumber=dataProjection%GLOBAL_NUMBER
        !Hard code the first mesh component since element topology is the same for all mesh components
        !\TODO: need to be changed once the elements topology is moved under meshTopologyType.
        elements=>mesh%TOPOLOGY(1)%PTR%ELEMENTS
        ALLOCATE(dataPointsTopology%elementDataPoint(elements%NUMBER_OF_ELEMENTS),STAT=ERR)     
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology element.",ERR,ERROR,*999)
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          dataPointsTopology%elementDataPoint(elementIdx)%elementNumber=elements%ELEMENTS(elementIdx)%GLOBAL_NUMBER
          dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData=0
        ENDDO        
        !Calculate number of projected data points on an element
        DO dataPointIdx=1,dataPoints%NUMBER_OF_DATA_POINTS
          dataProjectionResult=>dataProjection%DATA_PROJECTION_RESULTS(dataPointIdx)
          elementNumber=dataProjectionResult%ELEMENT_NUMBER
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
              dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData= &
                & dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData+1;
            ENDIF
          ENDDO !elementIdx
        ENDDO       
        !Allocate memory to store data indices and initialise them to be zero   
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          ALLOCATE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(dataPointsTopology% &
            & elementDataPoint(elementIdx)%numberOfProjectedData),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate data points topology element data points.",ERR,ERROR,*999)
          DO countIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=0
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=0
          ENDDO
        ENDDO     
        !Record the indices of the data that projected on the elements 
        globalCountIdx=0
        dataPointsTopology%totalNumberOfProjectedData=0
        DO dataPointIdx=1,dataPoints%NUMBER_OF_DATA_POINTS 
          dataProjectionResult=>dataProjection%DATA_PROJECTION_RESULTS(dataPointIdx)
          elementNumber=dataProjectionResult%ELEMENT_NUMBER
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            countIdx=1         
            IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
              globalCountIdx=globalCountIdx+1
              !Find the next data point index in this element
              DO WHILE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber/=0)
                countIdx=countIdx+1
              ENDDO
              dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=dataPointIdx
              dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=dataPointIdx!globalCountIdx (used this if only projected data are taken into account)
              dataPointsTopology%totalNumberOfProjectedData=dataPointsTopology%totalNumberOfProjectedData+1
            END IF            
          ENDDO !elementIdx
        ENDDO !dataPointIdx
        !Allocate memory to store total data indices in ascending order and element map
        ALLOCATE(dataPointsTopology%dataPoints(dataPointsTopology%totalNumberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology data points.",ERR,ERROR,*999)
        !The global number for the data points will be looping through elements.
        countIdx=1  
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          DO dataPointIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%dataPoints(countIdx)%userNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%userNumber
             dataPointsTopology%dataPoints(countIdx)%globalNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%globalNumber
             dataPointsTopology%dataPoints(countIdx)%elementNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
               & elementNumber
             countIdx=countIdx+1
          ENDDO !dataPointIdx
        ENDDO !elementIdx                      
      ELSE
        CALL FlagError("Data projection is not finished.",err,error,*999)
      END IF    
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyDataPointsCalculateProjection")
    RETURN
999 ERRORSEXITS("MeshTopologyDataPointsCalculateProjection",err,error)
    RETURN 1
  END SUBROUTINE MeshTopologyDataPointsCalculateProjection

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. \todo pass in the mesh topology
  SUBROUTINE MeshTopologyFinalise(mesh,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx

    ENTERS("MeshTopologyFinalise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
        CALL MeshTopologyComponentFinalise(mesh%topology(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      DEALLOCATE(mesh%topology)
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
 
    EXITS("MeshTopologyFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopologyFinalise

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. 
  SUBROUTINE MeshTopologyComponentFinalise(meshComponent,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: meshComponent !<A pointer to the mesh component to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyComponentFinalise",err,error,*999)

    IF(ASSOCIATED(meshComponent)) THEN
      CALL MeshTopologyNodesFinalise(meshComponent%nodes,err,error,*999)
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(meshComponent%elements,err,error,*999)
      CALL MeshTopologyDofsFinalise(meshComponent%dofs,err,error,*999)
      DEALLOCATE(meshComponent)
    ENDIF
 
    EXITS("MeshTopologyComponentFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyComponentFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopologyComponentFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given mesh. \todo finalise on error
  SUBROUTINE MeshTopologyInitialise(mesh,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to initialise the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    
    ENTERS("MeshTopologyInitialise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(mesh%topology)) THEN
        CALL FlagError("Mesh already has topology associated.",err,error,*999)
      ELSE
        !Allocate mesh topology
        ALLOCATE(mesh%topology(mesh%NUMBER_OF_COMPONENTS),STAT=err)
        IF(err/=0) CALL FlagError("Mesh topology could not be allocated.",err,error,*999)
        DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
          ALLOCATE(mesh%topology(componentIdx)%ptr,STAT=err)
          IF(err/=0) CALL FlagError("Mesh topology component could not be allocated.",err,error,*999)
          mesh%topology(componentIdx)%ptr%mesh=>mesh
          NULLIFY(mesh%topology(componentIdx)%ptr%elements)
          NULLIFY(mesh%topology(componentIdx)%ptr%nodes)
          NULLIFY(mesh%topology(componentIdx)%ptr%dofs)
          NULLIFY(mesh%topology(componentIdx)%ptr%dataPoints)
          !Initialise the topology components
          CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MeshTopologyNodesInitialise(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MeshTopologyDofsInitialise(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
        ENDDO !componentIdx
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyInitialise")
    RETURN
999 ERRORSEXITS("MeshTopologyInitialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopologyInitialise
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh component. 
  SUBROUTINE MeshTopologyElementCheckExistsMesh(mesh,meshComponentNumber,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to check the element exists on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component to check the element exits on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshElementsType), POINTER :: elements

    NULLIFY(elements)
    
    ENTERS("MeshTopologyElementCheckExistsMesh",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%MESH_FINISHED) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_GET(mesh,meshComponentNumber,elements,err,error,*999)
        CALL MeshTopologyElementCheckExistsMeshElements(elements,userElementNumber,elementExists,globalElementNumber,err,error,*999)
      ELSE
        CALL FlagError("Mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyElementCheckExistsMesh")
    RETURN
999 ERRORSEXITS("MeshTopologyElementCheckExistsMesh",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyElementCheckExistsMesh
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh elements. 
  SUBROUTINE MeshTopologyElementCheckExistsMeshElements(meshElements,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh elements, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("MeshTopologyElementCheckExistsMesh",err,error,*999)

    elementExists=.FALSE.
    globalElementNumber=0
    IF(ASSOCIATED(meshElements)) THEN
      NULLIFY(treeNode)
      CALL TREE_SEARCH(meshElements%ELEMENTS_TREE,userElementNumber,treeNode,err,error,*999)
      IF(ASSOCIATED(treeNode)) THEN
        CALL TREE_NODE_VALUE_GET(meshElements%ELEMENTS_TREE,treeNode,globalElementNumber,err,error,*999)
        elementExists=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyElementCheckExistsMeshElements")
    RETURN
999 ERRORSEXITS("MeshTopologyElementCheckExistsMeshElements",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyElementCheckExistsMeshElements
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh component. 
  SUBROUTINE MeshTopologyNodeCheckExistsMesh(mesh,meshComponentNumber,userNodeNumber,nodeExists,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to check the node exists on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component to check the node exits on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNodeNumber
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode

    NULLIFY(meshNodes)
    NULLIFY(nodes)
    NULLIFY(region)
    
    ENTERS("MeshTopologyNodeCheckExistsMesh",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%MESH_FINISHED) THEN
        CALL MeshTopologyNodesGet(mesh,meshComponentNumber,meshNodes,err,error,*999)
        CALL MeshRegionGet(mesh,region,err,error,*999)
        nodes=>region%nodes
        IF(ASSOCIATED(nodes)) THEN
          CALL NODE_CHECK_EXISTS(nodes,userNodeNumber,nodeExists,globalNodeNumber,err,error,*999)
          NULLIFY(treeNode)
          CALL TREE_SEARCH(meshNodes%nodesTree,globalNodeNumber,treeNode,err,error,*999)
          IF(ASSOCIATED(treeNode)) THEN
            CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
            nodeExists=.TRUE.
          ENDIF
        ELSE
          CALL FlagError("Region nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyNodeCheckExistsMesh")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeCheckExistsMesh",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodeCheckExistsMesh
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh nodes. 
  SUBROUTINE MeshTopologyNodeCheckExistsMeshNodes(meshNodes,userNodeNumber,nodeExists,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to check the node exists on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh nodes, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNodeNumber
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode

    NULLIFY(nodes)
    NULLIFY(region)
    
    ENTERS("MeshTopologyNodeCheckExistsMeshNodes",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(ASSOCIATED(meshNodes)) THEN
      meshComponentTopology=>meshNodes%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          IF(mesh%MESH_FINISHED) THEN
            CALL MeshRegionGet(mesh,region,err,error,*999)
            nodes=>region%nodes
            IF(ASSOCIATED(nodes)) THEN
              CALL NODE_CHECK_EXISTS(nodes,userNodeNumber,nodeExists,globalNodeNumber,err,error,*999)
              NULLIFY(treeNode)
              CALL TREE_SEARCH(meshNodes%nodesTree,globalNodeNumber,treeNode,err,error,*999)
              IF(ASSOCIATED(treeNode)) THEN
                CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
                nodeExists=.TRUE.
              ENDIF
            ELSE
              CALL FlagError("Region nodes is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh has not been finished.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyNodeCheckExistsMeshNodes")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeCheckExistsMeshNodes",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodeCheckExistsMeshNodes
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MeshTopologyNodeFinalise(node,err,error,*)

    !Argument variables
    TYPE(MeshNodeType) :: node !<The mesh node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx

    ENTERS("MeshTopologyNodeFinalise",err,error,*999)

    IF(ALLOCATED(node%derivatives)) THEN
      DO derivativeIdx=1,node%numberOfDerivatives
        CALL MeshTopologyNodeDerivativeFinalise(node%derivatives(derivativeIdx),err,error,*999)
      ENDDO !derivativeIdx
      DEALLOCATE(node%derivatives)
    ENDIF
    IF(ASSOCIATED(node%surroundingElements)) DEALLOCATE(node%surroundingElements)
  
    EXITS("MeshTopologyNodeFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MeshTopologyNodeInitialise(node,err,error,*)

    !Argument variables
    TYPE(MeshNodeType) :: node !<The mesh node to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyNodeInitialise",err,error,*999)

    node%userNumber=0
    node%globalNumber=0
    node%numberOfSurroundingElements=0
    NULLIFY(node%surroundingElements)
    node%numberOfDerivatives=0
    node%boundaryNode=.FALSE.
    
    EXITS("MeshTopologyNodeInitialise")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeInitialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopologyNodeInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the nodes used the mesh identified by a given mesh topology.
  SUBROUTINE MeshTopologyNodesCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,elementIdx,insertStatus,localNodeIdx,globalNode,meshNodeIdx,meshNode,numberOfNodes 
    INTEGER(INTG), POINTER :: globalNodeNumbers(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(TREE_TYPE), POINTER :: globalNodesTree
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(globalNodeNumbers)
    NULLIFY(globalNodesTree)
    
    ENTERS("MeshTopologyNodesCalculate",err,error,*998)

    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements
      IF(ASSOCIATED(elements)) THEN
        meshNodes=>topology%nodes
        IF(ASSOCIATED(meshNodes)) THEN
          mesh=>topology%mesh
          IF(ASSOCIATED(mesh)) THEN
            NULLIFY(nodes)
            CALL MeshGlobalNodesGet(mesh,nodes,err,error,*999)
            IF(ALLOCATED(meshNodes%nodes)) THEN
              localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
                & " already has allocated mesh topology nodes."
              CALL FlagError(localError,err,error,*998)
            ELSE
              !Work out what nodes are in the mesh
              CALL TREE_CREATE_START(globalNodesTree,err,error,*999)
              CALL TREE_INSERT_TYPE_SET(globalNodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
              CALL TREE_CREATE_FINISH(globalNodesTree,err,error,*999)
              DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
                basis=>elements%elements(elementIdx)%basis
                DO localNodeIdx=1,basis%NUMBER_OF_NODES
                  globalNode=elements%elements(elementIdx)%GLOBAL_ELEMENT_NODES(localNodeIdx)
                  CALL TREE_ITEM_INSERT(globalNodesTree,globalNode,globalNode,insertStatus,err,error,*999)
                ENDDO !localNodeIdx
              ENDDO !elementIdx
              CALL TREE_DETACH_AND_DESTROY(globalNodesTree,numberOfNodes,globalNodeNumbers,err,error,*999)
              !Set up the mesh nodes.
              ALLOCATE(meshNodes%nodes(numberOfNodes),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate mesh topology nodes nodes.",err,error,*999)
              CALL TREE_CREATE_START(meshNodes%nodesTree,err,error,*999)
              CALL TREE_INSERT_TYPE_SET(meshNodes%nodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
              CALL TREE_CREATE_FINISH(meshNodes%nodesTree,err,error,*999) 
              DO meshNodeIdx=1,numberOfNodes
                CALL MeshTopologyNodeInitialise(meshNodes%nodes(meshNodeIdx),err,error,*999)
                meshNodes%nodes(meshNodeIdx)%meshNumber=meshNodeIdx
                meshNodes%nodes(meshNodeIdx)%globalNumber=globalNodeNumbers(meshNodeIdx)
                meshNodes%nodes(meshNodeIdx)%userNumber=nodes%nodes(globalNodeNumbers(meshNodeIdx))%USER_NUMBER
                CALL TREE_ITEM_INSERT(meshNodes%nodesTree,globalNodeNumbers(meshNodeIdx),meshNodeIdx,insertStatus,err,error,*999)
              ENDDO !nodeIdx
              meshNodes%numberOfNodes=numberOfNodes
              IF(ASSOCIATED(globalNodeNumbers)) DEALLOCATE(globalNodeNumbers)
              !Now recalculate the mesh element nodes
              DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
                basis=>elements%elements(elementIdx)%basis
                ALLOCATE(elements%elements(elementIdx)%MESH_ELEMENT_NODES(basis%NUMBER_OF_NODES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate mesh topology elements mesh element nodes.",err,error,*999)
                DO localNodeIdx=1,basis%NUMBER_OF_NODES
                  globalNode=elements%elements(elementIdx)%GLOBAL_ELEMENT_NODES(localNodeIdx)
                  NULLIFY(treeNode)
                  CALL TREE_SEARCH(meshNodes%nodesTree,globalNode,treeNode,err,error,*999)
                  IF(ASSOCIATED(treeNode)) THEN
                    CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNode,err,error,*999)
                    elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)=meshNode
                  ELSE
                    localError="Could not find global node "//TRIM(NumberToVString(globalNode,"*",err,error))//" (user node "// &
                      & TRIM(NumberToVString(nodes%nodes(globalNode)%USER_NUMBER,"*",err,error))//") in the mesh nodes."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDDO !localNodeIdx
              ENDDO !elementIdx
            ENDIF
          ELSE
            CALL FlagError("Mesh topology mesh is not associated.",err,error,*998)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh nodes = ",meshNodes%numberOfNodes,err,error,*999)
      DO meshNodeIdx=1,meshNodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh node number = ",meshNodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global node number = ",meshNodes%nodes(meshNodeIdx)%globalNumber, &
          & err,error,*999)        
      ENDDO !meshNodeIdx
    ENDIF
    
    EXITS("MeshTopologyNodesCalculate")
    RETURN
999 IF(ASSOCIATED(globalNodeNumbers)) DEALLOCATE(globalNodeNumbers)
    IF(ASSOCIATED(globalNodesTree)) CALL TREE_DESTROY(globalNodesTree,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopologyNodesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopologyNodesCalculate

  !
  !================================================================================================================================
  !

  !>Destroys the nodes in a mesh topology.
  SUBROUTINE MeshTopologyNodesDestroy(nodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: nodes !<A pointer to the mesh nodes to destroy 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyNodesDestroy",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      CALL MeshTopologyNodesFinalise(nodes,err,error,*999)
    ELSE
      CALL FlagError("Mesh topology nodes is not associated",err,error,*999)
    ENDIF

    EXITS("MeshTopologyNodesDestroy")
    RETURN
999 ERRORSEXITS("MeshTopologyNodesDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodesDestroy
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh nodes for a given mesh component.
  SUBROUTINE MeshTopologyNodesGet(mesh,meshComponentNumber,nodes,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the nodes for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to get the nodes for
    TYPE(MeshNodesType), POINTER :: nodes !<On return, a pointer to the mesh nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("MeshTopologyNodesGet",err,error,*998)
    
    IF(ASSOCIATED(mesh)) THEN
      IF(meshComponentNumber>0.AND.meshComponentNumber<=mesh%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(nodes)) THEN
          CALL FlagError("Nodes is already associated.",err,error,*998)
        ELSE
          IF(ASSOCIATED(mesh%topology(meshComponentNumber)%ptr)) THEN
            IF(ASSOCIATED(mesh%topology(meshComponentNumber)%ptr%nodes)) THEN
              nodes=>mesh%topology(meshComponentNumber)%ptr%nodes
            ELSE
              CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh topology is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NumberToVString(mesh%NUMBER_OF_COMPONENTS,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",err,error,*998)
    ENDIF

    EXITS("MeshTopologyNodesGet")
    RETURN
999 NULLIFY(nodes)
998 ERRORSEXITS("MeshTopologyNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodesGet

  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MeshTopologyNodeDerivativeFinalise(nodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: nodeDerivative !<The mesh node derivative to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyNodeDerivativeFinalise",err,error,*999)

    IF(ALLOCATED(nodeDerivative%userVersionNumbers)) DEALLOCATE(nodeDerivative%userVersionNumbers)
    IF(ALLOCATED(nodeDerivative%dofIndex)) DEALLOCATE(nodeDerivative%dofIndex)

    EXITS("MeshTopologyNodeDerivativeFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeDerivativeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodeDerivativeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MeshTopologyNodeDerivativeInitialise(nodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: nodeDerivative !<The mesh node derivative to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyNodeDerivativeInitialise",err,error,*999)

    nodeDerivative%numberOfVersions=0
    nodeDerivative%globalDerivativeIndex=0
    nodederivative%partialDerivativeIndex=0

    EXITS("MeshTopologyNodeDerivativeInitialise")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeDerivativeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodeDerivativeInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the number of derivatives at each node in a topology.
  SUBROUTINE MeshTopologyNodesDerivativesCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the derivates at each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,element,elementIdx,globalDerivative,localNodeIdx,maxNumberOfDerivatives,nodeIdx, &
      & numberOfDerivatives
    INTEGER(INTG), ALLOCATABLE :: derivatives(:)
    LOGICAL :: found
    TYPE(LIST_TYPE), POINTER :: nodeDerivativeList
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopologyNodesDerivativesCalculate",err,ERROR,*999)

     IF(ASSOCIATED(topology)) THEN
       elements=>topology%elements
       IF(ASSOCIATED(elements)) THEN
         nodes=>topology%nodes
         IF(ASSOCIATED(nodes)) THEN
          !Loop over the mesh nodes
          DO nodeIdx=1,nodes%numberOfNodes
            !Calculate the number of derivatives and versions at each node. This needs to be calculated by looking at the
            !mesh elements as we may have an adjacent element in another domain with a higher order basis also with versions.
            NULLIFY(nodeDerivativeList)
            CALL LIST_CREATE_START(nodeDerivativeList,err,error,*999)
            CALL LIST_DATA_TYPE_SET(nodeDerivativeList,LIST_INTG_TYPE,err,error,*999)
            CALL LIST_INITIAL_SIZE_SET(nodeDerivativeList,8,err,error,*999)
            CALL LIST_CREATE_FINISH(nodeDerivativeList,err,error,*999)
            maxNumberOfDerivatives=-1
            DO elementIdx=1,nodes%nodes(nodeIdx)%numberOfSurroundingElements
              element=nodes%nodes(nodeIdx)%surroundingElements(elementIdx)
              basis=>elements%elements(element)%basis
              !Find the local node corresponding to this node
              found=.FALSE.
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                IF(elements%elements(element)%MESH_ELEMENT_NODES(localNodeIdx)==nodeIdx) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !nn
              IF(found) THEN
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                  CALL LIST_ITEM_ADD(nodeDerivativeList,basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx),err,error,*999)
                ENDDO !derivativeIdx
                IF(basis%NUMBER_OF_DERIVATIVES(localNodeIdx)>maxNumberOfDerivatives) &
                  & maxNumberOfDerivatives=basis%NUMBER_OF_DERIVATIVES(localNodeidx)
              ELSE
                CALL FlagError("Could not find local node.",err,error,*999)
              ENDIF
            ENDDO !elem_idx
            CALL LIST_REMOVE_DUPLICATES(nodeDerivativeList,err,error,*999)
            CALL LIST_DETACH_AND_DESTROY(nodeDerivativeList,numberOfDerivatives,derivatives,err,error,*999)
            IF(numberOfDerivatives==maxNumberOfDerivatives) THEN
              !Set up the node derivatives.
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(maxNumberOfDerivatives),STAT=err)
              nodes%nodes(nodeIdx)%numberOfDerivatives=maxNumberOfDerivatives
              DO derivativeIdx=1,numberOfDerivatives
                CALL MeshTopologyNodeDerivativeInitialise(nodes%nodes(nodeIdx)%derivatives(derivativeIdx),err,error,*999)
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex = derivatives(derivativeIdx)
                globalDerivative=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(derivatives(derivativeIdx))
                IF(globalDerivative/=0) THEN
                  nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex=globalDerivative
                ELSE
                  localError="The partial derivative index of "//TRIM(NumberToVstring(derivatives(derivativeIdx),"*", &
                    & err,error))//" for derivative number "//TRIM(NumberToVstring(derivativeIdx,"*",err,error))// &
                    & " does not have a corresponding global derivative."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
              DEALLOCATE(derivatives)
            ELSE
              localError="Invalid mesh configuration. User node "// &
                & TRIM(NumberToVstring(nodes%nodes(nodeIdx)%userNumber,"*",err,error))// &
                & " has inconsistent derivative directions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !nodeIdx
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*999)
    ENDIF
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ",nodes%nodes(nodeIdx)%numberOfDerivatives, &
          & err,error,*999)
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          !TODO: change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
        ENDDO !derivativeIdx
      ENDDO !node_idx
    ENDIF
    
    EXITS("MeshTopologyNodesDerivativesCalculate")
    RETURN
999 IF(ALLOCATED(derivatives)) DEALLOCATE(derivatives)
    IF(ASSOCIATED(nodeDerivativeList)) CALL LIST_DESTROY(nodeDerivativeList,err,error,*998)
998 ERRORSEXITS("MeshTopologyNodesDerivativesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopologyNodesDerivativesCalculate

  !
  !================================================================================================================================
  !
  
  !>Returns the number of derivatives for a node in a mesh
  SUBROUTINE MeshTopologyNodeNumberOfDerivativesGet(meshNodes,userNumber,numberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of derivatives for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user node number to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfDerivatives !<On return, the number of global derivatives at the specified user node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber,meshNumber
    LOGICAL :: nodeExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopologyNodeNumberOfDerivativesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopologyNodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN
        numberOfDerivatives=meshNodes%nodes(meshNumber)%numberOfDerivatives
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyNodeNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeNumberOfDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopologyNodeNumberOfDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the global derivative numbers for a node in mesh nodes
  SUBROUTINE MeshTopologyNodeDerivativesGet(meshNodes,userNumber,derivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the derivatives for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the node to get the derivatives for
    INTEGER(INTG), INTENT(OUT) :: derivatives(:) !<On return, the global derivatives at the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,meshComponentNumber,meshNumber,numberOfDerivatives
    LOGICAL :: nodeExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopologyNodeDerivativesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopologyNodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN         
        numberOfDerivatives=meshNodes%nodes(meshNumber)%numberOfDerivatives
        IF(SIZE(derivatives,1)>=numberOfDerivatives) THEN
          DO derivativeIdx=1,numberOfDerivatives
            derivatives(derivativeIdx)=meshNodes%nodes(meshNumber)%derivatives(derivativeIdx)%globalDerivativeIndex
          ENDDO !derivativeIdx
        ELSE
          localError="The size of the supplied derivatives array of "// &
            & TRIM(NumberToVString(SIZE(derivatives,1),"*",err,error))// &
            & " is too small. The size should be >= "// &
            & TRIM(NumberToVString(numberOfDerivatives,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyNodeDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopologyNodeDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of versions for a derivative of a node in mesh nodes
  SUBROUTINE MeshTopologyNodeNumberOfVersionsGet(meshNodes,derivativeNumber,userNumber,numberOfVersions,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the node to get the number of versions for
    INTEGER(INTG), INTENT(OUT) :: numberOfVersions !<On return, the number of versions for the specified derivative of the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber,meshNumber
    LOGICAL :: nodeExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopologyNodeNumberOfVersionsGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopologyNodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN
        IF(derivativeNumber>=1.AND.derivativeNumber<=meshNodes%nodes(meshNumber)%numberOfDerivatives) THEN
          numberOfVersions=meshNodes%nodes(meshNumber)%derivatives(derivativeNumber)%numberOfVersions
        ELSE
          localError="The specified derivative index of "// &
            & TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
            & " is invalid. The derivative index must be >= 1 and <= "// &
            & TRIM(NumberToVString(meshNodes%nodes(meshNumber)%numberOfDerivatives,"*",err,error))// &
            & " for user node number "//TRIM(NumberToVString(userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyNodeNumberOfVersionsGet")
    RETURN
999 ERRORSEXITS("MeshTopologyNodeNumberOfVersionsGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopologyNodeNumberOfVersionsGet
  
  !
  !================================================================================================================================
  !

  !>Returns the number of nodes for a node in a mesh
  SUBROUTINE MeshTopologyNodesNumberOfNodesGet(meshNodes,numberOfNodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On return, the number of nodes in the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("MeshTopologyNodesNumberOfNodesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      numberOfNodes=meshNodes%numberOfNodes
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyNodesNumberOfNodesGet")
    RETURN
999 ERRORSEXITS("MeshTopologyNodesNumberOfNodesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopologyNodesNumberOfNodesGet
  
  !
  !================================================================================================================================
  !

  !>Calculates the number of versions at each node in a topology.
  SUBROUTINE MeshTopologyNodesVersionCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the versions at each node for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element,localNodeIdx,derivativeIdx,nodeIdx,numberOfVersions,versionIdx
    INTEGER(INTG), ALLOCATABLE :: versions(:)
    TYPE(LIST_PTR_TYPE), POINTER :: nodeVersionList(:,:)
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(BASIS_TYPE), POINTER :: basis
    
    ENTERS("MeshTopologyNodesVersionCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes
        IF(ASSOCIATED(nodes)) THEN
          !Loop over the mesh elements
          !Calculate the number of versions at each node. This needs to be calculated by looking at all the mesh elements
          !as we may have an adjacent elements in another domain with a higher order basis along with different versions
          !being assigned to its derivatives.
          !\todo : See if there are any constraints that can be applied to restrict the amount of memory being allocated here
          ALLOCATE(nodeVersionList(MAXIMUM_GLOBAL_DERIV_NUMBER,nodes%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node version list.",err,error,*999)
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              NULLIFY(nodeVersionList(derivativeIdx,nodeIdx)%ptr)
              CALL LIST_CREATE_START(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
              CALL LIST_DATA_TYPE_SET(nodeVersionList(derivativeIdx,nodeIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
              CALL LIST_INITIAL_SIZE_SET(nodeVersionList(derivativeIdx,nodeIdx)%ptr,8,err,error,*999)
              CALL LIST_CREATE_FINISH(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
            ENDDO!derivativeIdx
          ENDDO!nodeIdx
          DO element=1,elements%NUMBER_OF_ELEMENTS
            basis=>elements%elements(element)%basis
            DO localNodeIdx=1,BASIS%NUMBER_OF_NODES
              DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                CALL LIST_ITEM_ADD(nodeVersionList(derivativeIdx,elements%elements(element)% &
                  & MESH_ELEMENT_NODES(localNodeIdx))%ptr,elements%elements(element)%USER_ELEMENT_NODE_VERSIONS( &
                  & derivativeIdx,localNodeIdx),err,error,*999)
              ENDDO!derivativeIdx
            ENDDO!localNodeIdx
          ENDDO!element
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              CALL LIST_REMOVE_DUPLICATES(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
              CALL LIST_DETACH_AND_DESTROY(nodeVersionList(derivativeIdx,nodeIdx)%ptr,numberOfVersions,versions, &
                & err,error,*999)
              nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions = MAXVAL(versions(1:numberOfVersions))
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(nodes%nodes(nodeIdx)% &
                & derivatives(derivativeIdx)%numberOfVersions),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate node global derivative index.",err,error,*999)
              DO versionIdx=1,nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions 
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx) = versionIdx
              ENDDO !versionIdx
              DEALLOCATE(versions)
            ENDDO!derivativeIdx
          ENDDO!nodeIdx
          DEALLOCATE(nodeVersionList)
          NULLIFY(nodeVersionList)
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*999)
    ENDIF
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ", &
          & nodes%nodes(nodeIdx)%numberOfDerivatives,err,error,*999)
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          !\todo : change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions,8,8, &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers, & 
            & '("    User Version index(derivativeIdx,:) :",8(X,I2))','(36X,8(X,I2))',err,error,*999)
        ENDDO!derivativeIdx
      ENDDO !nodeIdx
    ENDIF
    
    EXITS("MeshTopologyNodesVersionCalculate")
    RETURN
999 IF(ALLOCATED(versions)) DEALLOCATE(versions)
    IF(ASSOCIATED(nodeVersionList)) THEN
      DO nodeIdx=1,nodes%numberOfNodes
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          CALL LIST_DESTROY(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*998)
        ENDDO !derivativeIdx
      ENDDO !nodeIdx
      DEALLOCATE(nodeVersionList)
    ENDIF
998 ERRORSEXITS("MeshTopologyNodesVersionCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopologyNodesVersionCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a mesh.
  SUBROUTINE MeshTopologySurroundingElementsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the elements surrounding each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: element,elementIdx,insertPosition,localNodeIdx,node,surroundingElementNumber
    INTEGER(INTG), POINTER :: newSurroundingElements(:)
    LOGICAL :: foundElement
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes

    NULLIFY(newSurroundingElements)

    ENTERS("MeshTopologySurroundingElementsCalculate",err,error,*999)
    
    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements      
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes       
        IF(ASSOCIATED(nodes)) THEN
          IF(ALLOCATED(nodes%nodes)) THEN
            DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
              basis=>elements%elements(elementIdx)%basis
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                node=elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)
                foundElement=.FALSE.
                element=1
                insertPosition=1
                DO WHILE(element<=nodes%nodes(node)%numberOfSurroundingElements.AND..NOT.foundElement)
                  surroundingElementNumber=nodes%nodes(node)%surroundingElements(element)
                  IF(surroundingElementNumber==elementIdx) THEN
                    foundElement=.TRUE.
                  ENDIF
                  element=element+1
                  IF(elementIdx>=surroundingElementNumber) THEN
                    insertPosition=element
                  ENDIF
                ENDDO
                IF(.NOT.foundElement) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(newSurroundingElements(nodes%nodes(node)%numberOfSurroundingElements+1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate new surrounding elements.",err,error,*999)
                  IF(ASSOCIATED(nodes%nodes(node)%surroundingElements)) THEN
                    newSurroundingElements(1:insertPosition-1)=nodes%nodes(node)%surroundingElements(1:insertPosition-1)
                    newSurroundingElements(insertPosition)=elementIdx
                    newSurroundingElements(insertPosition+1:nodes%nodes(node)%numberOfSurroundingElements+1)= &
                      & nodes%nodes(node)%surroundingElements(insertPosition:nodes%nodes(node)%numberOfSurroundingElements)
                    DEALLOCATE(nodes%nodes(node)%surroundingElements)
                  ELSE
                    newSurroundingElements(1)=elementIdx
                  ENDIF
                  nodes%nodes(node)%surroundingElements=>newSurroundingElements
                  nodes%nodes(node)%numberOfSurroundingElements=nodes%nodes(node)%numberOfSurroundingElements+1
                ENDIF
              ENDDO !localNodeIdx
            ENDDO !elementIdx
          ELSE
            CALL FlagError("Mesh topology nodes nodes have not been allocated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology nodes are not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology not associated.",err,error,*999)
    ENDIF

    EXITS("MeshTopologySurroundingElementsCalculate")
    RETURN
999 IF(ASSOCIATED(newSurroundingElements)) DEALLOCATE(newSurroundingElements)
    ERRORSEXITS("MeshTopologySurroundingElementsCalculate",err,error)
    RETURN 1   
  END SUBROUTINE MeshTopologySurroundingElementsCalculate
  
  !
  !===============================================================================================================================
  !

  !>Finalises the nodes data structures for a mesh topology and deallocates any memory.
  SUBROUTINE MeshTopologyNodesFinalise(nodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: nodes !<A pointer to the mesh topology nodes to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx

    ENTERS("MeshTopologyNodesFinalise",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      IF(ALLOCATED(nodes%nodes)) THEN
        DO nodeIdx=1,SIZE(nodes%nodes,1)
          CALL MeshTopologyNodeFinalise(nodes%nodes(nodeIdx),err,error,*999)
        ENDDO !nodesIdx
        DEALLOCATE(nodes%nodes)
      ENDIF
      IF(ASSOCIATED(nodes%nodesTree)) CALL TREE_DESTROY(nodes%nodesTree,err,error,*999)
      DEALLOCATE(nodes)
    ENDIF
 
    EXITS("MeshTopologyNodesFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyNodesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyNodesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given mesh topology. \todo finalise on errors
  SUBROUTINE MeshTopologyNodesInitialise(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyNodesInitialise",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      IF(ASSOCIATED(topology%nodes)) THEN
        CALL FlagError("Mesh already has topology nodes associated.",err,error,*999)
      ELSE
        ALLOCATE(topology%nodes,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate topology nodes.",err,error,*999)
        topology%nodes%numberOfNodes=0
        topology%nodes%meshComponentTopology=>topology
        NULLIFY(topology%nodes%nodesTree)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyNodesInitialise")
    RETURN
999 ERRORSEXITS("MeshTopologyNodesInitialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopologyNodesInitialise

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given list of MESHES. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,MESHES,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(MESHES_TYPE), POINTER :: MESHES !<The list of meshes containing the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx

    ENTERS("MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        mesh_idx=1
        DO WHILE(mesh_idx<=MESHES%NUMBER_OF_MESHES.AND..NOT.ASSOCIATED(MESH))
          IF(MESHES%MESHES(mesh_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            MESH=>MESHES%MESHES(mesh_idx)%PTR
          ELSE
            mesh_idx=mesh_idx+1
          ENDIF
        ENDDO
      ENDIF
    ELSE
      CALL FlagError("Meshes is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_USER_NUMBER_FIND_GENERIC")
    RETURN
999 ERRORSEXITS("MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_USER_NUMBER_FIND_GENERIC

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given INTERFACE. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_INTERFACE(USER_NUMBER,INTERFACE,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    ENTERS("MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_USER_NUMBER_FIND_INTERFACE")
    RETURN
999 ERRORSEXITS("MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MESH_USER_NUMBER_FIND_INTERFACE

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given REGION. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_REGION(USER_NUMBER,REGION,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: REGION !<The region containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("MESH_USER_NUMBER_FIND_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_USER_NUMBER_FIND_REGION")
    RETURN
999 ERRORSEXITS("MESH_USER_NUMBER_FIND_REGION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_USER_NUMBER_FIND_REGION

  !
  !================================================================================================================================
  !

  !>Finalises the meshes and deallocates all memory
  SUBROUTINE MESHES_FINALISE(MESHES,ERR,ERROR,*)

   !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to finalise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH
 
    ENTERS("MESHES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      DO WHILE(MESHES%NUMBER_OF_MESHES>0)
        MESH=>MESHES%MESHES(1)%PTR
        CALL MESH_DESTROY(MESH,ERR,ERROR,*999)
      ENDDO !mesh_idx
      DEALLOCATE(MESHES)
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("MESHES_FINALISE")
    RETURN
999 ERRORSEXITS("MESHES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESHES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the generic meshes.
  SUBROUTINE MESHES_INITIALISE_GENERIC(MESHES,ERR,ERROR,*)

    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("MESHES_INITIALISE_GENERIC",ERR,ERROR,*998)

    IF(ASSOCIATED(MESHES)) THEN
      CALL FlagError("Meshes is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(MESHES,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Meshes could not be allocated",ERR,ERROR,*999)
      NULLIFY(MESHES%REGION)
      NULLIFY(MESHES%INTERFACE)
      MESHES%NUMBER_OF_MESHES=0
      NULLIFY(MESHES%MESHES)
    ENDIF
    
    EXITS("MESHES_INITIALISE_GENERIC")
    RETURN
999 CALL MESHES_FINALISE(MESHES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("MESHES_INITIALISE_GENERIC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_GENERIC

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given interface.
  SUBROUTINE MESHES_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESHES)) THEN
        LOCAL_ERROR="Interface number "//TRIM(NumberToVString(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(INTERFACE%MESHES,ERR,ERROR,*999)
        INTERFACE%MESHES%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_INTERFACE")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_INTERFACE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_INTERFACE

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given region.
  SUBROUTINE MESHES_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NumberToVString(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(REGION%MESHES,ERR,ERROR,*999)
        REGION%MESHES%REGION=>REGION
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_REGION")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_REGION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_REGION

  !
  !================================================================================================================================
  !

    !>Gets the domain for a given node in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OPENCMISS::CMISSDecompositionNodeDomainGet
  SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,USER_NODE_NUMBER,MESH_COMPONENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables`
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_NODE_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MeshNodesType), POINTER :: MESH_NODES
    TYPE(DOMAIN_TYPE), POINTER :: MESH_DOMAIN

    ENTERS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    GLOBAL_NODE_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_NODES=>MESH_TOPOLOGY%nodes
            IF(ASSOCIATED(MESH_NODES)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_NODES%nodesTree,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_NODES%nodesTree,TREE_NODE,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_NODE_NUMBER>0.AND.GLOBAL_NODE_NUMBER<=MESH_TOPOLOGY%NODES%numberOfNodes) THEN
                  IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
                    MESH_DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
                    IF(ASSOCIATED(MESH_DOMAIN)) THEN
                      DOMAIN_NUMBER=MESH_DOMAIN%NODE_DOMAIN(GLOBAL_NODE_NUMBER)
                    ELSE
                      CALL FlagError("Decomposition domain is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Mesh Component number "//TRIM(NumberToVString(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. The limits are 1 to "// &
                      & TRIM(NumberToVString(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NumberToVString(GLOBAL_NODE_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NumberToVString(MESH_TOPOLOGY%NODES%numberOfNodes,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Decomposition mesh node corresponding to user number not found.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Decomposition mesh nodes are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("DECOMPOSITION_NODE_DOMAIN_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET

  !
  !================================================================================================================================
  !

  !>Initialises the embedded meshes.
  SUBROUTINE EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*)

    !Argument variables
    !TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to initialise
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EMBEDDED_MESH_INITIALISE",ERR,ERROR,*998)
    
    ALLOCATE(MESH_EMBEDDING,STAT=ERR)
    NULLIFY(MESH_EMBEDDING%PARENT_MESH)
    NULLIFY(MESH_EMBEDDING%CHILD_MESH)
    
    EXITS("EMBEDDED_MESH_INITIALISE")
    RETURN
!999 CALL EMBEDDED_MESH_FINALISE(MESH_EMBEDDING,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EMBEDDED_MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EMBEDDED_MESH_INITIALISE

  !
  !================================================================================================================================
  !

  !>Creates an embedding of one mesh in another
  SUBROUTINE MESH_EMBEDDING_CREATE(MESH_EMBEDDING, PARENT_MESH, CHILD_MESH,ERR,ERROR,*)
!    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: PARENT_MESH !<The parent mesh
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: CHILD_MESH  !<The child mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: NGP = 0, ne

    ENTERS("MESH_EMBEDDING_CREATE",ERR,ERROR,*999) 
    
    WRITE(*,*) 'parent mesh', PARENT_MESH%NUMBER_OF_ELEMENTS
    WRITE(*,*) 'child mesh', child_MESH%NUMBER_OF_ELEMENTS
    CALL EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*999)

    DO ne=1,PARENT_MESH%NUMBER_OF_ELEMENTS
      NGP = MAX(NGP,PARENT_MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%&
        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS)
    ENDDO !ne

    MESH_EMBEDDING%PARENT_MESH => PARENT_MESH
    MESH_EMBEDDING%CHILD_MESH  => CHILD_MESH
    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate child node positions.",ERR,ERROR,*999)
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(NGP,PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate gauss point positions.",ERR,ERROR,*999)
    
    EXITS("MESH_EMBEDDING_CREATE")
    RETURN 

999 ERRORSEXITS("MESH_EMBEDDING_CREATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_CREATE

  !
  !================================================================================================================================
  !

  !>Sets the positions of nodes in the child mesh for one element in the parent mesh
  SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION(MESH_EMBEDDING, ELEMENT_NUMBER, NODE_NUMBERS, XI_COORDS,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER  !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBERS(:) !<NODE_NUMBERS(node_idx) Node numbers in child mesh for the node_idx'th embedded node in the ELEMENT_NUMBER'th element of the parent mesh
    REAL(DP), INTENT(IN)      :: XI_COORDS(:,:)  !<XI_COORDS(:,node_idx) Xi coordinates of the node_idx'th embedded node in the ELEMENT_NUMBER'th

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR,*999)

    IF(ELEMENT_NUMBER<1 .OR. ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Element number out of range",ERR,ERROR,*999)
    ENDIF

    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NUMBER_OF_NODES = SIZE(NODE_NUMBERS)

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(SIZE(NODE_NUMBERS)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(1:SIZE(NODE_NUMBERS)) = NODE_NUMBERS(1:SIZE(NODE_NUMBERS))

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(SIZE(XI_COORDS,1),SIZE(XI_COORDS,2)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2)) = &
      & XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2))
    
    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  !
  !================================================================================================================================
  !

  !>Sets the positions of a Gauss point of the parent mesh in terms of element/xi coordinate in the child mesh
  SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA(MESH_EMBEDDING, PARENT_ELEMENT_NUMBER, GAUSSPT_NUMBER,&
    & PARENT_XI_COORD, CHILD_ELEMENT_NUMBER, CHILD_XI_COORD,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING   !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: PARENT_ELEMENT_NUMBER           !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: GAUSSPT_NUMBER                  !<Gauss point number in this element
    REAL(DP), INTENT(IN) :: PARENT_XI_COORD(:)              !<Xi coordinate in parent element

    INTEGER(INTG), INTENT(IN) :: CHILD_ELEMENT_NUMBER !<Element number in the child mesh
    REAL(DP), INTENT(IN) :: CHILD_XI_COORD(:)    !<Xi coordinate in child element

    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string

    ENTERS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR,*999)

    IF(PARENT_ELEMENT_NUMBER<1 .OR. PARENT_ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Parent element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(CHILD_ELEMENT_NUMBER<1 .OR. CHILD_ELEMENT_NUMBER > MESH_EMBEDDING%CHILD_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Child element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(GAUSSPT_NUMBER<1 .OR. GAUSSPT_NUMBER > SIZE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION,1)) THEN
      CALL FlagError("Gauss point number out of range",ERR,ERROR,*999)
    ENDIF

    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %PARENT_XI_COORD(SIZE(PARENT_XI_COORD)))
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %CHILD_XI_COORD(SIZE(CHILD_XI_COORD)))


    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD)) = &
      & PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD)) = &
      & CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%ELEMENT_NUMBER = CHILD_ELEMENT_NUMBER

    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA


  !
  !================================================================================================================================
  !

  !> Find the mesh with the given user number, or throw an error if it does not exist.
  SUBROUTINE MESH_USER_NUMBER_TO_MESH( USER_NUMBER, REGION, MESH, ERR, ERROR, * )
    !Arguments
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: REGION !<The region containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the mesh with the specified user number.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_USER_NUMBER_TO_MESH", ERR, ERROR, *999 )

    NULLIFY( MESH )
    CALL MESH_USER_NUMBER_FIND( USER_NUMBER, REGION, MESH, ERR, ERROR, *999 )
    IF( .NOT.ASSOCIATED( MESH ) ) THEN
      LOCAL_ERROR = "A mesh with an user number of "//TRIM(NumberToVString( USER_NUMBER, "*", ERR, ERROR ))// &
        & " does not exist on region number "//TRIM(NumberToVString( REGION%USER_NUMBER, "*", ERR, ERROR ))//"."
      CALL FlagError( LOCAL_ERROR, ERR, ERROR, *999 )
    ENDIF

    EXITS( "MESH_USER_NUMBER_TO_MESH" )
    RETURN
999 ERRORSEXITS( "MESH_USER_NUMBER_TO_MESH", ERR, ERROR )
    RETURN 1

  END SUBROUTINE MESH_USER_NUMBER_TO_MESH

  !
  !================================================================================================================================
  !
!!\todo THIS SHOULD REALLY BE MESH_USER_NUMBER_TO_DECOMPOSITION
  
  !> Find the decomposition with the given user number, or throw an error if it does not exist.
  SUBROUTINE DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION( USER_NUMBER, MESH, DECOMPOSITION, ERR, ERROR, * )
    !Arguments
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to find
    TYPE(MESH_TYPE), POINTER :: MESH !<The mesh containing the decomposition
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On exit, a pointer to the decomposition with the specified user number.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION", ERR, ERROR, *999 )

    NULLIFY( DECOMPOSITION )
    CALL DECOMPOSITION_USER_NUMBER_FIND( USER_NUMBER, MESH, DECOMPOSITION, ERR, ERROR, *999 )
    IF( .NOT.ASSOCIATED( DECOMPOSITION ) ) THEN
      LOCAL_ERROR = "A decomposition with an user number of "//TRIM(NumberToVString( USER_NUMBER, "*", ERR, ERROR ))// &
        & " does not exist on mesh number "//TRIM(NumberToVString( MESH%USER_NUMBER, "*", ERR, ERROR ))//"."
      CALL FlagError( LOCAL_ERROR, ERR, ERROR, *999 )
    ENDIF

    EXITS( "DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION" )
    RETURN
999 ERRORS( "DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION", ERR, ERROR )
    EXITS( "DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION")
    RETURN 1

  END SUBROUTINE DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION

 
  !================================================================================================================================
  !> Following routine sets tpwgts
  !===============================================================================================================================
  SUBROUTINE DECOMPOSITION_TPWGT_SET(DECOMPOSITION, TPWGT, ERR, error, *)
    !Argument
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the vertex domains for.
    REAL(DP), INTENT(IN)	      :: TPWGT(:) !<An array that specifies the distribution of processor load.
    INTEGER(INTG), INTENT(OUT)        :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !local variable
    INTEGER(INTG)                     :: number_computational_nodes

      
    ENTERS("DECOMPOSITION_TPWGT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN

      number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
      IF(ERR/=0) GOTO 999

      IF(ALLOCATED(DECOMPOSITION%TPWGT)) DEALLOCATE(DECOMPOSITION%TPWGT)

      ALLOCATE(DECOMPOSITION%TPWGT(SIZE(TPWGT,1)),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate tpwgt.",ERR,ERROR,*999)
      DECOMPOSITION%TPWGT = TPWGT

    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    END IF
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TPWGT_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TPWGT_SET


  !================================================================================================================================
  !> Following routine activates vertex-wise decomposition. 
  !===============================================================================================================================


  SUBROUTINE DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET(DECOMPOSITION, NODE_BASED_DECOMPOSITION, ERR, error, *)
    !Argument
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    LOGICAL, INTENT(IN)	              :: NODE_BASED_DECOMPOSITION
    INTEGER(INTG), INTENT(OUT)        :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

      
    ENTERS("DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN

      DECOMPOSITION%NODE_BASED_DECOMPOSITION = NODE_BASED_DECOMPOSITION

    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    END IF
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET

 !==================================================================================================================================
 !> Following routine sets ubvec parameter for vertex-wise decomposition
 !==================================================================================================================================

  SUBROUTINE DECOMPOSITION_UBVEC_SET(DECOMPOSITION, UBVEC, ERR, error, *)
    !Argument
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the vertex domains for.
    REAL(DP), INTENT(IN)	      :: UBVEC(:) ! SPecifies error in vertex load distribution. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

 
    ENTERS("DECOMPOSITION_UBVEC_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
       
      IF(ALLOCATED(DECOMPOSITION%UBVEC)) DEALLOCATE(DECOMPOSITION%UBVEC) ! In case UBVEC has been previously allocated before.

      ALLOCATE(DECOMPOSITION%UBVEC(SIZE(UBVEC)), STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate UBVEC.",ERR,ERROR,*999)

      DECOMPOSITION%UBVEC = UBVEC

    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    END IF
    RETURN
999 ERRORSEXITS("DECOMPOSITION_UBVEC_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_UBVEC_SET

  !==================================================================================================================================
  !> Following routine sets number of constraints  for decomposition
  !==================================================================================================================================

  SUBROUTINE DECOMPOSITION_NUMBER_OF_CONSTRAINTS_SET(DECOMPOSITION, NUMBER_OF_CONSTRAINTS, ERR, error, *)
    !Argument
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(IN)	      :: NUMBER_OF_CONSTRAINTS !< Number of vertex load constraints to balance. 
    INTEGER(INTG), INTENT(OUT)        :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

 
    ENTERS("DECOMPOSITION_NUMBER_OF_CONSTRAINTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN

      DECOMPOSITION%NUMBER_OF_CONSTRAINTS = NUMBER_OF_CONSTRAINTS

    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    END IF
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_CONSTRAINTS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_CONSTRAINTS_SET



  !==================================================================================================================================
  !> Following routine sets node weights on a coupled mesh graph G_i.
  !========================================================================================================================
  SUBROUTINE DECOMPOSITION_NODE_WEIGHT_SET(DECOMPOSITION, NODE_WEIGHT_SET, ERR, error, *)
    !Argument
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(IN)	      :: NODE_WEIGHT_SET(:)
    INTEGER(INTG), INTENT(OUT)        :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

 
    ENTERS("DECOMPOSITION_NODE_WEIGHT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN


      IF(ALLOCATED(DECOMPOSITION%NODE_WEIGHT_SET)) DEALLOCATE(DECOMPOSITION%NODE_WEIGHT_SET) 
      ALLOCATE(DECOMPOSITION%NODE_WEIGHT_SET(SIZE(NODE_WEIGHT_SET,1)),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate NODE_WEIGHT_SET.",ERR,ERROR,*999)
      DECOMPOSITION%NODE_WEIGHT_SET = NODE_WEIGHT_SET

    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    END IF
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NODE_WEIGHT_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NODE_WEIGHT_SET

!================================================================================================================
! Calculate graph information for a mesh.
!================================================================================================================
  SUBROUTINE GET_SURROUNDING_NODES(MESH, ERR, ERROR, *)
    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    TYPE(BASIS_TYPE), POINTER                :: BASIS
 
    ENTERS("GET_SURROUNDING_NODES",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN

      BASIS=>MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(1)%BASIS

      IF(ASSOCIATED(BASIS)) THEN

        IF(BASIS%NUMBER_OF_NODES==2_INTG) THEN 
  
          CALL GET_SURROUNDING_NODES_LINEAR_LINE_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==3_INTG .AND. BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE & 
          & .AND. BASIS%INTERPOLATION_XI(1)==BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN 

          CALL GET_SURROUNDING_NODES_QUADRATIC_LINE_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==4_INTG .AND. BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE & 
          & .AND. BASIS%INTERPOLATION_XI(1)==BASIS_CUBIC_LAGRANGE_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_CUBIC_LINE_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==3_INTG .AND. BASIS%TYPE==BASIS_SIMPLEX_TYPE .AND. & 
            &  BASIS%INTERPOLATION_XI(1)==BASIS_LINEAR_SIMPLEX_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_LINEAR_TRIANGULAR_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==6_INTG .AND. BASIS%TYPE==BASIS_SIMPLEX_TYPE .AND. & 
          &  BASIS%INTERPOLATION_XI(1)==BASIS_QUADRATIC_SIMPLEX_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_QUADRATIC_TRIANGULAR_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==4_INTG .AND. BASIS%TYPE==BASIS_SIMPLEX_TYPE .AND. & 
          &  BASIS%INTERPOLATION_XI(1)==BASIS_LINEAR_SIMPLEX_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_LINEAR_TETRAHEDRAL_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==10_INTG .AND. BASIS%TYPE==BASIS_SIMPLEX_TYPE .AND. & 
          &  BASIS%INTERPOLATION_XI(1)==BASIS_QUADRATIC_SIMPLEX_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_QUADRATIC_TETRAHEDRAL_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==4_INTG .AND. BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE & 
          & .AND. BASIS%INTERPOLATION_XI(1)==BASIS_LINEAR_LAGRANGE_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_LINEAR_QUADRILATERAL_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==9_INTG .AND. BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE & 
          & .AND. BASIS%INTERPOLATION_XI(1)==BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_QUADRATIC_QUADRILATERAL_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==16_INTG .AND. BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE & 
          & .AND. BASIS%INTERPOLATION_XI(1)==BASIS_CUBIC_LAGRANGE_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_CUBIC_QUADRILATERAL_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==8_INTG .AND. BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE & 
          & .AND. BASIS%INTERPOLATION_XI(1)==BASIS_LINEAR_LAGRANGE_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_LINEAR_HEXAHEDRAL_MESH(MESH, ERR, error, *999)

        ELSE IF(BASIS%NUMBER_OF_NODES==27_INTG .AND. BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE & 
          & .AND. BASIS%INTERPOLATION_XI(1)==BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN

          CALL GET_SURROUNDING_NODES_QUADRATIC_HEXAHEDRAL_MESH(MESH, ERR, error, *999)

        END IF 

      ELSE
        CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
      END IF      
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    END IF
    RETURN
    EXITS("GET_SURROUNDING_NODES")
    RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES",ERR,ERROR)
999 RETURN 1 
  END SUBROUTINE GET_SURROUNDING_NODES
 !=====================================================================================================

  SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_QUADRILATERAL_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_QUADRATIC_QUADRILATERAL_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)

      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 9

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 4) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              ELSE IF(NOdeIdx == 5) THEN 
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)             

              ELSE IF(NOdeIdx == 6) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)               

              ELSE IF(NOdeIdx == 7) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 8) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 9) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnodee

    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
    EXITS("GET_SURROUNDING_NODES_QUADRATIC_QUADRILATERAL_MESH.")
    RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_QUADRATIC_QUADRILATERAL_MESH.",ERR,ERROR)
999 RETURN 1 

  END SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_QUADRILATERAL_MESH
 !====================================================================================================

  SUBROUTINE GET_SURROUNDING_NODES_CUBIC_QUADRILATERAL_MESH(MESH,ERR,error,*)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_CUBIC_QUADRILATERAL_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh.
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)
      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)

      DO MeshNOde = 1, MeshnumberOfNodes

        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 16

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1) THEN 
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 4) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 5) THEN 
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)                

              ELSE IF(NOdeIdx == 6) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 7) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 8) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(12)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 9) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(13)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 10) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(10)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 11) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(11)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 12) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(16)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 13) THEN 
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(14)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)                

              ELSE IF(NOdeIdx == 14) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(13)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
            
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(15)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              
              ELSE IF(NOdeIdx == 15) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(14)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(16)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 16) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(15)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(12)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
             
              END  IF

            END  IF !NOdeIdx

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnodee
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
    EXITS("GET_SURROUNDING_NODES_CUBIC_QUADRILATERAL_MESH.")
    RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_CUBIC_QUADRILATERAL_MESH.",ERR,ERROR)
999 RETURN 1 

  END SUBROUTINE GET_SURROUNDING_NODES_CUBIC_QUADRILATERAL_MESH

 !====================================================================================================

  SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_HEXAHEDRAL_MESH(MESH,ERR,error,*)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_QUADRATIC_HEXAHEDRAL_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh.
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)
      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)

      DO MeshNOde = 1, MeshnumberOfNodes

        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 27

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(10)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(12)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 4) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 5) THEN 
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
           
              ELSE IF(NOdeIdx == 6) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 7) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(16)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 8) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 9) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(18)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 10) THEN 


                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(19)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 11) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(11)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 12) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(21)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 13) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(13)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              ELSE IF(NOdeIdx == 14) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(14)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
            
              
              ELSE IF(NOdeIdx == 15) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(15)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 16) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(25)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 17) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(17)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 18) THEN 


                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(27)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 19) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(20)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(22)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(10)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 20) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(19)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(21)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 21) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(12)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(20)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(24)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
           
              ELSE IF(NOdeIdx == 22) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(19)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(25)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 23) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(23)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 24) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(21)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(27)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 25) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(22)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(26)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(16)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 26) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(25)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(27)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 27) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(24)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(26)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(18)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

             
              END IF

            END IF !NOdeIdx

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)


      END DO ! meshnodee
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
    EXITS("GET_SURROUNDING_NODES_QUADRATIC_HEXAHEDRAL_MESH.")
    RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_QUADRATIC_HEXAHEDRAL_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_HEXAHEDRAL_MESH
 !=====================================================================================================

  SUBROUTINE GET_SURROUNDING_NODES_LINEAR_HEXAHEDRAL_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_LINEAR_HEXAHEDRAL_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)

      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 8

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

 
              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              ELSE IF(NOdeIdx == 4) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              
              
                EXIT

              ELSE IF(NOdeIdx == 5) THEN 
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              ELSE IF(NOdeIdx == 6) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
                
              ELSE IF(NOdeIdx == 7) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              ELSE IF(NOdeIdx == 8) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              
              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnodee
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
    EXITS("GET_SURROUNDING_NODES_LINEAR_HEXAHEDRAL_MESH.")
    RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_LINEAR_HEXAHEDRAL_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_LINEAR_HEXAHEDRAL_MESH
 !====================================================================================================
 SUBROUTINE GET_SURROUNDING_NODES_CUBIC_LINE_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !!<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_LINEAR_LINE_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)

      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 4
 
            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 4) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              END  IF


            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)


        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnode
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
   EXITS("GET_SURROUNDING_NODES_CUBIC_LINE_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_CUBIC_LINE_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_CUBIC_LINE_MESH
 !=====================================================================================================
  SUBROUTINE GET_SURROUNDING_NODES_LINEAR_LINE_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_LINEAR_LINE_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)


      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 2
 
            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              
 
              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnode
      ELSE
        CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
      END IF
   EXITS("GET_SURROUNDING_NODES_LINEAR_LINE_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_LINEAR_LINE_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_LINEAR_LINE_MESH

  !======================================================================================================================================

  SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_LINE_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_QUADRATIC_LINE_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)

      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)

      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 3
 
            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
               
              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)


        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnode
      RETURN
      ELSE
        CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
      END IF
   EXITS("GET_SURROUNDING_NODES_QUADRATIC_LINE_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_QUADRATIC_LINE_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_LINE_MESH


 !====================================================================================================

  SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_TRIANGULAR_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)

    ENTERS("GET_SURROUNDING_NODES_QUADRATIC_TRIANGULAR_MESH",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES

      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)


      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          DO NOdeIdx = 1, 6
 
            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              
                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
 
              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)



              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 4) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)



              ELSE IF(NOdeIdx == 5) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)


              ELSE IF(NOdeIdx == 6) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)


      END DO ! meshnode
      ELSE
        CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
      END IF
   EXITS("GET_SURROUNDING_NODES_QUADRATIC_TRIANGULAR_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_QUADRATIC_TRIANGULAR_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_TRIANGULAR_MESH


 !====================================================================================================
  SUBROUTINE GET_SURROUNDING_NODES_LINEAR_QUADRILATERAL_MESH(MESH,ERR,ERROR,*)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)



    ENTERS("GET_SURROUNDING_NODES_LINEAR_QUADRILATERAL_MESH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES
      PROCID=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,ERR,ERROR,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)
      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST.",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes

        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes,ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          Flag = .FALSE.

          DO NOdeIdx = 1, 4

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) THEN
        
              IF(NOdeIdx == 1 .OR. NOdeIdx == 4) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

 
              ELSE IF(NOdeIdx == 2 .OR. NOdeIdx == 3 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              
              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)
 
        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        
        DEALLOCATE(SURROUNDING_NODES)
      END DO ! meshnode 
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
   EXITS("GET_SURROUNDING_NODES_LINEAR_QUADRILATERAL_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_LINEAR_QUADRILATERAL_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_LINEAR_QUADRILATERAL_MESH

  !======================================================================================================================================
   SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_TETRAHEDRAL_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)


    ENTERS("GET_SURROUNDING_NODES_QUADRATIC_TETRAHEDRAL_MESH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN

      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES
      PROCID=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)

      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          Flag = .FALSE.

          DO NOdeIdx = 1, 10

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

 
              ELSE IF(NOdeIdx == 2 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(5)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(10)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)              
              
              ELSE IF(NOdeIdx == 3 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(8)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(6)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)   
              
              
              ELSE IF(NOdeIdx == 4 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(9)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(7)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(10)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)   
              
              ELSE IF(NOdeIdx == 5 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
           
              
              ELSE IF(NOdeIdx == 6 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
  
                           
              ELSE IF(NOdeIdx == 7 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              ELSE IF(NOdeIdx == 8 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
             
              
              ELSE IF(NOdeIdx == 9 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
 
                      
              ELSE IF(NOdeIdx == 10 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

              
              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnode
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
   EXITS("GET_SURROUNDING_NODES_QUADRATIC_TETRAHEDRAL_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_QUADRATIC_TETRAHEDRAL_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_QUADRATIC_TETRAHEDRAL_MESH


  !======================================================================================================================================
   SUBROUTINE GET_SURROUNDING_NODES_LINEAR_TETRAHEDRAL_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)



    ENTERS("GET_SURROUNDING_NODES_LINEAR_TETRAHEDRAL_MESH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES
      PROCID=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)

      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          Flag = .FALSE.

          DO NOdeIdx = 1, 4

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              
 
              ELSE IF(NOdeIdx == 2 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)              
              
              ELSE IF(NOdeIdx == 3 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(4)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)   
              
              ELSE IF(NOdeIdx == 4 ) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)   
              
              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnode

    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
   EXITS("GET_SURROUNDING_NODES_LINEAR_TETRAHEDRAL_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_LINEAR_TETRAHEDRAL_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_LINEAR_TETRAHEDRAL_MESH

!====================================================================================================
  SUBROUTINE GET_SURROUNDING_NODES_LINEAR_TRIANGULAR_MESH(MESH, ERR, error, *)

    ! Arguments
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT)  :: MESH !<A pointer to the mesh to calculate grap information for.
    INTEGER(INTG), INTENT(OUT)               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)        :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG)                            :: AdjacentNode,ElementIdx,ElementNUmberOfNOdes,MeshNUmberOfNOdes, &
                                                & MeshNOde,NOdeIdx,Node,NUMBER_OF_SURROUNDING_NODES,PROCID
    TYPE(MeshNodesType), pointer             :: meshNOdes
    INTEGER(INTG), ALLOCATABLE               :: ElementNodes(:,:), SURROUNDING_NODES(:)
    LOGICAL                                  :: Flag
    TYPE(LIST_PTR_TYPE), ALLOCATABLE         :: SURROUNDING_NODE_LIST(:)



    ENTERS("GET_SURROUNDING_NODES_LINEAR_TRIANGULAR_MESH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      meshNOdes=>MESH%TOPOLOGY(1)%PTR%NODES
      PROCID=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      ! Total nodes in a mesh
      CALL MeshTopologyNodesNumberOfNodesGet(meshNodes,MeshnumberOfNodes,err,error,*999)

      ALLOCATE(SURROUNDING_NODE_LIST(MeshnumberOfNodes), STAT=ERR)

      IF(ERR/=0)  CALL FlagError("Could not allocate SURROUNDING_NODE_LIST .",ERR,ERROR,*999)


      DO MeshNOde = 1, MeshnumberOfNodes


        NULLIFY(SURROUNDING_NODE_LIST(MeshNOde)%PTR)
        CALL LIST_CREATE_START(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(SURROUNDING_NODE_LIST(MeshNOde)%PTR,MeshnumberOfNodes, ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)

        DO ElementIdx = 1, SIZE(MESH%TOPOLOGY(1)%ptr%Elements%Elements,1)

          Flag = .FALSE.

          DO NOdeIdx = 1, 3

            IF(MESH%TOPOLOGY(1)%ptr%Elements%Elements(ElementIdx)%GLOBAL_ELEMENT_NODES(NOdeIdx) == MeshNOde) then
        
              IF(NOdeIdx == 1) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
 
              ELSE IF(NOdeIdx == 2) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(3)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              ELSE IF(NOdeIdx == 3) THEN 

                AdjacentNode = MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(1)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)

                AdjacentNode= MESH%TOPOLOGY(1)%ptr%Elements%Elements(ELementIdx)%GLOBAL_ELEMENT_NODES(2)
                CALL LIST_ITEM_ADD(SURROUNDING_NODE_LIST(MeshNOde)%PTR,AdjacentNode,ERR,ERROR,*999)
              

              END  IF

            END  IF 

          END DO!DO column = 1, SIZE(CoupledMesh%ConnectivityMatrix,2)
          
        END  DO ! DO row = 1, SIZE(CoupledMesh%ConnectivityMatrix,1)    

        CALL LIST_REMOVE_DUPLICATES(SURROUNDING_NODE_LIST(MeshNOde)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(SURROUNDING_NODE_LIST(MeshNOde)%PTR, &
          & NUMBER_OF_SURROUNDING_NODES, SURROUNDING_NODES, ERR,ERROR,*999)

        IF(ALLOCATED(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)) &
          & DEALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes)

        ALLOCATE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes & 
            & (NUMBER_OF_SURROUNDING_NODES), STAT= ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate surrounding node data structure.",ERR,ERROR,*999)

        MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(MeshNOde)%surroundingNodes(1:NUMBER_OF_SURROUNDING_NODES) = &
          & SURROUNDING_NODES(1:NUMBER_OF_SURROUNDING_NODES)

        DEALLOCATE(SURROUNDING_NODES)

      END DO ! meshnode

    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    END IF
   EXITS("GET_SURROUNDING_NODES_LINEAR_TRIANGULAR_MESH.")
   RETURN
!999 ERRORSEXITS("GET_SURROUNDING_NODES_LINEAR_TRIANGULAR_MESH.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GET_SURROUNDING_NODES_LINEAR_TRIANGULAR_MESH

!========================================================================================================================
! This subroutine builds Adjncy and xadj arrays for vertex-based decomposition  
!========================================================================================================================
 SUBROUTINE GET_ADJNCY_AND_XADJ(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the node domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    INTEGER(INTG)              :: adjncy_size, counter, PROC_ID, vtx_dist_idx, xadj_counter 
    TYPE(MESH_TYPE), POINTER   :: MESH

    ENTERS("GET_ADJNCY_AND_XADJ",ERR,ERROR,*999) 

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN

          PROC_ID=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) ! fetch processor id

          IF(ALLOCATED(DECOMPOSITION%XADJ)) DEALLOCATE(DECOMPOSITION%XADJ)
          ALLOCATE(DECOMPOSITION%XADJ(0:(DECOMPOSITION%VTX_DIST(PROC_ID+1)-DECOMPOSITION%VTX_DIST(PROC_ID) ) ))

          DECOMPOSITION%XADJ(0) = 0 
          xadj_counter= 0 
          adjncy_size = 0  

          DO vtx_dist_idx = DECOMPOSITION%VTX_DIST(PROC_ID),DECOMPOSITION%VTX_DIST(PROC_ID+1)-1
    
            xadj_counter= xadj_counter+ 1 
       
            adjncy_size = adjncy_size + SIZE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes,1)
 
            DECOMPOSITION%XADJ(xadj_counter)      = adjncy_size 

          END DO !vtx_dist_idx

          IF(ALLOCATED(DECOMPOSITION%ADJNCY)) DEALLOCATE(DECOMPOSITION%ADJNCY)
          ALLOCATE(DECOMPOSITION%ADJNCY(0:adjncy_size-1), STAT= ERR)
          IF(ERR/=0) CALL FlagError("Could not ADJNCY array.",ERR,ERROR,*999)

          counter= 0 
          DO vtx_dist_idx = DECOMPOSITION%VTX_DIST(PROC_ID),DECOMPOSITION%VTX_DIST(PROC_ID+1)-1

            DECOMPOSITION%ADJNCY( & 
              & counter: COUNTER+SIZE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes)-1)= & 
                & MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes
 
            counter= counter+ SIZE(MESH%TOPOLOGY(1)%ptr%Nodes%Nodes(vtx_dist_idx+1)%surroundingNodes,1) 

          END DO !vtx_dist_idx
          DECOMPOSITION%ADJNCY = DECOMPOSITION%ADJNCY - 1 ! to make numbering start from 0 

        ELSE
          CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("GET_ADJNCY_AND_XADJ.")
   RETURN
999 ERRORSEXITS("GET_ADJNCY_AND_XADJ.",ERR,ERROR)
   RETURN 1

 END SUBROUTINE GET_ADJNCY_AND_XADJ


! ===========================================================================================================================
! The following subroutine calculates the parameters related to vertex-based decomposition.
! ===========================================================================================================================
  SUBROUTINE DECOMPOSITION_NODE_BASED_DECOMPOSITION_PARAMETERS_SET(DECOMPOSITION,ERR,ERROR,*)

   !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(INOUT) :: DECOMPOSITION !<A pointer to the decomposition to calculate the node domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    REAL(DP), ALLOCATABLE      ::  SUM_OF_TPWGT(:)
    INTEGER(INTG) :: constraint_idx, PROC_ID, tpwgt_idx 
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_NODE_BASED_DECOMPOSITION_PARAMETERS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        NULLIFY(MESH)
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN

          PROC_ID=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) ! Fetches the processor id.

          DECOMPOSITION%WEIGHT_FLAG = 3_INTG ! Edge and vertex eights activated.
          DECOMPOSITION%NUM_FLAG    = 0_INTG ! Node numbering starts from 0.

          ! As initial partitioning assigning the whole domain to processor 0.          
          DECOMPOSITION%NODE_DOMAIN = 0.0_INTG 
 
          IF(DECOMPOSITION%NUMBER_OF_CONSTRAINTS > 1) THEN ! Which means the user defined value is to be used now. 
         
            DEALLOCATE(DECOMPOSITION%TPWGT)
            DEALLOCATE(DECOMPOSITION%UBVEC)

            ALLOCATE(DECOMPOSITION%TPWGT(DECOMPOSITION%NUMBER_OF_CONSTRAINTS*DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
            IF(ERR/=0)  CALL FlagError("Could not allocate DECOMPOSITION%TPWGT data structure",ERR,ERROR,*999) 

            ALLOCATE(DECOMPOSITION%UBVEC(DECOMPOSITION%NUMBER_OF_CONSTRAINTS), STAT=ERR)
            IF(ERR/=0)  CALL FlagError("Could not allocate DECOMPOSITION%UBVEC data structure",ERR,ERROR,*999) 

            DECOMPOSITION%TPWGT = REAL(1./DECOMPOSITION%NUMBER_OF_DOMAINS,RP)
            DECOMPOSITION%UBVEC = 1.000001_RP

          END IF 

          ! Check if the sum of tpwgts for each constrain is equal to 1. 
          ALLOCATE(SUM_OF_TPWGT(DECOMPOSITION%NUMBER_OF_CONSTRAINTS), STAT=ERR)
          IF(ERR/=0)  CALL FlagError("Could not allocate SUM_OF_TPWGT array",ERR,ERROR,*999) 
          SUM_OF_TPWGT = 0 

          DO constraint_idx = 1, DECOMPOSITION%NUMBER_OF_CONSTRAINTS

            DO tpwgt_idx = constraint_idx , SIZE(DECOMPOSITION%TPWGT,1), DECOMPOSITION%NUMBER_OF_CONSTRAINTS

               SUM_OF_TPWGT(constraint_idx) = SUM_OF_TPWGT(constraint_idx) + DECOMPOSITION%TPWGT(tpwgt_idx)

            END DO !tpwgt_idx

            IF(ABS(SUM_OF_TPWGT(constraint_idx)-1) > 1E-5) THEN

              LOCAL_ERROR = &
                & " The sum of tpwgts for constraint # "//TRIM(NumberToVString( constraint_idx, "*", ERR, ERROR ))// &
                  & " is not equal to 1. "
              CALL FlagError( LOCAL_ERROR, ERR, ERROR, *999 )

            END IF !IF(ABS(SUM_OF_TPWGT(constraint_idx)-1) > 1E-5) THEN
          END DO !constraint_idx
  

          !Setting values of vtx_dist array
          ALLOCATE(DECOMPOSITION%VTX_DIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
          IF(ERR==0) THEN  ! If not already allocated by the user.

            CALL DECOMPOSITION_SET_DEFAULT_VTX_DIST(DECOMPOSITION,ERR,ERROR,*999)
           
          END IF


          ALLOCATE(DECOMPOSITION%XADJ(0:(DECOMPOSITION%VTX_DIST(PROC_ID+1)-DECOMPOSITION%VTX_DIST(PROC_ID))), STAT=ERR)
          IF(ERR==0) THEN  ! If not already allocated by the user. 

            CALL GET_ADJNCY_AND_XADJ(DECOMPOSITION,ERR,error,*999)
           
          END IF

          ALLOCATE(DECOMPOSITION%NODE_WEIGHT_SET( & 
            & 1:DECOMPOSITION%NUMBER_OF_CONSTRAINTS*(DECOMPOSITION%VTX_DIST(PROC_ID+1) &  
              & - DECOMPOSITION%VTX_DIST(PROC_ID))), STAT=ERR)
          IF(ERR==0)  THEN ! If not already allocated by the user.
            DECOMPOSITION%NODE_WEIGHT_SET   = 1
          END IF 

          ALLOCATE(DECOMPOSITION%ADJWT(SIZE(DECOMPOSITION%ADJNCY,1)), STAT=ERR)
          IF(ERR==0) THEN ! If not already allocated by the user.
           DECOMPOSITION%ADJWT             = 1
          END IF

          DEALLOCATE(SUM_OF_TPWGT)
        ELSE
          CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("DECOMPOSITION_NODE_BASED_DECOMPOSITION_PARAMETERS_SET.")
   RETURN
! 999 ERRORSEXITS("DECOMPOSITION_NODE_BASED_DECOMPOSITION_PARAMETERS_SET.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE DECOMPOSITION_NODE_BASED_DECOMPOSITION_PARAMETERS_SET


! ==========================================================================================================================
! Set the default values of VTX_DIST array.
! ==========================================================================================================================

  SUBROUTINE DECOMPOSITION_SET_DEFAULT_VTX_DIST(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(INOUT) :: DECOMPOSITION !<A pointer to the decomposition to calculate the node domains for.
    INTEGER(INTG), INTENT(OUT)                       :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)                :: ERROR !<The error string
    !Local Variables

    INTEGER(INTG)                                    :: vtx_dist_idx
    TYPE(MESH_TYPE), POINTER                         :: MESH   
 
    ENTERS("DECOMPOSITION_SET_DEFAULT_VTX_DIST",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN

       MESH=>DECOMPOSITION%MESH
 
       !Setting default values of vtx_dist array
       IF(ALLOCATED(DECOMPOSITION%VTX_DIST)) DEALLOCATE(DECOMPOSITION%VTX_DIST)
       ALLOCATE(DECOMPOSITION%VTX_DIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)

       DECOMPOSITION%VTX_DIST(0)                               = 0
       DECOMPOSITION%VTX_DIST(DECOMPOSITION%NUMBER_OF_DOMAINS) = MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes 

       DO vtx_dist_idx = 1, DECOMPOSITION%NUMBER_OF_DOMAINS-1 

          DECOMPOSITION%VTX_DIST(vtx_dist_idx) =  DECOMPOSITION%VTX_DIST(vtx_dist_idx-1) + &
            & INT(MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes/DECOMPOSITION%NUMBER_OF_DOMAINS, INTG)                  

       END DO !vtx_dist_idx
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("DECOMPOSITION_SET_DEFAULT_VTX_DIST.")
   RETURN
999 ERRORSEXITS("DECOMPOSITION_SET_DEFAULT_VTX_DIST.",ERR,ERROR)
   RETURN 1

  END SUBROUTINE DECOMPOSITION_SET_DEFAULT_VTX_DIST

! ==========================================================================================================================
! THis subroutine initializes the VTX_DIST array by user defined value
! ==========================================================================================================================

  SUBROUTINE DECOMPOSITION_VTX_DIST_SET(DECOMPOSITION,VTX_DIST,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(INOUT) :: DECOMPOSITION !<A pointer to the decomposition to calculate the vertex domains for.
    INTEGER, INTENT(IN)                              :: VTX_DIST(:) !< Array that stores distribution of vertex between the processor i.e. VTX_DIST(2)-VTX_DIST(1) shows vertex Ids initially distributed to the processor 1. 
    INTEGER(INTG), INTENT(OUT)                       :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)                :: ERROR !<The error string

    ENTERS("DECOMPOSITION_VTX_DIST_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN

       !Setting values of vtx_dist array
       DECOMPOSITION%VTX_DIST = VTX_DIST
 
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("DECOMPOSITION_VTX_DIST_SET.")
   RETURN
999 ERRORSEXITS("DECOMPOSITION_VTX_DIST_SET.",ERR,ERROR)
   RETURN 1

  END SUBROUTINE DECOMPOSITION_VTX_DIST_SET
! ============================================================================================================================
  !>Calculates the vertex domains for a decomposition of a mesh. \see OPENCMISS::CMISSDecompositionNOdeDomainCalculate
  SUBROUTINE DECOMPOSITION_NODE_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the node domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    INTEGER(INTG) :: PARMETIS_OPTIONS(3), PROC_ID, proc_idx_receive, proc_idx_send, STATUS(MPI_STATUS_SIZE)
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_NODE_DOMAIN_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
    
      
      DECOMPOSITION%NODE_DOMAIN = 0 
      IF(DECOMPOSITION%NUMBER_OF_DOMAINS > 1) THEN 
        PROC_ID=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) ! fetch the processor id.


        PARMETIS_OPTIONS(1) = 1
        PARMETIS_OPTIONS(2) = 7      
        PARMETIS_OPTIONS(3) = 99999 

        CALL ParMETIS_V3_PartKway(DECOMPOSITION%VTX_DIST,DECOMPOSITION%XADJ, &
          & DECOMPOSITION%ADJNCY,DECOMPOSITION%NODE_WEIGHT_SET, & 
            & DECOMPOSITION%ADJWT, DECOMPOSITION%WEIGHT_FLAG, & 
              & DECOMPOSITION%NUM_FLAG, DECOMPOSITION%NUMBER_OF_CONSTRAINTS, & 
                & DECOMPOSITION%NUMBER_OF_DOMAINS, DECOMPOSITION%TPWGT, & 
                  & DECOMPOSITION%UBVEC, PARMETIS_OPTIONS, & 
                    & DECOMPOSITION%NUMBER_OF_EDGES_CUT, & 
                      & DECOMPOSITION%NODE_DOMAIN( &
                        & DECOMPOSITION%VTX_DIST(PROC_ID): &
                          & DECOMPOSITION%VTX_DIST(PROC_ID+1)-1), COMPUTATIONAL_ENVIRONMENT%MPI_COMM)


         DO proc_idx_receive = 0, DECOMPOSITION%NUMBER_OF_DOMAINS-1   
       
           IF(PROC_ID /= proc_idx_receive) THEN

             CALL MPI_SEND(DECOMPOSITION%NODE_DOMAIN( &
               & DECOMPOSITION%VTX_DIST(PROC_ID): &
                 & DECOMPOSITION%VTX_DIST(PROC_ID+1)-1), & 
                   & SIZE(DECOMPOSITION%NODE_DOMAIN( &
                     & DECOMPOSITION%VTX_DIST(PROC_ID): &
                       & DECOMPOSITION%VTX_DIST(PROC_ID+1)-1)), & 
                         & MPI_INT, proc_idx_receive, PROC_ID ,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, ERR) 


           END IF
         END DO !!ROC_idx_RECEIVE

         DO proc_idx_send = 0, DECOMPOSITION%NUMBER_OF_DOMAINS-1      
         
           IF(proc_idx_send /= PROC_ID) THEN

              CALL MPI_RECV( &
                & DECOMPOSITION%NODE_DOMAIN( &
                  & DECOMPOSITION%VTX_DIST(proc_idx_send): &
                    & DECOMPOSITION%VTX_DIST(proc_idx_send+1)-1), & 
                      & SIZE(DECOMPOSITION%NODE_DOMAIN( &
                      & DECOMPOSITION%VTX_DIST(proc_idx_send): &
                          & DECOMPOSITION%VTX_DIST(proc_idx_send+1)-1)), & 
                            & MPI_INT, proc_idx_send,proc_idx_send, & 
                              & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,STATUS, Err) 

 
            END IF !proc_idx_send
     

         END DO !procidx_send = 0, number_parts-1
       END IF  !IF(DECOMPOSTION%NUMBER_OF_DOMAINS > 1).  
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    RETURN 
!999 ERRORSEXITS("DECOMPOSITION_NODE_DOMAIN_CALCULATE.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE DECOMPOSITION_NODE_DOMAIN_CALCULATE

 !=====================================================================================================
 
  SUBROUTINE GET_MAXIMUM_EDGE_LENGTH(COUPLED_DECOMPOSITION, MAXIMUM_EDGE_LENGTH, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN)     :: COUPLED_DECOMPOSITION !<A pointer to the decomposition to calculate the node domains for.
    REAL(RP), INTENT(OUT)                     :: MAXIMUM_EDGE_LENGTH   !Maximum edge length in a mesh given by object MESH.
    INTEGER(INTG), INTENT(OUT)                :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)         :: ERROR !<The error string

    !Local Variables
    REAL(RP), ALLOCATABLE                     :: EDGE_LENGTHS(:)
    INTEGER(INTG)                             :: ADJACENCY_NODE, adjacency_idx, adjncy_idx, component_idx, DOMAIN_NUMBER, & 
      & edge_length_idx, MY_COMPUTATIONAL_NODE, node_idx, numberOfComponents 
    REAL(RP)                                  :: coordinate_node_idx(3), coordinate_adjncy_idx(3), VECTOR_ALONG_THE_EDGE(3)
    TYPE(DECOMPOSITION_TYPE), POINTER         :: DECOMPOSITION   
    TYPE(FIELD_VARIABLE_TYPE), POINTER        :: FIELD_VARIABLE 
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
                                         COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(ADJACENCY_NODE,1)   
             VECTOR_ALONG_THE_EDGE(2) =  COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(node_idx,2) - &
                                         COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(ADJACENCY_NODE,2)   
             VECTOR_ALONG_THE_EDGE(3) =  COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(node_idx,3) - &
                                         COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(ADJACENCY_NODE,3)   
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
   ! The following subroutine builds interedges I_{iI} between coupled mesh graphs G_i and the interface graphs G_I
 !================================================================================================================

  SUBROUTINE COUPLED_DECOMPOSITION_GET_INTER_EDGES(COUPLED_DECOMPOSITION, MAXIMUM_INTERFACE_EDGE_LENGTH, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)   :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    REAL(RP), INTENT(IN)                                       :: MAXIMUM_INTERFACE_EDGE_LENGTH !<Maximum edge length in interface graph mesh G^I.   
    INTEGER(INTG), INTENT(OUT)                                 :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                          :: ERROR !<The error string.

    !Local Variables
    TYPE(LIST_TYPE),  POINTER :: BOUNDARY_NODE_LIST, INTER_EDGE_LIST_INTERFACE_MESH, INTER_EDGE_LIST_COUPLED_MESH
    TYPE(LIST_PTR_TYPE), ALLOCATABLE:: ADJNCY_RESTRICTED_TO_INTERFACE_LIST(:)
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION_COUPLED_MESH, DECOMPOSITION_INTERFACE_MESH   
    TYPE(MESH_TYPE), POINTER:: COUPLED_MESH, INTERFACE_MESH 
    TYPE(FIELD_VARIABLE_TYPE), POINTER  :: FIELD_VARIABLE_COUPLED_MESH, FIELD_VARIABLE_INTERFACE_MESH
    INTEGER(INTG) :: boundary_node_idx,COUNTER,coupled_mesh_normal_vector_idx, interface_mesh_normal_vector_idx, & 
      & coupled_mesh_nodes_idx,coupled_mesh_node_idx,component_idx, & 
        & DOMAIN_NUMBER,distance_idx,GLOBAL_NUMBER_OF_INTER_EDGES,INTERFACE_MESH_NODE_TO_ADD,interface_node_idx, & 
          & interface_mesh_node_idx,inter_edge_idx,MY_COMPUTATIONAL_NODE,NUMBER_OF_BOUNDARY_NODES,NUMBER_OF_INTER_EDGES, &
            & NUMBER_OF_NODES,surrounding_node_idx
    REAL(RP) :: coordinate_coupled_mesh_node_idx(3), coordinate_interface_mesh_node_idx(3), DIFFERENCE_IN_ANGLE
    REAL(RP), ALLOCATABLE :: DISTANCE(:), NORMAL_VECTOR_COUPLED_MESH(:,:), NORMAL_VECTOR_INTERFACE_MESH(:,:)				
    INTEGER(INTG), ALLOCATABLE  :: ADJNCY_RESTRICTED_TO_INTERFACE(:,:),BOUNDARY_NODES(:), & 
      & COUPLED_MESH_INTER_EDGE_NODES(:),INTERFACE_MESH_INTER_EDGE_NODES(:),LIST_OF_NODES(:), & 
        & NUMBER_ADJNCY_RESTRICTED_TO_INTERFACE(:) 
    LOGICAL :: INTERFACE_NODE_DETECTED  

   ENTERS("GET_INTER_EDGES",ERR,ERROR,*999)

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
           IF(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(coupled_mesh_node_idx)%boundaryNode .EQV. .TRUE. ) THEN
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
        
           ! to determine the size of COupledMeshNodeAdjancyRestrictedToBoundary
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

             INTERFACE_MESH_NODE_TO_ADD= SQRT(REAL(DOT_PRODUCT(&
               & MINLOC(PACK(DISTANCE,DISTANCE >= -1E-5)-MINVAL(PACK(DISTANCE,DISTANCE >= -1E-5))),&
                 & MINLOC(PACK(DISTANCE,DISTANCE >= -1E-5)-MINVAL(PACK(DISTANCE,DISTANCE >= -1E-5))))))

             ! Calculate normal at vertex coupled_mesh_node_idx
             CALL COUPLED_MESH_CALCULATE_NORMAL(COUPLED_DECOMPOSITION,  & 
               & boundary_node_idx, NORMAL_VECTOR_COUPLED_MESH, ERR, ERROR, *999)

             ! Calculate normal at vertex interface_mesh_node_idx
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

                 IF((DIFFERENCE_IN_ANGLE .GE. 0 .AND. DIFFERENCE_IN_ANGLE .LE. 5) .OR. & 
                   & (DIFFERENCE_IN_ANGLE .GE. 175 .AND. DIFFERENCE_IN_ANGLE .LE. 180)) THEN

                   CALL LIST_ITEM_ADD(INTER_EDGE_LIST_INTERFACE_MESH,INTERFACE_MESH_NODE_TO_ADD, ERR, ERROR, *999)  

                   CALL LIST_ITEM_ADD(INTER_EDGE_LIST_COUPLED_MESH,boundary_node_idx, ERR, ERROR, *999)  

                   EXIT
                    
                 END IF

               END DO !interface_mesh_normal_vector_idx

                 IF((DIFFERENCE_IN_ANGLE .GE. 0 .AND. DIFFERENCE_IN_ANGLE .LE. 5) .OR. & 
                   & (DIFFERENCE_IN_ANGLE .GE. 175 .AND. DIFFERENCE_IN_ANGLE .LE. 180)) EXIT
                
             END DO ! coupled_mesh_normal_vector_idx

             DEALLOCATE(NORMAL_VECTOR_COUPLED_MESH)
             DEALLOCATE(NORMAL_VECTOR_INTERFACE_MESH)

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
! Calculate all possible normal angles on node NODE_ID of the coupled mesh graph G_i
!==================================================================================================================
  SUBROUTINE COUPLED_MESH_CALCULATE_NORMAL(COUPLED_DECOMPOSITION, NODE_ID, NORMAL_VECTOR, ERR, ERROR, *) 

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER ,INTENT(INOUT) :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    INTEGER(INTG), INTENT(OUT)                               :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                        :: ERROR !<The error string.
    INTEGER(INTG)                                            :: NODE_ID !<The global node id to find normal at. 
    REAL(RP), ALLOCATABLE                                    :: NORMAL_VECTOR(:,:) !< 2D array containing all possible normal vectors at NODE_ID.
    !Local variables
    INTEGER(INTG)  :: coordinate_idx, counter, idx, MY_COMPUTATIONAL_NODE, node_idx, normal_vector_idx, & 
      & normal_vector_idx_1,normal_vector_idx_2, NUMBER_OF_SURROUNDING_NODES, surrouding_node_idx, & 
        & surrouding_node, surrounding_node_idx, vector_idx
    REAL(RP), ALLOCATABLE                                    :: TANGENT_VECTOR(:,:) 
    TYPE(MESH_TYPE), POINTER                                 :: COUPLED_MESH


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
            DO node_idx = 1, SIZE(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES,1)

              IF(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES(node_idx, coordinate_idx)==0) counter= counter+ 1 

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
          ALLOCATE(NORMAL_VECTOR(3,0.5*(NUMBER_OF_SURROUNDING_NODES-1)*(NUMBER_OF_SURROUNDING_NODES)+1), STAT=ERR) ! (n-1)n/2 + 1
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
! Calculate all possible normal angles on node NODE_ID of the interface mesh graph G_I
!==================================================================================================================
  SUBROUTINE INTERFACE_MESH_CALCULATE_NORMAL(COUPLED_DECOMPOSITION, NODE_ID, NORMAL_VECTOR, ERR, ERROR, *) 

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER ,INTENT(INOUT) :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    INTEGER(INTG), INTENT(OUT)                               :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                        :: ERROR !<The error string.
    INTEGER(INTG)                                            :: NODE_ID !<The global node id to find normal at. 
    REAL(RP), ALLOCATABLE                                    :: NORMAL_VECTOR(:,:) !< 2D array containing all possible normal vectors at NODE_ID.
    !Local variables
    INTEGER(INTG)  :: coordinate_idx, counter, idx, MY_COMPUTATIONAL_NODE, node_idx, normal_vector_idx, & 
      & normal_vector_idx_1,normal_vector_idx_2, NUMBER_OF_SURROUNDING_NODES, surrouding_node_idx, & 
        & surrouding_node, surrounding_node_idx, vector_idx
    REAL(RP), ALLOCATABLE                                    :: TANGENT_VECTOR(:,:) 
    TYPE(MESH_TYPE), POINTER                                 :: INTERFACE_MESH


    ENTERS(" INTERFACE_MESH_CALCULATE_NORMAL",ERR,ERROR,*999)
    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH
      IF(ASSOCIATED(INTERFACE_MESH)) THEN

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


        IF(INTERFACE_MESH%NUMBER_OF_DIMENSIONS == 2) THEN 


          DO coordinate_idx = 1, 3
            counter= 0 
            DO node_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES,1)

              IF(COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES(node_idx, coordinate_idx)==0) counter= counter+ 1 

            END DO !node_idx
            IF(counter== SIZE(COUPLED_DECOMPOSITION%INTERFACE_MESH_COORDINATES,1)) EXIT  
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

        ELSE IF(INTERFACE_MESH%NUMBER_OF_DIMENSIONS == 3) THEN 
             
          IF(ALLOCATED(NORMAL_VECTOR)) DEALLOCATE(NORMAL_VECTOR)
          ALLOCATE(NORMAL_VECTOR(3,0.5*(NUMBER_OF_SURROUNDING_NODES-1)*(NUMBER_OF_SURROUNDING_NODES)+1), STAT=ERR) ! (n-1)n/2 + 1
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
! In the following subroutine algorithm B of fixed vertex partitioning is implemented.
!=================================================================================================================
  SUBROUTINE COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_B(COUPLED_DECOMPOSITION, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER ,INTENT(INOUT) :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    INTEGER(INTG), INTENT(OUT)                               :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)                        :: ERROR !<The error string

    !Local Variables
    REAL(RP)                                     :: INTERFACE_MESH_MAXIMUM_EDGE_LENGTH 
    TYPE(MESH_TYPE), POINTER                     :: COUPLED_MESH, INTERFACE_MESH, NEW_COUPLED_MESH
    TYPE(FIELD_TYPE), POINTER                    :: FIELD_COUPLED_MESH, FIELD_INTERFACE_MESH 
    REAL(RP), ALLOCATABLE                        :: NODE_WEIGHTS_TO_BE_PROJECTED(:) 
    TYPE(DECOMPOSITION_TYPE), POINTER            :: COUPLED_MESH_DECOMPOSITION, INTERFACE_MESH_DECOMPOSITION, & 
      & NEW_COUPLED_MESH_DECOMPOSITION
    INTEGER(INTG)                                :: element_idx, MESH_ELEMENTS, MESH_NODES,  MY_COMPUTATIONAL_NODE,  node_idx,  & 
      & NUM, NUMBER_OF_COMPUTATIONAL_NODES, PARMETIS_OPTIONS(3), proc_idx_send, proc_idx_receive, proc_idx, PROC_ID, & 
        & STATUS(MPI_STATUS_SIZE), TOTAL_NODES, vtx_dist_idx   
    INTEGER(INTG), ALLOCATABLE                   :: FLIP_PARTITION(:,:), NEW_TO_OLD_INDEX_MAPPING(:,:), & 
      & NODES_TO_IMPOSE_WEIGHTS_ON(:), TEMP_ARRAY(:)
    INTEGER(INTG)                                :: xadj_counter, adjncy_size, counter, new_node_idx 
    TYPE(REGION_TYPE), POINTER                   :: COUPLED_MESH_REGION, PARENT_REGION
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER        :: COUPLED_MESH_COORDINATE_SYSTEM
    TYPE(NODES_TYPE), POINTER                    :: COUPLED_MESH_NODE
    TYPE(meshElementsType), POINTER              :: COUPLED_MESH_ELEMENTS
   

    ENTERS("COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_B",ERR,ERROR,*999)

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

                  TOTAL_NODES = SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES,1) !maybe remove this.

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
                ! Initialite the coordinate system.
                  NULLIFY(COUPLED_MESH_COORDINATE_SYSTEM)
                  CALL COORDINATE_SYSTEM_CREATE_START(& 
                    & 9999+COUPLED_DECOMPOSITION%mesh_idx,COUPLED_MESH_COORDINATE_SYSTEM,ERR,ERROR,*999)
                  CALL COORDINATE_SYSTEM_DIMENSION_SET(&
                    & COUPLED_MESH_COORDINATE_SYSTEM,COUPLED_MESH%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  CALL COORDINATE_SYSTEM_CREATE_FINISH(COUPLED_MESH_COORDINATE_SYSTEM,ERR,ERROR,*999)

                  ! Initialize the mesh region. 
                  NULLIFY(COUPLED_MESH_REGION)
                  CALL REGION_INITIALISE(COUPLED_MESH_REGION,ERR,ERROR,*999)
                  NULLIFY(COUPLED_MESH_REGION)
                  CALL REGION_CREATE_START(777+COUPLED_DECOMPOSITION%mesh_idx, & 
                    & COUPLED_MESH%REGION%PARENT_REGION,COUPLED_MESH_REGION,ERR,ERROR,*999)
                  CALL REGION_COORDINATE_SYSTEM_SET(COUPLED_MESH_REGION,COUPLED_MESH_COORDINATE_SYSTEM,ERR,ERROR,*999)
                  CALL REGION_CREATE_FINISH(COUPLED_MESH_REGION,ERR,ERROR,*999) 
                  ! Initialize the mesh coupled mesh object from here onward. 
                  CALL MESH_NUMBER_OF_ELEMENTS_GET(COUPLED_MESH,MESH_ELEMENTS,ERR,ERROR,*999)        
                  NULLIFY(COUPLED_MESH_NODE)
                  CALL NODES_CREATE_START(COUPLED_MESH_REGION,9999+COUPLED_DECOMPOSITION%mesh_idx,COUPLED_MESH_NODE,ERR,ERROR,*999)
                  CALL NODES_CREATE_FINISH(COUPLED_MESH_NODE,ERR,ERROR,*999)
                  NULLIFY(NEW_COUPLED_MESH)
                  CALL MESH_INITIALISE(NEW_COUPLED_MESH,ERR,ERROR,*999)
                  NULLIFY(NEW_COUPLED_MESH)
                  CALL MESH_CREATE_START(9999+COUPLED_DECOMPOSITION%mesh_idx,COUPLED_MESH_REGION, & 
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
                    & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE))))
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
 
                    NUM = SQRT(REAL(DOT_PRODUCT(FLIP_PARTITION &
                      & (MINLOC(ABS(FLIP_PARTITION(:,2)-COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_idx-1))),1), & 
         
                    FLIP_PARTITION(MINLOC(ABS(FLIP_PARTITION(:,2)-COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_idx-1))),1))))
 
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

                  !DO node_idx= 0, NUMBER_OF_COMPUTATIONAL_NODES-1
                  !   print *, node_idx, "==", count_integer(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN,node_idx)
                  !END DO
 
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
    EXITS("COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_B.")
    RETURN 
!!999 ERRORSEXITS("COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_B.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_B


!=================================================================================================================
! In the following subroutine algorithm A of fixed vertex partitioning is implemented.
!=================================================================================================================
  SUBROUTINE COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A(COUPLED_DECOMPOSITION, ERR, ERROR, *)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER ,INTENT(INOUT) :: COUPLED_DECOMPOSITION !< Coupled mesh decomposition type object.
    INTEGER(INTG), INTENT(OUT)                            :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT)                     :: ERROR !<The error string

    !Local Variables
    REAL(RP) :: INTERFACE_MESH_MAXIMUM_EDGE_LENGTH 
    TYPE(MESH_TYPE), POINTER :: COUPLED_MESH, INTERFACE_MESH
    TYPE(FIELD_TYPE), POINTER:: FIELD_COUPLED_MESH, FIELD_INTERFACE_MESH 
    INTEGER(INTG), ALLOCATABLE :: INTER_EDGES(:,:)
    REAL(RP), ALLOCATABLE ::  NODE_WEIGHTS_TO_BE_PROJECTED(:) 
    TYPE(DECOMPOSITION_TYPE), POINTER :: COUPLED_MESH_DECOMPOSITION, INTERFACE_MESH_DECOMPOSITION
    INTEGER(INTG):: inter_edge_idx,I, NUMBER_OF_COMPUTATIONAL_NODES, MY_COMPUTATIONAL_NODE, PROC_ID, & 
      & STATUS(MPI_STATUS_SIZE), PARMETIS_OPTIONS(3), proc_idx_send, proc_idx_receive, proc_idx

    ENTERS("COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A",ERR,ERROR,*999)

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

                  MY_COMPUTATIONAL_NODE         = COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  NUMBER_OF_COMPUTATIONAL_NODES = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  ! Gather the coupled mesh coordinates in one data structure.
                  CALL GATHER_MESH_COORDINATES(FIELD_COUPLED_MESH, COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES, &
                    & ERR, Error, *999)
                  ! Calculate vertex adjacencies
                  CALL GET_SURROUNDING_NODES(COUPLED_MESH, ERR, error, *999)
                  ! Calculate maximum edge length of the interfac mesh
                  CALL GET_MAXIMUM_EDGE_LENGTH(COUPLED_DECOMPOSITION, INTERFACE_MESH_MAXIMUM_EDGE_LENGTH, ERR, ERROR, *999)

                  ! IF(ALLOCATED(COUPLED_DECOMPOSITION%INTER_EDGES)) 
                  ! BUild intermediate edges between the coupled mesh graph and the iterface graph
                  CALL COUPLED_DECOMPOSITION_GET_INTER_EDGES(COUPLED_DECOMPOSITION, INTERFACE_MESH_MAXIMUM_EDGE_LENGTH, & 
                    & ERR, ERROR, *999)
                  ! Calculate vertex weights to be projected on the coupled mesh graph
                  CALL COUPLED_DECOMPOSITION_NODE_WEIGHT_SCHEME_A_GET(COUPLED_DECOMPOSITION, &
                    & NODE_WEIGHTS_TO_BE_PROJECTED, ERR, Error, *999)
                  ! Initialize, create start and create finish the decomposition object
                  NULLIFY(COUPLED_MESH%DECOMPOSITIONS)
              
                  CALL DECOMPOSITIONS_INITIALISE(COUPLED_MESH,ERR,ERROR,*999)
  
                  CALL DECOMPOSITION_CREATE_START(COUPLED_DECOMPOSITION%mesh_idx,COUPLED_MESH,COUPLED_MESH_DECOMPOSITION, & 
                    & ERR,ERROR,*999)

                  CALL DECOMPOSITION_NODE_BASED_DECOMPOSITION_SET(COUPLED_MESH_DECOMPOSITION, .TRUE. , ERR, error, *999)

                  CALL DECOMPOSITION_TYPE_SET(COUPLED_MESH_DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)

                  CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(COUPLED_MESH_DECOMPOSITION,NUMBER_OF_COMPUTATIONAL_NODES,err,error,*999)

                  CALL DECOMPOSITION_SET_DEFAULT_VTX_DIST(COUPLED_MESH_DECOMPOSITION,ERR,ERROR,*999)

                  CALL DECOMPOSITION_NUMBER_OF_CONSTRAINTS_SET(COUPLED_MESH_DECOMPOSITION, 1_INTG, err, error, *999)

                  CALL GET_ADJNCY_AND_XADJ(COUPLED_MESH_DECOMPOSITION,ERR,error,*999)

                  COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION=>&
                    & COUPLED_MESH_DECOMPOSITION

                  ! Prescribe the vertex weight on the coupled mesh graph.
                  CALL COUPLED_DECOMPOSITION_SET_NODE_WEIGHT(COUPLED_DECOMPOSITION, NINT(NODE_WEIGHTS_TO_BE_PROJECTED), & 
                    & ERR, Error, *999)
                  ! Prescribe the edge weight on the coupled mesh graph.
                  CALL COUPLED_DECOMPOSITION_SET_EDGE_WEIGHT(COUPLED_DECOMPOSITION, ERR, Error, *999)

                  COUPLED_MESH_DECOMPOSITION=> &
                    & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION

                  CALL DECOMPOSITION_CREATE_FINISH(COUPLED_MESH_DECOMPOSITION, ERR, Error, *999)
                  ! Flip parition if required. 
                  CALL COUPLED_DECOMPOSITION_FLIP_PARTITION(COUPLED_DECOMPOSITION,ERR,Error,*999)

                  DEALLOCATE(COUPLED_DECOMPOSITION%COUPLED_MESH_COORDINATES)
                  DEALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES)
                  DEALLOCATE(NODE_WEIGHTS_TO_BE_PROJECTED) 

                  COUPLED_DECOMPOSITION%mesh_idx = COUPLED_DECOMPOSITION%mesh_idx + 1 

                ELSE
                  CALL FlagError(" Interface mesh decomposition  is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError(" Coupled mesh decomposition  is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError(" Interface mesh field  is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError(" Coupled mesh decomposition  is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError(" Coupled mesh field  is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError(" Coupled mesh  is not associated.",ERR,ERROR,*999)
      ENDIF                
    ELSE
     CALL FlagError(" Coupled decomposition  is not associated.",ERR,ERROR,*999)
    ENDIF
   EXITS("COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A")
   RETURN
!999 ERRORSEXITS("COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A",err,error)
999 RETURN 1


  END SUBROUTINE COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A

!========================================================================================
 ! The following subroutine trivially decomposes the interface graph G_I.
!========================================================================================


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

        CALL DECOMPOSITION_CREATE_START(INT(CMISS_RANDOM_SEEDS(1),INTG),INTERFACE_MESH,INTERFACE_DECOMPOSITION,ERR,ERROR,*999)

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
!999 ERRORSEXITS("DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMPOSITION_SET",err,error)
999 RETURN 1

  END SUBROUTINE DECOMPOSITION_INTERFACE_MESH_TRIVIAL_DECOMPOSITION_SET
!===================================================================================================================!
! The following subroutine adds the geometric field information of the interface graph G_{i}  in COUPLED_DECOMPOSITION object.
!===================================================================================================================!
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
! The following subroutine adds the geometric field information of the coupled mesh graph G_{i} in COUPLED_DECOMPOSITION object.
!===================================================================================================================!

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
! The following subroutine initialize members of COUPLED_DECOMPOSITION object.
!===================================================================================================================!

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
 ! IN the following subroutine the "coupling aware" mesh partitioning is implemented.
!====================================================================================!

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

                        IF(.FALSE.) THEN ! Scheme A starts here for fixed partitioning.
                          CALL COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A(COUPLED_DECOMPOSITION,ERR,ERROR,*999)
                          CALL COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_A(COUPLED_DECOMPOSITION,ERR,ERROR,*999)
                        ELSE   ! Scheme B starts here for fixed partitioning.
                          CALL COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_B(COUPLED_DECOMPOSITION,ERR,ERROR,*999)  
                          CALL COUPLED_DECOMPOSITION_FIXED_PARTITIONING_SCHEME_B(COUPLED_DECOMPOSITION,ERR,ERROR,*999) 
                        END IF 

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
!======================================================================================================================
! The following subroutine generates the set of weights supposed to be imposed on the vertices of the sub-graph G_i^I.
!======================================================================================================================
  SUBROUTINE COUPLED_DECOMPOSITION_NODE_WEIGHT_SCHEME_A_GET(COUPLED_DECOMPOSITION, &
   NODE_WEIGHTS_TO_BE_PROJECTED, ERR, Error, *)
  !Argument variables
   TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)  :: COUPLED_DECOMPOSITION !<Coupled Decomposition type object to create finish.
   REAL(RP), INTENT(OUT), ALLOCATABLE                        :: NODE_WEIGHTS_TO_BE_PROJECTED(:) !<NODE_WEIGHTS_TO_BE_PROJECTED(node_idx) contains vertex weight to be imposed on a vertex of subgraph G_i^I.
   INTEGER(INTG), INTENT(OUT)                                :: ERR !<The error code.
   TYPE(VARYING_STRING), INTENT(OUT)                         :: ERROR !<The error string.

  !Local variables
   INTEGER(INTG)                                             :: domain_idx,idx,node_idx
   TYPE(MESH_TYPE), POINTER                                  :: COUPLED_MESH, INTERFACE_MESH
   TYPE(DECOMPOSITION_TYPE), POINTER                         :: INTERFACE_DECOMPOSITION
   INTEGER(INTG), ALLOCATABLE                                :: NUMBER_OF_NODES_PER_PROCESSOR(:), &
                                                                & PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION(:)
   REAL(RP)                                                  :: NODAL_WEIGHTS_TO_BE_DISTRIBUTED

   ENTERS("COUPLED_DECOMPOSITION_NODE_WEIGHT_SCHEME_A_GET",ERR,ERROR,*999)

   IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
     
     INTERFACE_DECOMPOSITION=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION

     IF(ASSOCIATED(INTERFACE_DECOMPOSITION)) THEN

     INTERFACE_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION%MESH

       IF(ASSOCIATED(INTERFACE_MESH)) THEN

         COUPLED_MESH=>&
           & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%MESH
         IF(ASSOCIATED(COUPLED_MESH)) THEN
       
           NODAL_WEIGHTS_TO_BE_DISTRIBUTED = 1.*COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfnodes

           ALLOCATE(NODE_WEIGHTS_TO_BE_PROJECTED(SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)), STAT=ERR)
           IF(ERR /=0) CALL FlagError(" Cannot allocate NODE_WEIGHTS_TO_BE_PROJECTED .",ERR,ERROR,*999) 

    
           ALLOCATE(PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION(SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)), STAT=ERR)
           IF(ERR /=0) CALL FlagError(" Cannot allocate PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION.",ERR,ERROR,*999) 

           ALLOCATE(NUMBER_OF_NODES_PER_PROCESSOR(INTERFACE_DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
           IF(ERR /=0) CALL FlagError(" Cannot allocate NUMBER_OF_NODES_PER_PROCESSOR.",ERR,ERROR,*999) 


           DO idx = 1, SIZE(PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION,1)

           PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION(idx) = &
             & INTERFACE_DECOMPOSITION%NODE_DOMAIN(COUPLED_DECOMPOSITION%INTER_EDGES(idx,1)-1)

           END DO ! Idx


           NUMBER_OF_NODES_PER_PROCESSOR = 0 

    
           DO node_idx = 1 , SIZE(PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION,1) 

             domain_idx = PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION(node_idx)
             NUMBER_OF_NODES_PER_PROCESSOR(domain_idx+1) = NUMBER_OF_NODES_PER_PROCESSOR(domain_idx+1)  +1 

           END DO


           DO idx = 1, SIZE(PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION,1)

             node_idx   = COUPLED_DECOMPOSITION%INTER_EDGES(idx,1)
             domain_idx = INTERFACE_DECOMPOSITION%NODE_DOMAIN(node_idx-1)
             NODE_WEIGHTS_TO_BE_PROJECTED(idx) = NODAL_WEIGHTS_TO_BE_DISTRIBUTED/NUMBER_OF_NODES_PER_PROCESSOR(domain_idx+1)

           END DO 

           IF(ALLOCATED(NUMBER_OF_NODES_PER_PROCESSOR)) DEALLOCATE(NUMBER_OF_NODES_PER_PROCESSOR)
           IF(ALLOCATED(PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION)) & 
             & DEALLOCATE(PROCESSORS_TO_BE_IMPOSED_ON_COUPLED_PARTITION)

         ELSE
           CALL FlagError(" Coupled mesh object is not associated.",ERR,ERROR,*999)
         ENDIF
       ELSE
         CALL FlagError(" Interface mesh is not associated.",ERR,ERROR,*999)
       ENDIF
     ELSE
       CALL FlagError(" Interface decomposition is not associated.",ERR,ERROR,*999)
     ENDIF
   ELSE
     CALL FlagError(" Coupled decomposition is not associated.",ERR,ERROR,*999)
   ENDIF
   EXITS("COUPLED_DECOMPOSITION_NODE_WEIGHT_SCHEME_A_GET.")
   RETURN
!999 ERRORSEXITS("COUPLED_DECOMPOSITION_NODE_WEIGHT_SCHEME_A_GET.",ERR,ERROR)
999 RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_NODE_WEIGHT_SCHEME_A_GET
!=========================================================================================================================

  SUBROUTINE GATHER_MESH_COORDINATES(GEOMETRIC_FIELD, MESH_COORDINATES,  ERR, Error, * )

   !Argument variables
   TYPE(FIELD_TYPE), POINTER, INTENT(INOUT)             :: GEOMETRIC_FIELD ! A pointer to the geometric field of the mesh whose coordinates are to be gathered.
   REAL(RP), INTENT(OUT), ALLOCATABLE                   :: MESH_COORDINATES(:,:) ! Data structure where the geometric field coordinates are gathered.
   INTEGER(INTG), INTENT(OUT)                           :: ERR !<The error code
   TYPE(VARYING_STRING), INTENT(OUT)                    :: ERROR !<The error string
    
   !Local variables
   INTEGER(INTG)                                        :: component_idx, DOMAIN_NUMBER, node_idx, MY_COMPUTATIONAL_NODE
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
      
           ALLOCATE(MESH_COORDINATES(MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes,3), STAT=ERR)
           IF(ERR /= 0)   CALL FlagError(" Cannot allocate MESH_COORDINATES.",ERR,ERROR,*999)

           MESH_COORDINATES=0
      
           DO node_idx = 1, MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes

             CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,node_idx,1_INTG, &
               & DOMAIN_NUMBER,ERR,ERROR,*999)
              
             IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE)  THEN                

               DO component_idx = 1, 3

                 CALL FIELD_PARAMETER_SET_GET_NODE(&
                   & GEOMETRIC_FIELD,FIELD_VARIABLE_MESH%VARIABLE_TYPE,1_INTG, 1_INTG, 1_INTG, &
                     & node_idx,component_idx,MESH_COORDINATES(node_idx,component_idx),ERR,ERROR,*997)  

                 IF(.FALSE.) THEN
997                MESH_COORDINATES(node_idx,component_idx) = 0._RP    
                 END IF

               END DO ! component_idx
  
             END IF ! IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE)        
    
           END DO ! node_idx

           ALLOCATE(MESH_COORDINATES_NEW(MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes*3), STAT=ERR)
           IF(ERR /= 0)   CALL FlagError(" Cannot allocate MESH_COORDINATES_NEW.",ERR,ERROR,*999)
           MESH_COORDINATES_NEW = 0.  

 
           CALL MPI_ALLREDUCE(MESH_COORDINATES, MESH_COORDINATES_NEW, MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes*3, &
             & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,Err)

           MESH_COORDINATES = 0.
           MESH_COORDINATES = RESHAPE(MESH_COORDINATES_NEW , (/MESH%TOPOLOGY(1)%PTR%NODES%NUmberOfNOdes,3/))

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
!999 ERRORSEXITS("GATHER_MESH_COORDINATES.",ERR,ERROR)
999 RETURN 1

  END SUBROUTINE GATHER_MESH_COORDINATES
!=========================================================================================================================
! The following subroutines impose node weights NODE_WEIGHTS_TO_BE_SET on the coupled mesh. 
!=========================================================================================================================

  SUBROUTINE COUPLED_DECOMPOSITION_SET_NODE_WEIGHT(COUPLED_DECOMPOSITION,NODE_WEIGHTS_TO_BE_SET,ERR,Error,*)
  !Argument variables
   TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)          :: COUPLED_DECOMPOSITION !<Derive type object that stores information about coupled decomposition. 
   INTEGER(INTG), INTENT(IN)                                         :: NODE_WEIGHTS_TO_BE_SET(:) !NOde weights to be set on ...
                                                                        !... the interface region of the coupled mesh graph G_i^I.
   INTEGER(INTG), INTENT(OUT)                                        :: ERR !<The error code
   TYPE(VARYING_STRING), INTENT(OUT)                                 :: ERROR !<The error string

  !Local variables
   INTEGER(INTG)                                                     :: inter_edge_idx, MY_COMPUTATIONAL_NODE
   TYPE(DECOMPOSITION_TYPE), POINTER                                 :: COUPLED_MESH_DECOMPOSITION
 
   ENTERS("COUPLED_DECOMPOSITION_SET_NODE_WEIGHT.",ERR,ERROR,*999)

   IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
     
     COUPLED_MESH_DECOMPOSITION=>&
       & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION

     IF(ASSOCIATED(COUPLED_MESH_DECOMPOSITION)) THEN

       MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)

       ALLOCATE(COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET( & 
         & 1:COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS* &
           & (COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1) &  
             & - COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE))), STAT=ERR)

       IF(ERR/=0) &
         & CALL FlagError("Cannot allocate COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET data structure.", ERR,ERROR,*999) 

       COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET = 1 
    
       DO inter_edge_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)

         IF(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2) <= &
           & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1) .AND. &
             & COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2) > &
               & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE)) THEN    


           COUPLED_MESH_DECOMPOSITION%NODE_WEIGHT_SET( &
             & COUPLED_MESH_DECOMPOSITION%NUMBER_OF_CONSTRAINTS*(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2)- &
               & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE))) = NODE_WEIGHTS_TO_BE_SET(inter_edge_idx)

         ENDIF
       END DO ! inter_edge_idx
     ELSE
       CALL FlagError(" COUPLED_MESH_DECOMPOSITION.",ERR,ERROR,*999)
     ENDIF
   ELSE
     CALL FlagError(" COUPLED DECOMPOSITION is not associated.",ERR,ERROR,*999)
   ENDIF
   EXITS("COUPLED_DECOMPOSITION_SET_NODE_WEIGHT.")
   RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_SET_NODE_WEIGHT.",ERR,ERROR)
   RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_SET_NODE_WEIGHT


!=========================================================================================================================
! The following subroutine impose weights on the inter edges. THe values of the weights is based on fact whether the domain ...
! .. assigned to ends of the inter edges is similar or not. If similar then imposed a huge edge weight and if not similar ....
! .. then imposed an edge weight of 0 on the inter edge.

  SUBROUTINE COUPLED_DECOMPOSITION_SET_EDGE_WEIGHT(COUPLED_DECOMPOSITION,ERR,Error,*)
    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)          :: COUPLED_DECOMPOSITION !<Object that stores information about coupled decomposition. 
    INTEGER(INTG), INTENT(OUT)                                        :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                                 :: ERROR !<The error string.

    !Local variables
    INTEGER(INTG) :: adjacency_idx, CORRECT_INDEX, inter_edge_idx, INTER_EDGE_COUPLED_MESH_NODE_CORRECTED, &
      & INTER_EDGE_COUPLED_MESH_NODE, MY_COMPUTATIONAL_NODE, PROCESSOR_I, PROCESSOR_J, VTX_DIST, &
        & x_adj_idx, x_adj_idx_final, x_adj_idx_initial
    TYPE(DECOMPOSITION_TYPE), POINTER  :: COUPLED_MESH_DECOMPOSITION, INTERFACE_MESH_DECOMPOSITION
    TYPE(MESH_TYPE), POINTER           :: COUPLED_MESH, INTERFACE_MESH
    INTEGER(INTG)  , PARAMETER         :: LARGE_EDGE_WEIGHT = 1E5_INTG
    INTEGER(INTG)  , PARAMETER         :: SMALL_EDGE_WEIGHT = 0_INTG

 
    ENTERS("COUPLED_DECOMPOSITION_SET_EDGE_WEIGHT",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
     
      COUPLED_MESH_DECOMPOSITION=> &
        & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION

      IF(ASSOCIATED(COUPLED_MESH_DECOMPOSITION)) THEN
    
        INTERFACE_MESH_DECOMPOSITION=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION
  
        IF(ASSOCIATED(INTERFACE_MESH_DECOMPOSITION)) THEN 
          
          COUPLED_MESH=>COUPLED_MESH_DECOMPOSITION%MESH

          IF(ASSOCIATED(COUPLED_MESH)) THEN
         
            INTERFACE_MESH=>INTERFACE_MESH_DECOMPOSITION%MESH 
  
            IF(ASSOCIATED(INTERFACE_MESH)) THEN   
    
              MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)

              ALLOCATE(COUPLED_MESH_DECOMPOSITION%ADJWT( & 
                & SIZE(COUPLED_MESH_DECOMPOSITION%ADJNCY,1)), STAT=ERR)
              IF(ERR/=0) &
                & CALL FlagError("COUPLED_MESH_DECOMPOSITION%ADJWT data structure is already allocated.", & 
                  & ERR,ERROR,*999)
 
              COUPLED_MESH_DECOMPOSITION%ADJWT = 1

              DO inter_edge_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)


                IF(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2) <= &
                  & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE+1) .AND. &
                    & COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2) > &
                      & COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE)) THEN 

                  PROCESSOR_I = & ! Domain assigned to vertex i in inter_edge_idx (i,j)
                    & INTERFACE_MESH_DECOMPOSITION%NODE_DOMAIN(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,1)-1)    

                  INTER_EDGE_COUPLED_MESH_NODE = COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2)

                  DO adjacency_idx=1, &
                    & SIZE(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%NODES(INTER_EDGE_COUPLED_MESH_NODE)%surroundingNOdes)

                    IF(ANY(COUPLED_MESH%TOPOLOGY(1)%PTR%&  !If any of the surrounding node is connected with the inter edges.
                      & NODES%NODES(INTER_EDGE_COUPLED_MESH_NODE)%surroundingNOdes(adjacency_idx)==&
                        & COUPLED_DECOMPOSITION%INTER_EDGES(:,2))) THEN

                      VTX_DIST =  COUPLED_MESH_DECOMPOSITION%VTX_DIST(MY_COMPUTATIONAL_NODE)

                      INTER_EDGE_COUPLED_MESH_NODE_CORRECTED = COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2)-1 ! to make the node nubering starting from 0.

                      ! x_adj_idx_initial and x_adj_idx_final defines the range of nodes to find the node j linked to ... 
                      ! ... inter edge inter_edge_idx (i,j).
                      x_adj_idx_initial = COUPLED_MESH_DECOMPOSITION%XADJ(INTER_EDGE_COUPLED_MESH_NODE_CORRECTED-VTX_DIST)

                      x_adj_idx_final = COUPLED_MESH_DECOMPOSITION%XADJ(INTER_EDGE_COUPLED_MESH_NODE_CORRECTED-VTX_DIST+1)-1

                      DO x_adj_idx =  x_adj_idx_initial, x_adj_idx_final 

                        IF(ANY(COUPLED_MESH_DECOMPOSITION%ADJNCY(x_adj_idx) == &  
                          & COUPLED_DECOMPOSITION%INTER_EDGES(:,2)-1)) THEN

                          CORRECT_INDEX  = MINLOC(ABS(COUPLED_DECOMPOSITION%INTER_EDGES(:,2)-1- &
                            & COUPLED_MESH_DECOMPOSITION%ADJNCY(x_adj_idx)),1)

                          ! Domain assigned to vertex j in inter_edge_idx (i,j).

                          PROCESSOR_J = INTERFACE_MESH_DECOMPOSITION%NODE_DOMAIN( &
                            & COUPLED_DECOMPOSITION%INTER_EDGES(CORRECT_INDEX,1)-1)

                          IF(PROCESSOR_J == PROCESSOR_I) THEN ! If domains assigned to vertex  i and j are the same....
                                                              ! ... then impose a large edge weight         
                            COUPLED_MESH_DECOMPOSITION%ADJWT(x_adj_idx+1) = LARGE_EDGE_WEIGHT

                          ELSE IF(PROCESSOR_J /= PROCESSOR_I) THEN ! If domains assigned to vertex  i and j are the same....
                                                              ! ... then impose a small edge weight 

                            COUPLED_MESH_DECOMPOSITION%ADJWT(x_adj_idx+1) = SMALL_EDGE_WEIGHT

                          END IF !(PROCESSOR_J == PROCESSOR_I)
                
                       END IF ! (ANY(COUPLED_MESH_DECOMPOSITION%ADJNCY(x_adj_idx) == &

                     END DO ! x_adj_idx

                   END IF ! (ANY(COUPLED_MESH%TOPOLOGY(1)%PTR%&

                 END DO !adjacency_idx

               END IF  

             END DO ! inter_edge_idx
           ELSE
             CALL FlagError("COUPLED_MESH is not associated.",ERR,ERROR,*999)
           END IF 
         ELSE
           CALL FlagError("INTERFACE_MESH is not associated.",ERR,ERROR,*999)
         END IF 
        ELSE
          CALL FlagError("INTERFACE_MESH_DECOMPOSITION is not associated.",ERR,ERROR,*999)
        END IF 
      ELSE
        CALL FlagError(" COUPLED_MESH_DECOMPOSITION is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError(" COUPLED_DECOMPOSITION is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_DECOMPOSITION_SET_EDGE_WEIGHT")
    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_SET_EDGE_WEIGHT",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_SET_EDGE_WEIGHT


! ====================================================================================================================
! THe following subroutine flips the processor ids (or sub-domain ids) based on the subdomain ids of the interface. 
! THis routine makes sure that processor ids on both side of the inter edges are identical. 

  SUBROUTINE COUPLED_DECOMPOSITION_FLIP_PARTITION(COUPLED_DECOMPOSITION,ERR,Error,*)
    !Argument variables.
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(INOUT) :: COUPLED_DECOMPOSITION !<Object that stores information about coupled decomposition. 
    INTEGER(INTG), INTENT(OUT)                               :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT)                        :: ERROR !<The error string.

    ! Local variables.

    TYPE(DECOMPOSITION_TYPE), POINTER         :: COUPLED_MESH_DECOMPOSITION, INTERFACE_MESH_DECOMPOSITION

    TYPE(LIST_PTR_TYPE), ALLOCATABLE          :: FLIP_PARTITION_COUPLED_MESH_LIST(:), &
      & FLIP_PARTITION_INTERFACE_MESH_LIST(:)

    INTEGER(INTG)               :: domain_idx, domain_idx_2, DOMAIN_ID, idx, inter_edge_idx,NUMBER_OF_PROCESSORS, & 
                                   & node_domain_idx, proc_idx, LOC

    INTEGER(INTG),ALLOCATABLE   :: FLIP_PARTITION_COUPLED_MESH_PROCESSOR(:), FLIP_PARTITION(:,:) , &
      & NUMBER_OF_PROCESSOR_IDS_IN_LIST(:)   

    ENTERS("COUPLED_DECOMPOSITION_FLIP_PARTITION",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN

      COUPLED_MESH_DECOMPOSITION=> &
        & COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION

      IF(ASSOCIATED(COUPLED_MESH_DECOMPOSITION)) THEN

        INTERFACE_MESH_DECOMPOSITION=> &
        & COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION

        IF(ASSOCIATED(INTERFACE_MESH_DECOMPOSITION)) THEN       
        
          ALLOCATE(FLIP_PARTITION(COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS,2), STAT=ERR)
          IF(ERR/=0) CALL FlagError("FLIP_PARTITION is already allocated.", & 
           & ERR,ERROR,*999)     

          ALLOCATE(FLIP_PARTITION_INTERFACE_MESH_LIST(COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate FLIP_PARTITION_INTERFACE_MESH_LIST array.",ERR,ERROR,*999)

          ALLOCATE(FLIP_PARTITION_COUPLED_MESH_LIST(COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate FLIP_PARTITION_COUPLED_MESH_LIST array.",ERR,ERROR,*999)

          ALLOCATE(NUMBER_OF_PROCESSOR_IDS_IN_LIST(COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS), STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_PROCESSOR_IDS_IN_LIST array.",ERR,ERROR,*999)

          DO domain_idx = 1, COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS

            NULLIFY(FLIP_PARTITION_INTERFACE_MESH_LIST(domain_idx)%PTR)

            CALL LIST_CREATE_START(FLIP_PARTITION_INTERFACE_MESH_LIST(domain_idx)%PTR, &
              & ERR, ERROR, *999)

            CALL LIST_DATA_TYPE_SET(FLIP_PARTITION_INTERFACE_MESH_LIST(domain_idx)%PTR, &
              & LIST_INTG_TYPE, ERR, ERROR, *999)

            CALL LIST_INITIAL_SIZE_SET(FLIP_PARTITION_INTERFACE_MESH_LIST(domain_idx)%PTR, &
              & SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1), ERR, ERROR, *999)

            CALL LIST_CREATE_FINISH(FLIP_PARTITION_INTERFACE_MESH_LIST(domain_idx)%PTR, &
              & ERR, ERROR, *999)


            NULLIFY(FLIP_PARTITION_COUPLED_MESH_LIST(domain_idx)%PTR)

            CALL LIST_CREATE_START(FLIP_PARTITION_COUPLED_MESH_LIST(domain_idx)%PTR, &
              & ERR, ERROR, *999)

            CALL LIST_DATA_TYPE_SET(FLIP_PARTITION_COUPLED_MESH_LIST(domain_idx)%PTR, &
              & LIST_INTG_TYPE, ERR, ERROR, *999)

            CALL LIST_INITIAL_SIZE_SET(FLIP_PARTITION_COUPLED_MESH_LIST(domain_idx)%PTR, &
              & SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1), ERR, ERROR, *999)

            CALL LIST_CREATE_FINISH(FLIP_PARTITION_COUPLED_MESH_LIST(domain_idx)%PTR, &
              & ERR, ERROR, *999)

          END DO ! domain_idx


          DO proc_idx = 1, COUPLED_MESH_DECOMPOSITION%NUMBER_OF_DOMAINS
          
            DO inter_edge_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)

              DOMAIN_ID = INTERFACE_MESH_DECOMPOSITION%NODE_DOMAIN&
                & (COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,1)-1)

              IF(DOMAIN_ID==proc_idx-1) THEN

                CALL LIST_ITEM_ADD(FLIP_PARTITION_COUPLED_MESH_LIST(proc_idx)%PTR, &
                  & COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2)-1), &
                    &  ERR,ERROR,*999)

                CALL LIST_ITEM_ADD(FLIP_PARTITION_INTERFACE_MESH_LIST(proc_idx)%PTR, &
                  & INTERFACE_MESH_DECOMPOSITION%NODE_DOMAIN(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,1)-1), &
                    &  ERR,ERROR,*999)

              END IF 

            END DO !inter_edge_idx

            CALL LIST_DETACH_AND_DESTROY(FLIP_PARTITION_COUPLED_MESH_LIST(proc_idx)%PTR, &
            & NUMBER_OF_PROCESSORS, FLIP_PARTITION_COUPLED_MESH_PROCESSOR, ERR,ERROR,*999)  

            NUMBER_OF_PROCESSOR_IDS_IN_LIST = 0       
 
            DO idx = 1,  NUMBER_OF_PROCESSORS
       
              DOMAIN_ID = FLIP_PARTITION_COUPLED_MESH_PROCESSOR(idx)
 
              NUMBER_OF_PROCESSOR_IDS_IN_LIST(DOMAIN_ID+1) =  NUMBER_OF_PROCESSOR_IDS_IN_LIST(DOMAIN_ID+1) +1     

            END DO ! idx

            FLIP_PARTITION(proc_idx,1) = MAXLOC(NUMBER_OF_PROCESSOR_IDS_IN_LIST,1)-1
            FLIP_PARTITION(proc_idx,2) = proc_idx - 1 

            DEALLOCATE(FLIP_PARTITION_COUPLED_MESH_PROCESSOR)

          END DO  !proc_idx

          DO node_domain_idx = 0, SIZE(COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN,1) -1

            LOC = MINLOC(ABS(FLIP_PARTITION(:,1)-COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_domain_idx)), 1)

            COUPLED_MESH_DECOMPOSITION%NODE_DOMAIN(node_domain_idx) =  FLIP_PARTITION(LOC,2)

          END DO !node_domain_idx

    
          IF(ALLOCATED(FLIP_PARTITION_COUPLED_MESH_PROCESSOR)) DEALLOCATE(FLIP_PARTITION_COUPLED_MESH_PROCESSOR)
          IF(ALLOCATED(FLIP_PARTITION)) DEALLOCATE(FLIP_PARTITION)
          IF(ALLOCATED(NUMBER_OF_PROCESSOR_IDS_IN_LIST)) DEALLOCATE(NUMBER_OF_PROCESSOR_IDS_IN_LIST)

        ELSE
          CALL FlagError(" INTERFACE_MESH_DECOMPOSITION is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError(" COUPLED_MESH_DECOMPOSITION is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError(" COUPLED_DECOMPOSITION is not associated.",ERR,ERROR,*999) 
    END IF
    EXITS("COUPLED_DECOMPOSITION_FLIP_PARTITION")
    RETURN
999 ERRORSEXITS("COUPLED_DECOMPOSITION_FLIP_PARTITION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_FLIP_PARTITION

!========================================================================================================
! The following subroutine updates the new geometric field.
!========================================================================================================
  SUBROUTINE DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD(FIELD,DECOMPOSITION,VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER, INTENT(IN)         :: FIELD !<A pointer to the field to update the geometric parameters for
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: DECOMPOSITION !<The mesh which is generated by the generated mesh \todo is this necessary???
    INTEGER(INTG), INTENT(IN)  :: VARIABLE_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING)               :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER           :: MESH
    INTEGER(INTG)                      :: DOMAIN,MY_COMPUTATIONAL_NODE, node_idx, TOTAL_NODES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING)               :: localError

    ENTERS("DECOMPOSITION_ASSIGN_DECOMPOSITION_FIELD",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN     
     
            MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
            TOTAL_NODES = SIZE(MESH%TOPOLOGY(1)%PTR%NODES%NODES)

            DO node_idx = 1, TOTAL_NODES

              CALL DECOMPOSITION_NODE_DOMAIN_GET( DECOMPOSITION, node_idx, 1_INTG, DOMAIN, ERR, ERROR, *999)

              IF(DOMAIN==MY_COMPUTATIONAL_NODE) THEN      

                CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1_INTG,1_INTG, &
                  & node_idx,1, DECOMPOSITION%NODE_DOMAIN(node_idx-1), ERR,Error,*999)

              END IF
    

            END DO !node_idx

          ELSE
            localError="The specified variable type of "//TRIM(NumberToVString(VARIABLE_TYPE,"*",err,error))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
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
! The following subroutine updates the decomposition field of the coupled mesh graph G_i.
!========================================================================================================

  SUBROUTINE COUPLED_DECOMPOSITION_UPDATE_DECOMPOSITION(COUPLED_DECOMPOSITION,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)         :: DECOMPOSITION !< A decomposition type object to be updated.
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !<Coupled decomposition type object that stores information about the coupled mesh partitioning.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING)            :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER        :: MESH
    INTEGER(INTG)                   :: DOMAIN, mapping_idx, NUMBER_OF_COMPUTATIONAL_NODES, new_node_idx, node_idx, & 
      & TOTAL_NODES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
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

          CALL DECOMPOSITION_CREATE_START(COUPLED_DECOMPOSITION%mesh_idx,MESH,NEW_DECOMPOSITION,ERR,ERROR,*999)

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
! The following subroutine updates the decomposition field of the interface mesh graph G_I.
!========================================================================================================

  SUBROUTINE COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION(COUPLED_DECOMPOSITION,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(INOUT)         :: DECOMPOSITION !<A pointer to the field to update the geometric parameters for
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !<The mesh which is generated by the generated mesh \todo is this necessary???
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING)               :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER           :: MESH
    INTEGER(INTG)                      :: DOMAIN,MY_COMPUTATIONAL_NODE, node_idx, TOTAL_NODES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING)               :: localError

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
!ERRORSEXITS("COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION",ERR,ERROR)
999     RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION


!=====================================================================================================================
!====================================================================================================
 SUBROUTINE COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING(COUPLED_DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !< Coupled decomposition type objects.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING)               :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER           :: MESH
    INTEGER(INTG)                      :: inter_edge_idx, MY_COMPUTATIONAL_NODES, proc_idx, PROC_ID, TOTAL_COMPUTATIONAL_NODES 
    INTEGER(INTG), ALLOCATABLE         :: SUB_DOMAIN_COUNTER(:), TEMP_ARRAY(:,:)
    TYPE(DECOMPOSITION_TYPE), POINTER  :: INTERFACE_MESH_DECOMPOSITION

    ENTERS("COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING",ERR,ERROR,*999)

    IF(ASSOCIATED(COUPLED_DECOMPOSITION)) THEN
      
      INTERFACE_MESH_DECOMPOSITION=> &
        & COUPLED_DECOMPOSITION%COUPLED_FIELDS(3)%PTR%DECOMPOSITION      
    
      IF(ASSOCIATED(INTERFACE_MESH_DECOMPOSITION)) THEN

        TOTAL_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR) ! Get rank.
        MY_COMPUTATIONAL_NODES=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) ! Get id of the processor

        ALLOCATE(SUB_DOMAIN_COUNTER(TOTAL_COMPUTATIONAL_NODES), STAT=ERR)
        IF(ERR/=0)   CALL FlagError("SUB_DOMAIN_countercannot be allocated.",ERR,ERROR,*999)
    
        ALLOCATE(TEMP_ARRAY(SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1),TOTAL_COMPUTATIONAL_NODES), STAT=ERR) !TEMP_ARRAY(:,PROC_idx) contains all the vertices that are assigned sub-domain PROC_idx
        IF(ERR/=0)   CALL FlagError("TEMP_ARRAY cannot be allocated.",ERR,ERROR,*999)
  
        TEMP_ARRAY =  0
        SUB_DOMAIN_counter= 1 

        DO proc_idx = 1, TOTAL_COMPUTATIONAL_NODES

          DO inter_edge_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)
       
            PROC_ID = INTERFACE_MESH_DECOMPOSITION%NODE_DOMAIN(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,1)-1)

            IF(PROC_ID == proc_idx-1) THEN 
             
              TEMP_ARRAY(SUB_DOMAIN_COUNTER(proc_idx),proc_idx)=COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,2)

              SUB_DOMAIN_COUNTER(proc_idx) = SUB_DOMAIN_COUNTER(proc_idx) + 1 
           
            END IF 

          END DO !inter_edge_idx

        END DO !proc_idx
     
        DEALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES) ! Store the vertices in COUPLED_DECOMPOSITION%INTER_EDGES data structure ... 
                                                      !... such that !COUPLED_DECOMPOSITION%INTER_EDGES(:,PROC_idx) contains all the vertices that are assigned sub-domain PROC_idx
        ALLOCATE(COUPLED_DECOMPOSITION%INTER_EDGES(SIZE(TEMP_ARRAY,1),TOTAL_COMPUTATIONAL_NODES)) 
        COUPLED_DECOMPOSITION%INTER_EDGES = TEMP_ARRAY
        DEALLOCATE(TEMP_ARRAY)
        DEALLOCATE(SUB_DOMAIN_COUNTER)
      ELSE
        CALL FlagError("Interface mesh decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Coupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING")
    RETURN
!ERRORSEXITS("COUPLED_DECOMPOSITION_UPDATE_INTERFACE_DECOMPOSITION",ERR,ERROR)
999 RETURN 1
  END SUBROUTINE COUPLED_MESH_VERTICES_TO_IMPOSE_FIXED_PARTITIONING

!=====================================================================================================================
 ! This subroutine builds mapping between numbering of a vertices of the new graph G_{i,merged} and the orginal graph G_i. 
!=====================================================================================================================
 SUBROUTINE COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING(COUPLED_DECOMPOSITION, NEW_TO_OLD_INDEX_MAPPING, ERR,ERROR,*)

    !Argument variables
    TYPE(COUPLED_DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: COUPLED_DECOMPOSITION !<Coupled decomposition type objects.
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: NEW_TO_OLD_INDEX_MAPPING(:,:) !<NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,:). Cotains one-to-one mapping between new vertex ids and old vertex ids  
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
 
        NUMBER_OF_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)

        FLAG1=.FALSE.
        FLAG2=.FALSE.
    
        COUPLED_MESH=>COUPLED_DECOMPOSITION%COUPLED_FIELDS(COUPLED_DECOMPOSITION%mesh_idx)%PTR%DECOMPOSITION%MESH

        ALLOCATE(NEW_TO_OLD_INDEX_MAPPING(COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNOdes,2), STAT=ERR)
        new_vertex_idx = 1 

        DO old_vertex_idx = 1, COUPLED_MESH%TOPOLOGY(1)%PTR%NODES%numberOfNOdes

          FLAG1= .FALSE.
          FLAG2= .FALSE.
          NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,1) = old_vertex_idx 
 
          DO proc_idx = 1, NUMBER_OF_COMPUTATIONAL_NODES

            DO inter_edge_idx = 1, SIZE(COUPLED_DECOMPOSITION%INTER_EDGES,1)

              IF(COUPLED_DECOMPOSITION%INTER_EDGES(inter_edge_idx,proc_idx)==old_vertex_idx) THEN
                FLAG2= .TRUE.
                IF(inter_edge_idx==1) THEN
                  FLAG1=.FALSE.
                ELSE 
                  FLAG1=.TRUE.
                END IF 
                EXIT
              END IF 
            END DO  

            IF(FLAG2 .EQV. .TRUE.) EXIT

          END DO

          IF(FLAG1 .EQV. .FALSE.) THEN

            NEW_TO_OLD_INDEX_MAPPING(OLD_VERTEX_idx,2)=new_vertex_idx
            new_vertex_idx = new_vertex_idx + 1 

          ELSE 

            NEW_TO_OLD_INDEX_MAPPING(old_vertex_idx,2)= &
              NEW_TO_OLD_INDEX_MAPPING(COUPLED_DECOMPOSITION%INTER_EDGES(1,proc_idx),2)

          END IF
        END DO ! OLD_VERTEX_idx
      ELSE
      CALL FlagError("Coupled mesh decomposition is not associated.",ERR,ERROR,*999)
      ENDIF 
    ELSE
      CALL FlagError("Coupled Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    EXITS("COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING")
    RETURN
!ERRORSEXITS("COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING",ERR,ERROR)
999 RETURN 1
  END SUBROUTINE COUPLED_DECOMPOSITION_NEW_TO_OLD_VERTEX_MAPPING


! =========================================================================================================!
!THe following subroutine build adjacencies of the new coupled mesh graph G_{i,merged}.
! =========================================================================================================!
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
    LOGICAL :: FLAG1, FLAG2

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
!ERRORSEXITS("COUPLED_DECOMPOSITION_GET_NEW_GRAPH",ERR,ERROR)
999 RETURN 1 

  END SUBROUTINE COUPLED_DECOMPOSITION_GET_NEW_GRAPH

!=============================================================================================================

  integer FUNCTION count_integer(Array,NumberToBeCOunted)

    integer, INTENT(IN)   :: Array(:), NumberToBeCOunted


    ! local variable
    integer               :: counter
    
    ! action starts here
    count_integer                   = 0 
    DO counter= 1,SIZE(Array)
       
      IF(Array(counter) == NumberToBeCOunted) & 
        & count_integer = count_integer + 1 

    END  DO
  


  END  function count_integer


  
END MODULE MESH_ROUTINES



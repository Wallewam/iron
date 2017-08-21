  !>Calculates the local/global node and dof mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENT_NEW_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,no_computational_node,no_ghost_element,adjacent_element,ghost_element, &
    & NUMBER_OF_NODES_PER_DOMAIN,domain_idx,domain_idx2,domain_no,Element_idx,derivative_idx,version_idx,ny,NUMBER_OF_DOMAINS, &
    & MAX_NUMBER_DOMAINS,NUMBER_OF_GHOST_ELEMENTS,my_computational_node_number,number_computational_nodes,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_ELEMENT_NUMBERS(:),LOCAL_DOF_NUMBERS(:),ELEMENT_COUNT(:),NUMBER_INTERNAL_ELEMENTS(:), &
    & NUMBER_BOUNDARY_ELEMENTS(:)
    INTEGER(INTG), ALLOCATABLE :: DOMAINS(:),ALL_DOMAINS(:),GHOST_ELEMENTS(:)
    LOGICAL :: BOUNDARY_DOMAIN
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST,ALL_ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_ELEMENTS_LIST(:)
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_ELEMENT_NEW_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
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
                  my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  !Calculate the local and global numbers and set up the mappings
                  ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
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

    
    EXITS("DOMAIN_MAPPINGS_ELEMENT_NEW_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ALL_DOMAINS)) DEALLOCATE(ALL_DOMAINS)
    IF(ALLOCATED(GHOST_ELEMENTS)) DEALLOCATE(GHOST_ELEMENTS)
    IF(ALLOCATED(NUMBER_INTERNAL_ELEMENTS)) DEALLOCATE(NUMBER_INTERNAL_ELEMENTS)
    IF(ALLOCATED(NUMBER_BOUNDARY_ELEMENTS)) DEALLOCATE(NUMBER_BOUNDARY_ELEMENTS)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENT_NEW_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENT_NEW_CALCULATE
  






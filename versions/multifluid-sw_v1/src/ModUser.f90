!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModUser
  ! User module for Mercury, modified by Alex Cushen to support multifluid (H+ and He2+)
  ! Must compile with ModEquationMultiIonPeMercury equation set.
  ! Version 2.1 can work with both Cartesian and spherical grid.
  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest, iProc, nDim, &
       nI, nJ, nK, nG, IsCartesian, Xyz_DGB, FaceNormal_DDFB, CellFace_DFB
  use ModMain,       ONLY: UseSolidState,Coord1MinBc_,Coord1MaxBc_,solidBc_, &
       UseHyperbolicDivb, UseResistivePlanet, FaceBCType
  use ModVarIndexes
  use ModPhysics, ONLY: Si2No_V, No2Si_V, UnitRho_, UnitN_, UnitU_, &
       UnitP_, UnitX_, UnitT_, BodyP_I, &
       FaceState_VI, CellState_VI
  use ModMultiFluid, ONLY: iRho, iRhoUx, iRhoUz, iUx, iUz, iRho, iP, &
       iUx_I, iUz_I, nFluid, select_fluid, MassIon_I
  use ModNumConst, ONLY: cPi
  use ModConst, ONLY: cBoltzmann
  use ModUserEmpty,                          &
       IMPLEMENTED1 => user_init_session,    &
       IMPLEMENTED2 => user_set_ics,         &
       IMPLEMENTED3 => user_set_resistivity, &
       IMPLEMENTED4 => user_read_inputs,     &
       IMPLEMENTED5 => user_set_face_boundary,&
       IMPLEMENTED6 => user_set_cell_boundary,&
       IMPLEMENTED7 => user_calc_sources_expl

  include 'user_module.h' ! list of public methods

  real,             parameter :: VersionUserModule = 4.1
  character(len=*), parameter :: NameUserFile = "ModUserMercurySolid"
  character(len=*), parameter :: NameUserModule = "MercuryMultiFluid"

  real :: PlanetRadius=1.
  character(len=10) :: InitType = 'B1U1'
  integer :: nSolidBcType = 1
  !real :: HpHe2pRatio = 0.9
  !real :: HpSWTemp = 8.70e4
  !real :: He2pSWTemp = 8.70e4

  integer :: iLayer
  integer :: nLayer=0 ! Number of points in planet resistivity profile
  real, allocatable :: RadiusSi_I(:),Radius_I(:),&
       ResistivitySi_I(:),Resistivity_I(:)
  real, allocatable :: Resistivity(:)

  real :: nNeutral0, ScaleHeight, ScaleHeightIn, UUp
  real :: IonizationRate, CrossSection
  logical :: UseMassLoading = .false.

  ! The lower bound of pe/p at inner boundary when the electron
  ! pressure equation is used.
  real :: RatioPe2P = 0

contains
  !============================================================================
  subroutine user_read_inputs

    use ModReadParam

    character(len=100) :: NameCommand
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

	   !case('#PROTONHELIUMRATIO')
       !   call read_var('HpHe2pRatio', HpHe2pRatio)

       !case('#SWTEMPS')
       !   call read_var('HpSWTemp', HpSWTemp)
       !   call read_var('He2pSWTemp', He2pSWTemp)

       case('#INITTYPE')
          call read_var('InitType', InitType)

       case("#RESISTIVEPLANET")
          UseResistivePlanet = .true.
          call read_var('nResistivPoints', nLayer)
          if(nLayer == 1) then
             write(*,*) ' We need minimum 2 points for including resistivity &
                  &profile'
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if

          if(nLayer > 1) then
             allocate(Resistivity(nLayer-1),&
                  RadiusSi_I(nLayer),&
                  Radius_I(nLayer), &
                  ResistivitySi_I(nLayer),&
                  Resistivity_I(nLayer))

             do iLayer=1,nLayer
                call read_var('Radius', Radius_I(iLayer))
                call read_var('Resistivity', ResistivitySi_I(iLayer))
             end do

             ! Check values
             do iLayer=2,nLayer
                if(Radius_I(iLayer-1) < &
                     Radius_I(iLayer)) then
                   write(*,*) 'ERROR: Should be decreasing Radius.'
                   call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
                end if
             end do
          end if

       case("#MASSLOADING")
          UseMassLoading = .true.
          if(UseMassLoading) then
             call read_var('NeutralDensity (cm^-3)', nNeutral0)
             call read_var('ScaleHeight (km)', ScaleHeightIn)
             call read_var('IonizationRate (1/s)', IonizationRate)
             call read_var('CollisionCrossSection (cm2)', CrossSection)
          end if

       case("#SOLIDBOUNDARY")
          call read_var('nSolidBcType', nSolidBcType)

       case("#PEOVERP")
          call read_var('RatioPe2P', RatioPe2P)

       case('#USERINPUTEND')
          EXIT
       case default
          if(iProc==0) then
             write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) ' *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

    call test_stop(NameSub, DoTest)
  end subroutine user_read_inputs
  !============================================================================

  subroutine user_init_session

    use CON_planet,     ONLY: RadiusPlanet, MassPlanet
    use ModIO,          ONLY: write_myname
    use ModResistivity, ONLY: Si2NoEta
    use ModGeometry,    ONLY: TypeGeometry
    use ModAdvance,     ONLY: State_VGB
    use ModMain,        ONLY: xMinBc_, zMaxBc_

    !type(FaceBCType), intent(inout) :: FBC

    real :: UxUp,UyUp,UzUp
    logical:: DoTest

    integer :: iBoundary
    character(len=*), parameter:: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------
    !associate( VarsGhostFace_V => FBC%VarsGhostFace_V, &
    !           VarsTrueFace_V => FBC%VarsTrueFace_V, &
    !           FaceCoords_D => FBC%FaceCoords_D, &
    !           B0Face_D => FBC%B0Face_D, &
    !           iBoundary => FBC%iBoundary, &
    !           iFace => FBC%iFace, jFace => FBC%jFace, kFace => FBC%kFace, &
    !           iBlock => FBC%iBlockBc, iSide => FBC%iSide )
    
    call test_start(NameSub, DoTest)

    if(iProc==0)then
       select case(TypeGeometry)
       case('spherical_lnr','spherical')
          write(*,*) 'Using spherical grid...'
       case('cartesian','rotatedcartesian')
          write(*,*) 'Using Cartesian grid...'
       case default
          call stop_mpi('ERROR: need spherical/Cartesian grid!')
       end select
    end if

    ! Set upstream H+ and He2+ density and pressure
    !do iBoundary=xMinBc_, zMaxBc_
       ! I believe there are mass densities instead of number densities. Not 100% percent sure. --Yuxi.
       ! FaceState_VI(rho_,iBoundary)   = SWRho_total
       ! Using the user input for Hp and He2p temperatures to overwrite with p=nkT. 1e15 is from cc-> m^3 and Pa-> nPa conversion
       !FaceState_VI(Rho_,iBoundary) = FaceState_VI(rho_,iBoundary)*HpHe2pRatio
       !FaceState_VI(He2pRho_,iBoundary) = FaceState_VI(rho_,iBoundary)*(1-HpHe2pRatio)
       !FaceState_VI(P_,iBoundary) = FaceState_VI(rho_,iBoundary)*HpSWTemp*cBoltzmann * (1e15)
       !FaceState_VI(He2pP_,iBoundary) = FaceState_VI(He2pRho_,iBoundary)*He2pSWTemp*cBoltzmann * (1e15) / 4
    !end do

    
    if(nLayer > 1) then
       RadiusSi_I = Radius_I*No2Si_V(UnitX_)
       Resistivity_I = ResistivitySi_I*Si2NoEta
       do iLayer=2,nLayer
          Resistivity(iLayer-1) = &
               (Resistivity_I(iLayer) - Resistivity_I(iLayer-1))/&
               (Radius_I(iLayer) - Radius_I(iLayer-1))
       end do
    end if

    ! Convert to dimensionless quantities
    if(UseMassLoading) &
         ScaleHeight = ScaleHeightIn*1e3*Si2No_V(UnitX_)

    ! Get upstream plasma velocity
    UxUp = CellState_VI(rhoux_,Coord1MaxBc_) / CellState_VI(rho_,Coord1MaxBc_)
    UyUp = CellState_VI(rhouy_,Coord1MaxBc_) / CellState_VI(rho_,Coord1MaxBc_)
    UzUp = CellState_VI(rhouz_,Coord1MaxBc_) / CellState_VI(rho_,Coord1MaxBc_)
    UUp = sqrt(UxUp**2 + UyUp**2 + UzUp)


    if(iProc==0) then
       call write_myname
       write(*,*) ''
       write(*,*) '   Resistive Planet Model'
       write(*,*) '   ----------------------'
       write(*,*) ''
       write(*,"(A29,E10.3)") '  Planet radius   [m]      = ',RadiusPlanet
       if(nLayer > 0 ) then
          write(*,*) ''
          write(*,*) '   |-------- Planet Resistivity Profile -----|'
          write(*,*) '       Radius(SI)            Resistivity(SI)'
          do iLayer =1,nLayer
             write(*,"(A7,E10.3,A15,E10.3)") " ",RadiusSi_I(iLayer)," ",&
                  ResistivitySi_I(iLayer)
          end do
       else
          write(*,*) 'Conducting Planet (eta = 0)'
       end if
       write(*,*) ''
       if(UseMassLoading) then
          write(*,'(A24,E8.2,A)') 'neutral density =', nNeutral0, ' [cm^-3]'
          write(*,'(A24,E8.2,A)') 'scale height    =', ScaleHeight*No2Si_V(UnitX_)*1e-3, ' [km]'
          write(*,'(A24,E8.2,A)') 'ionization rate =', IonizationRate, ' [s^-1]'
          write(*,'(A24,E8.2,A)') 'cross section =', CrossSection, ' [cm^2]'
          write(*,'(A24,E8.2,A)') 'background density =', CellState_VI(rho_,Coord1MaxBc_), ' [amu/cc]'
          write(*,'(A24,E8.2,A)') 'upstream velocity =', UUp*No2Si_V(UnitU_), ' [m/s]'
          write(*,*) ''
       end if
    end if
    call test_stop(NameSub, DoTest)
    !end associate
  end subroutine user_init_session
  !============================================================================

  subroutine user_set_ics(iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: rMin_B, r_GB
    use ModMain,       ONLY: xMinBc_    

    integer, intent(in) :: iBlock

    integer :: i,j,k,iFluid
    integer :: Bc_
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(IsCartesian)then
       Bc_ = Coord1MinBc_
    else
       Bc_ = Coord1MaxBc_
    end if

    
    !State_VGB(Rho_,:,:,:,iBlock) = FaceState_VI(rho_,xMinBc_)
    !State_VGB(Rho_,:,:,:,iBlock) = FaceState_VI(Rho_,xMinBc_)
    !State_VGB(He2pRho_,:,:,:,iBlock) = FaceState_VI(He2pRho_,xMinBc_)
    
    select case(InitType)
    case('B1U1')
       ! Uniform density
       State_VGB(rho_,:,:,:,iBlock) = CellState_VI(rho_,Bc_)
       ! Uniform pressure
       State_VGB(p_ ,:,:,:,iBlock)  = CellState_VI(p_,Bc_)
       State_VGB(pe_,:,:,:,iBlock)  = CellState_VI(pe_,Bc_)

       ! Initialize mantle region
       if(rMin_B(iBlock) <= PlanetRadius) then
          do iFluid = 1, nFluid
             call select_fluid(iFluid)
             do k=1,nK; do j=1,nJ; do i=1,nI
                if(r_GB(i,j,k,iBlock) > PlanetRadius) CYCLE
                State_VGB(iRho,i,j,k,iBlock) = FaceState_VI(rho_,solidBc_)
                State_VGB(iP,i,j,k,iBlock)   = FaceState_VI(p_,solidBc_)
                State_VGB(pe_,i,j,k,iBlock)  = FaceState_VI(pe_,solidBc_)
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             end do; end do; end do
          end do
       end if

       ! Magnetic field inside the mantle is dipole field
       ! Note that here we store B1
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(r_GB(i,j,k,iBlock) > PlanetRadius)then
             State_VGB(bx_:bz_,i,j,k,iBlock) = &
                  CellState_VI(bx_:bz_,Bc_)
          else
             State_VGB(bx_:bz_,i,j,k,iBlock) = 0.0
          end if
       end do; end do; end do

       ! Init for velocity
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iFluid=1,nFluid
             call select_fluid(iFluid)
             if(r_GB(i,j,k,iBlock) < PlanetRadius)then
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             else
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     CellState_VI(iRhoUx:iRhoUz,Bc_)
             end if
          end do
       end do; end do; end do

    case('B1U01')
       ! Uniform density
       State_VGB(rho_,:,:,:,iBlock) = CellState_VI(rho_,Bc_)
       ! Uniform pressure
       State_VGB(p_ ,:,:,:,iBlock)  = CellState_VI(p_,Bc_)
       State_VGB(pe_,:,:,:,iBlock)  = CellState_VI(pe_,Bc_)

       ! Initialize mantle region
       if(rMin_B(iBlock) <= PlanetRadius) then
          do iFluid = 1, nFluid
             call select_fluid(iFluid)
             do k=1,nK; do j=1,nJ; do i=1,nI
                if(r_GB(i,j,k,iBlock) > PlanetRadius) CYCLE
                State_VGB(iRho,i,j,k,iBlock) = FaceState_VI(rho_,solidBc_)
                State_VGB(iP,i,j,k,iBlock)   = FaceState_VI(p_,solidBc_)
                State_VGB(pe_,i,j,k,iBlock)  = FaceState_VI(pe_,solidBc_)
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             end do; end do; end do
          end do
       end if

       ! Magnetic field inside the mantle is dipole
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(r_GB(i,j,k,iBlock) > PlanetRadius)then
             State_VGB(bx_:bz_,i,j,k,iBlock) = &
                  CellState_VI(bx_:bz_,Bc_)
          else
             State_VGB(bx_:bz_,i,j,k,iBlock) = 0.0
          end if
       end do; end do; end do

       ! Smooth init for velocity
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iFluid=1,nFluid
             call select_fluid(iFluid)
             if(r_GB(i,j,k,iBlock) < 2.5*PlanetRadius)then
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             elseif(r_GB(i,j,k,iBlock) < 5.0*PlanetRadius)then
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     CellState_VI(iRhoUx:iRhoUz,Bc_) *&
                     (r_GB(i,j,k,iBlock)/2.5 - 1)
             else
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     CellState_VI(iRhoUx:iRhoUz,Bc_)
             end if
          end do
       end do; end do; end do

    case('B0U0') ! Upstream propagation
       ! Magnetic field starts from dipole
       do k=1,nK; do j=1,nJ; do i=1,nI
          State_VGB(bx_:bz_,i,j,k,iBlock) = 0.0
       end do; end do; end do

       do k=1,nK; do j=1,nJ; do i=1,nI
          do iFluid=1,nFluid
             call select_fluid(iFluid)
             ! Velocity starts from zeros
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             ! Uniform density
             State_VGB(iRho,i,j,k,iBlock) = CellState_VI(rho_,Bc_)
             ! Uniform pressure
             State_VGB(p_,i,j,k,iBlock) = CellState_VI(p_,Bc_)
             State_VGB(pe_,i,j,k,iBlock) = CellState_VI(pe_,Bc_)
          end do
       end do; end do; end do

    case default
       call stop_MPI('ERROR: unknown initialization type!')
    end select

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_ics
  !============================================================================

  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, IsFound)

    use ModAdvance,      ONLY: State_VGB
    use ModCellBoundary, ONLY: iMin, iMax, jMin, jMax, kMin, kMax

    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i,j,k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_cell_boundary'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! This is only used for spherical grid with no outer face BC applied.
    ! If outer face BC is used, set the second parameter in command
    ! #OUTERBOUNDARY to none to avoid calling this function.

    ! Outer boundary condition
    if(iSide == 2)then
       do k=kMin,kMax; do j=jMin,jMax; do i=iMin,iMax
          if(Xyz_DGB(1,i,j,k,iBlock) < 0)then
             ! Upstream fixed
             State_VGB(:,i,j,k,iBlock) = CellState_VI(:,iSide)
          else
             ! Downstream float
             State_VGB(:,i,j,k,iBlock) = State_VGB(:,nI,j,k,iBlock)
          end if
       end do; end do; end do

       IsFound = .true.
       RETURN
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_cell_boundary
  !============================================================================

  subroutine user_set_face_boundary(FBC)

    use ModAdvance,      ONLY: State_VGB

    type(FaceBCType), intent(inout) :: FBC

    integer :: i, j, k, iFluid, iDir
    real :: bUnit_D(3), Normal_D(3), UdotR, r2Inv

    character(len=*), parameter:: NameSub = 'user_set_face_boundary'
    !--------------------------------------------------------------------------
    associate( VarsGhostFace_V => FBC%VarsGhostFace_V, &
               VarsTrueFace_V => FBC%VarsTrueFace_V, &
               FaceCoords_D => FBC%FaceCoords_D, &
               B0Face_D => FBC%B0Face_D, &
               iBoundary => FBC%iBoundary, &
               iFace => FBC%iFace, jFace => FBC%jFace, kFace => FBC%kFace, &
               iBlock => FBC%iBlockBc, iSide => FBC%iSide )

    ! Copy everything from physical face to ghost face
    ! rho, u, b, p, hyp
    VarsGhostFace_V = VarsTrueFace_V


    if(RatioPe2P>0)VarsGhostFace_V(Pe_) = &
         max(VarsTrueFace_V(Pe_), RatioPe2P*VarsTrueFace_V(P_))


    ! For 1st order floating BC of magnetic field, the true face value is
    ! the same as the first layer of physical cell center values.
    ! For 2nd order floating BC of magnetic field, extrapolation is needed,
    ! which requires two layers of physical cells' center values.

    if(nSolidBcType==1)then
       ! Fixed rho
       ! VarsGhostFace_V(rho_) = FaceState_VI(rho_,iBoundary)

       ! Fixed pressure
       ! VarsGhostFace_V(p_)   = FaceState_VI(p_,iBoundary)
       ! VarsGhostFace_V(pe_)  = FaceState_VI(pe_,iBoundary)

       ! Use true face B to get B unit direction
       bUnit_D = B0Face_D + VarsTrueFace_V(bx_:bz_)

       ! Boundary velocity based on B field geometry
       bUnit_D = bUnit_D/max(1e-30, sqrt(sum(bUnit_D**2)))
       VarsGhostFace_V(RhoUx_:RhoUz_) = VarsTrueFace_V(RhoUx_:RhoUz_) - &
            2*sum(bUnit_D*VarsTrueFace_V(RhoUx_:RhoUz_))*bUnit_D
       
       ! For outflow reflect radial velocity: uG = u - 2*(u.r)*r/r^2 
       UdotR = sum(VarsTrueFace_V(RhoUx_:RhoUz_)*FaceCoords_D)
       if(UdotR > 0.0) &
            VarsGhostFace_V(rho_) = FaceState_VI(rho_,iBoundary)

    elseif(nSolidBcType==2)then ! Absorb U
       ! Calculate 1/r^2
       r2Inv = 1.0/sum(FaceCoords_D**2)
       ! For outflow reflect radial velocity: uG = u - 2*(u.r)*r/r^2
       UdotR = sum(VarsTrueFace_V(RhoUx_:RhoUz_)*FaceCoords_D)
       if(UdotR > 0.0) &
            VarsGhostFace_V(RhoUx_:RhoUz_) = VarsTrueFace_V(RhoUx_:RhoUz_) &
            - 2*UdotR*r2Inv*FaceCoords_D
    elseif(nSolidBcType==3)then ! Reflect U
       iDir = (iSide+1)/2
       Normal_D(1:nDim) = FaceNormal_DDFB(:,iDir,iFace-1,jFace,kFace,iBlock)&
            / CellFace_DFB(iDir,iFace-1,jFace,kFace,iBlock)
       VarsGhostFace_V(RhoUx_:RhoUz_) = VarsTrueFace_V(RhoUx_:RhoUz_) &
            - 2*sum(VarsTrueFace_V(RhoUx_:RhoUz_)*Normal_D)*Normal_D
    end if

    end associate
  end subroutine user_set_face_boundary
  !============================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)

    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModGeometry,   ONLY: r_GB, rMin_B
    use ModResistivity, ONLY: Eta0

    integer, intent(in) :: iBlock
    real, intent(out) :: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)

    integer ::i,j,k,iLayer
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    Eta_G = Eta0

    if(nLayer < 2) RETURN
    if(rMin_B(iBlock) > Radius_I(1)) RETURN

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do iLayer=nLayer-1,1,-1
          if(r_GB(i,j,k,iBlock) < Radius_I(iLayer+1) ) CYCLE
          if(r_GB(i,j,k,iBlock) > Radius_I(iLayer) ) CYCLE
          ! Avoid eta jumps adding Eta_G
          Eta_G(i,j,k) = Eta_G(i,j,k) + Resistivity_I(iLayer+1)+&
               (r_GB(i,j,k,iBlock)-Radius_I(iLayer+1))* &
               Resistivity(iLayer)
       end do
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_resistivity
  !============================================================================
  subroutine user_calc_sources_expl(iBlock)

    use ModAdvance, ONLY: Source_VC, State_VGB
    use ModGeometry, ONLY: Xyz_DGB, r_GB, rMin_B
    use ModConst, ONLY: cBoltzmann

    integer, intent(in) :: iBlock

    real, dimension(1:nI,1:nJ,1:nK):: &
         SRho_C,SRhoUx_C,SRhoUy_C,SRhoUz_C,SP_C,SPe_C,SE_C

    ! Variable meanings:
    !   SRho_C: Source terms for the continuity equation
    !   SE_C,SP_C: Source terms for the energy (conservative) and
    !          presure (primative) equations
    !   SRhoUx_C,SRhoUy_C,SRhoUz_C: Source terms for the momentum equation
    !   SBx_C,SBy_C,SBz_C: Source terms for the magnetic field equations

    ! Local variables

    integer :: i,j,k
    ! To be consistent with Galileo calculations
    real, parameter :: IonMass = 32.0
    real, parameter :: rMax = 1.3, rMin = 1.0
    real :: Ux,Uy,Uz,U2
    real :: rhodot,nNeutral,alpha,Te
    real :: CollisionRate, RecombRate

    logical :: DoTest
    character(len=*), parameter:: NameSub = 'user_calc_sources_expl'
    !--------------------------------------------------------------------------
    call test_start(NameSub,DoTest,iBlock)

    if(.not. UseMassLoading) RETURN
    if(rMin_B(iBlock) > rMax) RETURN

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

      ! Calculate only in the region where the contribution is not
      ! negligible.

      if((r_GB(i,j,k,iBlock) < rMax).and.(r_GB(i,j,k,iBlock) > rMin)) then

         Ux = State_VGB(rhoUx_,i,j,k,iBlock) / State_VGB(rho_,i,j,k,iBlock)
         Uy = State_VGB(rhoUy_,i,j,k,iBlock) / State_VGB(rho_,i,j,k,iBlock)
         Uz = State_VGB(rhoUz_,i,j,k,iBlock) / State_VGB(rho_,i,j,k,iBlock)
         U2 = sqrt(Ux**2 + Uy**2 + Uz**2)

         ! Neutral mass density
         nNeutral = nNeutral0*exp((1.0 - r_GB(i,j,k,iBlock))/ScaleHeight)

         ! Collision loss
         CollisionRate = nNeutral*CrossSection*UUp*1e2/Si2No_V(UnitX_)

         ! Ionization
         rhodot = nNeutral*IonizationRate*IonMass

         ! Recombination
         ! Te = State_VGB(pe_,i,j,k,iBlock) / State_VGB(rho_,i,j,k,iBlock)* &
         !     MassIon_I(1)/cBoltzmann*No2Si_V(UnitT_)
         ! alpha = 1.6e-13*(300/Te)**0.55*Si2No_V(UnitX_)**3
         alpha = 7.8e-14*1e6 ! [cm^3/s]
         if(State_VGB(rho_,i,j,k,iBlock) > &
              CellState_VI(rho_,Coord1MaxBc_)) then
            RecombRate = alpha* &
                 (State_VGB(rho_,i,j,k,iBlock) - &
                 CellState_VI(rho_,Coord1MaxBc_))
         else
            RecombRate = 0.0
         end if

         ! Load the source terms into the right hand side

         SRho_C(i,j,k) = rhodot-RecombRate*State_VGB(rho_,i,j,k,iBlock)

         SRhoUx_C(i,j,k) = &
              - (CollisionRate + RecombRate)*State_VGB(rhoux_,i,j,k,iBlock)

         SRhoUy_C(i,j,k) = &
              - (CollisionRate + RecombRate)*State_VGB(rhouy_,i,j,k,iBlock)

         SRhoUz_C(i,j,k) = &
              - (CollisionRate + RecombRate)*State_VGB(rhouz_,i,j,k,iBlock)

         SE_C(i,j,k) = &
              - 0.5*(CollisionRate + RecombRate)*U2*State_VGB(rho_,i,j,k,iBlock) &
              + 1.5*(CollisionRate - RecombRate)*State_VGB(p_,i,j,k,iBlock)

         SP_C(i,j,k) = &
              + 0.5*(rhodot + CollisionRate*State_VGB(rho_,i,j,k,iBlock))*U2 &
              - 1.5*RecombRate*State_VGB(p_,i,j,k,iBlock)

         ! hyzhou: there are still issues with the pressure terms!

         ! Todo: adding electron pressure source terms
         ! SPe_C(i,j,k) = SPe_C(i,j,k)
      else
         SRho_C  (i,j,k) = 0.0
         SRhoUx_C(i,j,k) = 0.0
         SRhoUy_C(i,j,k) = 0.0
         SRhoUz_C(i,j,k) = 0.0
         SP_C    (i,j,k) = 0.0
         SPe_C   (i,j,k) = 0.0
         SE_C    (i,j,k) = 0.0
      end if ! end source term calc.
   end do; end do; end do ! end the i,j,k loops

   if(DoTest)then
      write(*,*) 'max(SRho_C) = ', maxval(SRho_C)
      write(*,*) 'max(State_VGB(rho_)) = ', maxval(State_VGB(rho_,:,:,:,iBlock))
      write(*,*) ''
   end if

   Source_VC(rho_   ,:,:,:) = SRho_C   + Source_VC(rho_   ,:,:,:)
   Source_VC(rhoux_ ,:,:,:) = SRhoUx_C + Source_VC(rhoux_ ,:,:,:)
   Source_VC(rhouy_ ,:,:,:) = SRhoUy_C + Source_VC(rhouy_ ,:,:,:)
   Source_VC(rhouz_ ,:,:,:) = SRhoUz_C + Source_VC(rhouz_ ,:,:,:)
   Source_VC(p_     ,:,:,:) = SP_C     + Source_VC(p_     ,:,:,:)
   ! Source_VC(pe_    ,:,:,:) = SPe_C    + Source_VC(pe_    ,:,:,:)
   Source_VC(Energy_,:,:,:) = SE_C     + Source_VC(Energy_,:,:,:)

   call test_stop(NameSub,DoTest,iBlock)

end subroutine user_calc_sources_expl
  !============================================================================

end module ModUser
!==============================================================================


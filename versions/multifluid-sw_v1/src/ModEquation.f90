!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!#NOTPUBLIC  email:xzjia@umich.edu  expires:12/31/2099
module ModVarIndexes

  use ModExtraVariables,        &
       Redefine => iPparIon_I,  &
       Redefine1 => Pe_, &
       Redefine2 => Hyp_

  implicit none

  save

  character (len=*), parameter :: &
       NameEquationFile = "ModEquationMultiIonPeMercury.f90"

  ! This equation module contains the standard MHD equations
  ! with electron pressure
  character (len=*), parameter :: &
       NameEquation = 'MultiFluid: H+, He2+, Hyp, and Pe DEVEL V1.4'

  integer, parameter :: nVar = 15

  integer, parameter :: nFluid    = 2
  integer, parameter :: nIonFluid = 2        ! Last individual ion fluid
  logical, parameter :: IsMhd     = .false.   ! First total ion fluid obeys MHD
  real               :: MassFluid_I(1:nFluid) = [1.0, 4.0]

  ! Fluids: total fluid, solar wind protons, cometary water ions
  character (len=4), parameter :: NameFluid_I(1:nFluid) = [ 'Hp  ', 'He2p' ]

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_        =  1, HpRho_ = 1, &
       RhoUx_      =  2, Ux_ = 2, HpRhoUx_ = 2, HpUx_ = 2, &
       RhoUy_      =  3, Uy_ = 3, HpRhoUy_ = 3, HpUy_ = 3, &
       RhoUz_      =  4, Uz_ = 4, HpRhoUz_ = 4, HpUz_ = 4, &
       Bx_         =  5, &
       By_         =  6, &
       Bz_         =  7, &
       Pe_         =  8, &
       P_          =  9, HpP_ = 9, &
       He2pRho_    = 10, &
       He2pRhoUx_  = 11, He2pUx_ = 11, &
       He2pRhoUy_  = 12, He2pUy_ = 12, &
       He2pRhoUz_  = 13, He2pUz_ = 13, &
       He2pP_      = 14, &
       Hyp_        = 15, &
       Energy_     = nVar+1, HpEnergy_ = nVar+1, &
       He2pEnergy_ = nVar+2

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)  =[HpRho_,   He2pRho_ ] ,&
       iRhoUx_I(nFluid)=[HpRhoUx_, He2pRhoUx_ ],&
       iRhoUy_I(nFluid)=[HpRhoUy_, He2pRhoUy_ ],&
       iRhoUz_I(nFluid)=[HpRhoUz_, He2pRhoUz_ ],&
       iP_I(nFluid)    =[HpP_,     He2pP_ ]

  integer, parameter :: iPparIon_I(nIonFluid) = [1,2]

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+nFluid) = [ &
       1.0, & ! HpRho_
       0.0, & ! HpRhoUx_
       0.0, & ! HpRhoUy_
       0.0, & ! HpRhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! Pe_
       1.0, & ! HpP_
       1.0, & ! He2pRho_
       0.0, & ! He2pRhoUx_
       0.0, & ! He2pRhoUy_
       0.0, & ! He2pRhoUz_
       1.0, & ! He2pP_
       0.0, & ! Hyp_
       1.0, & ! HpEnergy_
       1.0 ]  ! He2pEnergy_

  ! The names of the variables used in i/o
  character(len=7) :: NameVar_V(nVar+nFluid) = [ &
       'HpRho  ', & ! HpRho_
       'HpMx   ', & ! HpRhoUx_
       'HpMy   ', & ! HpRhoUy_
       'HpMz   ', & ! HpRhoUz_
       'Bx     ', & ! Bx_
       'By     ', & ! By_
       'Bz     ', & ! Bz_
       'Pe     ', & ! Pe_
       'HpP    ', & ! HpP_
       'He2pRho', & ! He2pRho_
       'He2pMx ', & ! He2pRhoUx_
       'He2pMy ', & ! He2pRhoUy_
       'He2pMz ', & ! He2pRhoUz_
       'He2pP  ', & ! He2pP_
       'Hyp    ', & ! Hyp_
       'HpE    ', & ! HpEnergy_
       'He2pE  ' ]  ! He2pEnergy_

  ! The space separated list of nVar primitive variables for plotting
    character(len=*), parameter :: NamePrimitiveVar = &
       'HpRho HpUx HpUy HpUz Bx By Bz Pe HpP '// &
       'He2pRho He2pUx He2pUy He2pUz He2pP '

  ! There are no extra scalars (Pe has its own flux)
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

end module ModVarIndexes
!==============================================================================


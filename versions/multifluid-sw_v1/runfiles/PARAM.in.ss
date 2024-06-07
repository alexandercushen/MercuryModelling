#ECHO
T

#DESCRIPTION
Multifluid solar wind (H+ and He2+)
Steady-state.
Sotuhward IMF standard run settings.

#PLANET
Mercury                NamePlanet

#ROTATION
F               UseRotation

#MAGNETICCENTER
0.0              MagCenterX
0.0             MagCenterY
0.2             MagCenterZ

#TIMEACCURATE
F               DoTimeAccurate

#LAYOUT
GM    0     -1   1   1
PC    0	    -1	 1   1	

#CHECKSTOP
T               DoCheckStop
1000            DnCheckStop
-0.015          DtCheckStop

#SAVERESTART
T                   DoSaveRestart
5000               DnSaveRestart
-5                 DtSaveRestart

#CPUTIMEMAX
70000.0         CpuTimeMax

#COMPONENT
PC                      NameComp
F                       UseComp

#COUPLETIME
GM                      NameComp
F                       DoCoupleOnTime

#BEGIN_COMP GM ---------------------------------------

INCLUDE
GM/restartIN/restart.H

#HYPERBOLICDIVB
T                       UseHyperbolicDivb
500.0                   SpeedHypDim
0.2                     HypDecay

#COORDSYSTEM
GSE                     TypeCoordinate

BODY
T                       UseBody
0.8                     rBody      [rPlanet]
3.5                     rCurrents  [rPlanet]
39                      BodyNDim   [/cc]
2.32e5                  BodyTDim   [K]

SOLARWIND
36.0	                SwNDim   [/cc]
8.70e4                  SwTDim   [K] (3.8nPa/(1.38e-23*2/cc))/1.2
-500.0                  SwUxDim  [km/s]
0.0                     SwUyDim  [km/s]
0.0                     SwUzDim  [km/s]
0.0                     SwBxDim  [nT]
0.0                     SwByDim  [nT]
23.0                    SwBzDim  [nT]

#PLASMA
1.0                     FluidMass [amu]
4.0                     FluidMass [amu] (He2+)
1.0                     AverageIonCharge [e]
2.0                     AverageIonCharge [e] (He2+)
1.0                 ElectronTemperatureRatio

TEST
krylov

#RESISTIVITY
T                       UseResistivity
user                    TypeResistivity
0.0                     Eta0Si

KRYLOV
GMRES                   TypeKrylov  (GMRES, BICGSTAB, CG)
nul                     TypeInitKrylov (nul, old, explicit, scaled)
0.001                   ErrorMaxKrylov
200                     MaxMatvecKrylov

#MINIMUMPRESSURE
1e-11                   pMinDim
1e-11
1e-11                   PeMinDim for electron pressure

#MINIMUMDENSITY
1e-5                     RhoMinDim
1e-5

#NONCONSERVATIVE
F                       UseNonConservative

#CONSERVATIVECRITERIA
0                       nConservCrit

#RESTARTOUTFILE
one                     TypeRestartOutFile

! Grid structure info
#INCLUDE
Grid

----------------BC-----------------
#OUTERBOUNDARY
fixedb1                 TypeCellBc1
none                    TypeCellBc2
periodic                TypeCellBc3
periodic                TypeCellBc4
periodic                TypeCellBc5
periodic                TypeCellBc6

#BOXBOUNDARY
float                   TypeBcXmin
fixed                   TypeBcXmax
float                   TypeBcYmin
float                   TypeBcYmax
float                   TypeBcZmin
float                   TypeBcZmax

#BOUNDARYSTATE
xmaxbox            StringBoundary
15.0                    BoundaryStateDim_V HpRho
-500.0                  BoundaryStateDim_V HpUx
0.0                     BoundaryStateDim_V HpUy
0.0                     BoundaryStateDim_V HpUz
0.0                     BoundaryStateDim_V Bx
0.0                     BoundaryStateDim_V By
-25.0                    BoundaryStateDim_V Bz
0.02                   BoundaryStateDim_V Pe
0.02                    BoundaryStateDim_V HpP
15.0            BoundaryStateDim_V He2pRho
-500.0            BoundaryStateDim_V He2pUx
0.0                     BoundaryStateDim_V He2pUy
0.0                     BoundaryStateDim_V He2pUz
0.02                     BoundaryStateDim_V He2pP
0.0                     BoundaryStateDim_V Hyp

#BOUNDARYSTATE
coord1min solid         StringBoundary
5.0                     BoundaryStateDim_V HpRho
0.0                     BoundaryStateDim_V HpUx
0.0                     BoundaryStateDim_V HpUy
0.0                     BoundaryStateDim_V HpUz
0.0                     BoundaryStateDim_V Bx
0.0                     BoundaryStateDim_V By
0.0                     BoundaryStateDim_V Bz
0.025                   BoundaryStateDim_V Pe
0.125                   BoundaryStateDim_V HpP
5.0                    BoundaryStateDim_V He2pRho
0.0                     BoundaryStateDim_V He2pUx
0.0                     BoundaryStateDim_V He2pUy
0.0                     BoundaryStateDim_V He2pUz
0.125                   BoundaryStateDim_V He2pP
0.0                     BoundaryStateDim_V Hyp

#SOLIDSTATE
T                       UseSolidState
user                    TypeBcSolid
sphere                  TypeSolidGeometry
1.0                     rSolid
5e-3                    SolidLimitDt

-------------end BC--------------

#USERSWITCH
+init +ic		StringSwitch

#USERINPUTBEGIN --------------------

#PEOVERP
0.2                     ***What does this do?

RESISTIVEPLANET
1.0                     PlanetRadius
5                       nResistivPoints
1.02                    Radius
0.0                     Resistivety
0.95                    Radius
1.0e13                  Resistivety
0.85                    Radius
1.0e13                  Resistivety
0.83                    Radius
1.0e7                   Resistivety
0.8                     Radius
0.0                     Resistivety

#RESISTIVEPLANET
5                      nResistivPoints
1.02                   Radius
0.0                    Resistivety
0.95                   Radius
1.0e13                  Resistivety
0.85                    Radius
1.0e13                  Resistivety
0.83                   Radius
0                      Resistivety
0.8                    Radius
0.0                   Resistivety

RESISTIVEPLANET
1.0                     PlanetRadius
5                       nResistivPoints
1.06                    Radius
0.0                     Resistivety
0.95                    Radius
1.0e13                  Resistivety
0.90                    Radius
1.0e13                  Resistivety
0.88                    Radius
0.0                     Resistivety
0.8                     Radius
0.0                     Resistivety

RESISTIVEPLANET
1.0                     PlanetRadius
5                       nResistivPoints
1.05                    Radius
0.0                     Resistivity
0.95                    Radius
6e11                    Resistivity
0.70                    Radius
6e11                    Resistivity
0.60                    Radius
6e9                     Resistivity
0.55                    Radius
0.0                     Resistivity

#USERINPUTEND ----------------------

MHDIONS
F            DoAddRho
T            DoAddRhoU

#MULTIION
1e-5            LowDensityRatio
1e-11            LowPressureRatio
F            DoRestrictMultiIon

#MULTIIONSTATE
F            UseSingleIonVelocity
T            UseSingleIonTemperature

#SCHEME
1                       nOrder (1 or 2)
Rusanov                   TypeFlux (Roe, Rusanov, Linde, Sokolov

UNUSED
mc3                     TypeLimiter
1.2                     LimiterBeta

#POINTIMPLICIT
T           UsePointImplicit
1.0         BetaPointImplicit (read if UsePointImplicit is true)
F           IsAsymmetric
T           DoNormalizeCell

ELECTRONENTROPY
F                       UseElectronEntropy
T                       UseElectronEnergy

#TIMESTEPPING
1                       nStage
0.6                     CflExlp

#HALLRESISTIVITY
F                       UseHallResist (rest of parameters read only if true)
4.0                     HallFactorMax
1.0                     HallCmaxFactor

#REGION
hallbox1                NameRegion
box tapered             NameHallRegion
-8                      xMinBox
-4                      yMinBox
-4                      zMinBox
2                       xMaxBox
4                       yMaxBox
4                       zMaxBox
0.5                     Taper

#REGION
hallsphere              NameRegion
sphere0 tapered         StringShape
1.15                     Radius
0.05                    Taper

#HALLREGION
+hallbox1 -hallsphere

#SAVELOGFILE
T                       DoSaveLogfile
RAW                     StringLogfile
1                       DnSaveLogfile
-1.                     DtSaveLogfile

#SAVEINITIAL
T

#SAVEPLOT
1                       nPlotFiles
y=0 VAR idl_ascii       StringPlot
100                    DnSavePlot
-1.                     DtSavePlot
-1.                     Dx
{MHD} b1x b1y b1z eta divb dt dx hall    NameVars
{default}                       NamePars
z=0 VAR idl_ascii       StringPlot
10                    DnSavePlot
-1.                     DtSavePlot
-1.                     Dx
{MHD} b1x b1y b1z eta divb dt dx hall    NameVars
{default}                       NamePars
3d VAR tec             plot_strin              StringPlot
5000                    DnSavePlot
-1                      DtSavePlot
rho ux uy uz b1x b1y b1z bx by bz p eta jx jy jz dt dtblk cons impl dx     NameVars
{default} rbody                 NamePars

#END_COMP GM ---------------------------------------

#STOP
1000000                  MaxIteration
-1000                   tSimulationMax

#RUN


#END




#ECHO
T

#DESCRIPTION
MESSENGER M2 flyby run using the new layered inner boundary condition

#PLANET
Ganymede                NamePlanet
2439000                 RadiusPlanet [m]
1.4819E+23              MassPlanet   [kg]
0.0                     OmegaPlanet  [radian/s]
0.0                     TiltRotation [degree]
DIPOLE                  TypeBField
0.0                     MagAxisThetaGeo [degree]
289.1550                MagAxisPhiGeo   [degree]
-200.0E-9            DipoleStrength  [T]

#ROTATION
F               UseRotation

INCLUDE
RESTART.in              NameRestartFile

MAGNETICAXIS
T                       IsMagAxisPrimary (rest of parameters read if true)
180.0
0.0

#MAGNETICCENTER
0.0              MagCenterX
0.0             MagCenterY
0.2             MagCenterZ

#TIMEACCURATE
F               DoTimeAccurate

CPUTIMEMAX
14000.                  CpuTimeMax

#LAYOUT
GM    0     -1   1   1
PC    0	    -1	 1   1	

#CHECKSTOP
T               DoCheckStop
1000            DnCheckStop
-0.015          DtCheckStop


#SAVERESTART
T                   DoSaveRestart
10000               DnSaveRestart
-5                 DtSaveRestart

#CPUTIMEMAX
70000.0         CpuTimeMax

COUPLE2TIGHT
GM                      NameCompMaster
PC                      NameCompSlave
T                       DoCouple

COUPLE2
GM                      NameCompMaster
PC                      NameCompSlave
-1                      DnCouple
0.005                   DtCouple

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

#BODY
T                       UseBody
0.8                     rBody      [rPlanet]
3.5                     rCurrents  [rPlanet]
39                      BodyNDim   [/cc]
2.32e5                  BodyTDim   [K]

! use BodyNDim, BodyTDim only, rbody as a parameter for plot
#BODY
F

SOLARWIND
35.8852                 SwNDim   [/cc]
8.70e4                  SwTDim   [K] (3.8nPa/(1.38e-23*2/cc))/1.2
-500.0                   SwUxDim  [km/s]
0.0                     SwUyDim  [km/s]
0.0                     SwUzDim  [km/s]
0.0                     SwBxDim  [nT]
0.0                     SwByDim  [nT]
-22.88                  SwBzDim  [nT]

#PLASMA
1.0                     FluidMass [amu]
1.0                     AverageIonCharge [e]
1.0                 ElectronTemperatureRatio

TEST
krylov

KRYLOV
GMRES                   TypeKrylov  (GMRES, BICGSTAB, CG)
nul                     TypeInitKrylov (nul, old, explicit, scaled)
0.001                   ErrorMaxKrylov
200                     MaxMatvecKrylov

#MINIMUMPRESSURE
0.00001                   pMinDim
0.00001                   PeMinDim for electron pressure

#MINIMUMDENSITY
0.001                     RhoMinDim

#NONCONSERVATIVE
F                       UseNonConservative

#CONSERVATIVECRITERIA
0                       nConservCrit

#RESTARTOUTFILE
one                     TypeRestartOutFile

! Grid structure info
#INCLUDE
Grid_init

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
coord1min solid	                StringBoundary
5.0                     BoundaryStateDim_V Rho
0.0                     BoundaryStateDim_V Ux
0.0                     BoundaryStateDim_V Uy
0.0                     BoundaryStateDim_V Uz
0.0                     BoundaryStateDim_V Bx
0.0                     BoundaryStateDim_V By
0.0                     BoundaryStateDim_V Bz
0.0                     BoundaryStateDim_V Hyp
0.025                   BoundaryStateDim_V Pe
0.125                   BoundaryStateDim_V p

#BOUNDARYSTATE
xmaxbox                 StringBoundary
36.0                     BoundaryStateDim_V Rho
-500.0                     BoundaryStateDim_V Ux
0.0                     BoundaryStateDim_V Uy
0.0                     BoundaryStateDim_V Uz
0.0                     BoundaryStateDim_V Bx
23.0                     BoundaryStateDim_V By
0.0                     BoundaryStateDim_V Bz
0.0                     BoundaryStateDim_V Hyp
0.043                   BoundaryStateDim_V Pe
0.043                   BoundaryStateDim_V p

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
0.2

RESISTIVEPLANET
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

#USERINPUTEND ----------------------

#ELECTRONENTROPY
F                       UseElectronEntropy
T                       UseElectronEnergy

#TIMESTEPPING
1                       nStage
0.8                     CflExlp

#SCHEME
1                       nOrder (1 or 2)
Rusanov                   TypeFlux (Roe, Rusanov, Linde, Sokolov
mc3                     TypeLimiter
1.2                     LimiterBeta


#UNIFORMAXIS
T                   UseUniformAxis

#COARSEAXIS
T                       UseCoarseAxis
1                       nCoarseLayer

#HALLRESISTIVITY
F                       UseHallResist (rest of parameters read only if true)
4.0                     HallFactorMax
1.0                     HallCmaxFactor

#REGION
hallbox1                NameRegion
box tapered             NameHallRegion
-8.0                      xMinBox
-4.0                      yMinBox
-4.0                      zMinBox
2                       xMaxBox
-4.0                       yMaxBox
4.0                       zMaxBox
0.5                     Taper

#REGION
hallsphere              NameRegion
sphere0 tapered         StringShape
1.05                     Radius
0.05                    Taper

#HALLREGION
+hallbox1 -hallsphere

#SAVELOGFILE
T                       DoSaveLogfile
RAW                     StringLogfile
1                       DnSaveLogfile
-1.                     DtSaveLogfile

#SAVEPLOT
2                       nPlotFiles
y=0 VAR idl_ascii       StringPlot
5000                    DnSavePlot
-1.                     DtSavePlot
-1.                     Dx
{MHD} b1x b1y b1z eta divb dt dx hall    NameVars
{default}                       NamePars
z=0 VAR idl_ascii       StringPlot
5000                    DnSavePlot
-1.                     DtSavePlot
-1.                     Dx
{MHD} b1x b1y b1z eta divb dt dx hall    NameVars
{default}                       NamePars
3d VAR tec             plot_strin              StringPlot
5000                    DnSavePlot
-5                      DtSavePlot
rho ux uy uz b1x b1y b1z bx by bz p eta jx jy jz dt dtblk cons impl dx     NameVars
{default} rbody                 NamePars

#END_COMP GM ---------------------------------------

#STOP
30000                  MaxIteration
-1000                   tSimulationMax

#RUN

#BEGIN_COMP GM ---------------------------------------
#TIMESTEPPING
2                       nStage
0.8                     CflExlp

#SCHEME
2                       nOrder (1 or 2)
Rusanov                   TypeFlux (Roe, Rusanov, Linde, Sokolov
mc3                     TypeLimiter
1.2                     LimiterBeta

#END_COMP GM ---------------------------------------

#STOP
50000                  MaxIteration
-1000                   tSimulationMax



#RUN

#BEGIN_COMP GM ---------------------------------------

#AMR
100
F


#GRIDLEVEL
3               nLevelArea
box_gen         East magnetopause
1.2             rmin            xMinBox
0.0             LonMin          yMinBox
-70.0           LatMin          zMinBox
1.8            rmax            xMaxBox
90.0           LonMax          yMaxBox
70.0            LatMax          zMaxBox

#GRIDLEVEL
3               nLevelArea
box_gen         West magnetopause
1.2             rmin            xMinBox
270.0             LonMin          yMinBox
-70.0           LatMin          zMinBox
1.8            rmax            xMaxBox
360.0           LonMax          yMaxBox
70.0            LatMax          zMaxBox

#GRIDLEVEL
3               nLevelArea
box             tail box
-5            xMinBox
-1.5            yMinBox
-0.6            zMinBox
-1.1            xMaxBox
1.5             yMaxBox
1.1             zMaxBox

GRIDLEVEL
4               nLevelArea
box             tail box fine
-4.0            xMinBox
-1.2            yMinBox
-0.4            zMinBox
-1.15            xMaxBox
1.2             yMaxBox
0.8            zMaxBox

GRIDLEVEL
4               nLevelArea
box_gen             tail box fill
1.15             rmin           * Corrected from run0/ss-15 EXPERIMENTAL
90.0             LonMin          yMinBox
-20.0           LatMin          * Corrected from run0/ss-15
1.6             rmax            * Corrected from run0/ss-15
270.0           LonMax          yMaxBox
35.0            LatMax          zMaxBox

#HYPERBOLICDIVB
T                       UseHyperbolicDivb
250.0                   SpeedHypDim
0.2                     HypDecay

#END_COMP GM ---------------------------------------

#STOP
55000                  MaxIteration
-1000                   tSimulationMax

#RUN

#BEGIN_COMP GM ---------------------------------------

#AMR
-100
F

#TEST
krylov

#RESISTIVITY
T                       UseResistivity
user                    TypeResistivity
0.0                     Eta0Si

#SEMIIMPLICIT
T                       UseSemiImplicit
resistivity             TypeSemiImplicit

#END_COMP GM ---------------------------------------

#STOP
100000                  MaxIteration
-1000                   tSimulationMax

#RUN

#END



#!/bin/ksh
##     ***  script  config_monarc.ksh  ***
##  Set environment variables for DMONTOOLS/ARC-DIAGS scripts
## =====================================================================
## History : 1.0  !  2008     J.M. Molines    Original  code
##           2.0  !  2012     all contrib     Rationalization
##            -   !  2018     C. Talandier    Original version for Arctic
## ----------------------------------------------------------------------
#####################################################################################
#####################################################################################

CONFIG=<CONFIG>               # The configuration name
CASE=<CASE>                   # The reference experiment
CASE2=                        # The experiment to compare with the reference CASE one
REF_YEAR=                     # The starting year of the experiment for time-series
XIOS=1m   	              # [1m/5d] depending the output frequency directory
GETFILE=1                     # [0/1] to set if input files required both obs. and model 
OUTNCDF=1                     # [0/1] to output diagnostics into NetCDF 

# Specify input files name (should be stored in the GRID_DIR directory to set below)
INIT_FILE=                    # Initial file with both temp/sal initial state store in GRID_DIR directory
BATH_FILE=                    # Bathymetry file store in GRID_DIR directory
GRDH_FILE=                    # Horizontal mesh grid file
GRDZ_FILE=                    # Vertical mesh grid file
MASK_FILE=                    # Mask grid file 
COOR_FILE=                    # Coordinates grid file 

#####################################################################################
#####################################################################################
#####################################################################################
# >>>>>>>>>      FIRST ORDER SELECTION OF PLOTS & DIAGS
########################################################
########### MAPS ###############
# To perform horizontal maps and diagnostics as FWC 
# --------------------------------------------------
DOMAPS=1
########## SECTIONS ############
# To perform section plots of Temp/Sal & velocity as quantitative calculation 
# ---------------------------------------------------------------------------
DOSECTIONS=1
########## MOORINGS ############
# To perform Tem/sal annual mean or/and time evolution profiles 
# ---------------------------------------------------------------------------------------
DOMOORINGS=1
########## INTEGRATED QUANTITIES ############
# To perform the calculation of sepcific diagnostics as integrated quantities over subdomain
# ---------------------------------------------------------------------------------------
DOINTQUANT=1

#####################################################################################
# >>>>>>>>>      SECOND ORDER SELECTION OF PLOTS & DIAGS
########################################################
########### MAPS ###############
################################
# List of maps to plot
# AWT keyword for temperature max & depth
# FWC keyword for freshwater content & ssh
# ICE keyword for ice thickness & march/september sea-ice concentration
# MLD keyword for mixed layer depth in march & september
# DYN keyword for surface & 100m depth EKE & stream function
# TSD keyword for surface & 100m depth temperature & salinity differences with initial state
# MTS keyword for mean temperature & salinity in the ML in March and September

MAPS='AWT FWC ICE MLD DYN MTS'

# Compute a mean over a period spanning given years in argument
MKCLIM=0

# To plot (to do only once) a map that summarizes sections, moorings and boxes location
PLT_MAP=1

#####################################################################################
################################
########## SECTIONS ############
################################
# List of sections for which time series of net volume, heat, salt & Ice can be computed
MKDIAGSSEC=[FramS,Davis,Bering]

# List of sections for which salinity, temperature & cross section velocity are plotted
MKSECPLOT=[FramS,Beauf,Kara]

#####################################################################################
################################
########## MOORINGS ############
################################
# Select the mooring [which is currently a grid point NOT a BOX]
# Two moorings by default: EURA within the Eurasian basin & ARCB within the Beaufort Gyre
BOX="['ARCB','EURA']"

#####################################################################################
#############################################
########## INTEGRATED QUANTITIES ############
#############################################
# To compute diags such as FWC within the CRF box centered over the Beaufort Gyre

#####################################################################################
MACHINE=<MACHINE>
INITDIR=<INITDIR>

# set xiosid according to XIOS FLAG
if [ $XIOS ] ; then  xiosid=.$XIOS ; else  xiosid= ; fi 


case $MACHINE in
    ( datarmor )

    WORKDIR=$SCRATCH
    USER=`whoami`
    REMOTE_USER=`whoami`
    DATA_DIR=/home/datawork-lops-drakkarcom/SIMULATION-OUTPUTS/FREDY/CONFIGS/${CONFIG}/${CONFIG}-${CASE}-MEAN/       # input directory (main experiment)
    GRID_DIR=$WORKDIR/${CONFIG}/${CONFIG}-I/                  					  # grid directory
    DATA_DIR2=/home/datawork-lops-drakkarcom/SIMULATION-OUTPUTS/FREDY/CONFIGS/${CONFIG}/${CONFIG}-${CASE2}-MEAN/     # input directory (2nd experiment)
    WPDIR=$WORKDIR/WRUN_${CONFIG}/${CONFIG}-${CASE}/CTL/CDF   					  # Working Pdir for templates scripts
    DMONTOOLS=$HOME/DEV/DMONTOOLS                             					  # DMONTOOLS location
    OBS_DIR=/home/datawork-lops-drakkarcom/DATA-REFERENCE/FOR-MONITORING/    			  # observations directory 
    BATCH=qsub 

    ;;

    ( occigen )

    USER=`whoami`
    REMOTE_USER=`whoami`
    DATA_DIR=$WORKDIR/${CONFIG}/${CONFIG}-${CASE}-MEAN/      # input directory (main experiment)
    DATA_DIR2=$WORKDIR/${CONFIG}/${CONFIG}-${CASE2}-MEAN/    # input directory (2nd experiment)
    GRID_DIR=$SCRATCHDIR/${CONFIG}/${CONFIG}-I/              # grid directory
    WPDIR=$WORKDIR/WRUN_${CONFIG}/${CONFIG}-${CASE}/CTL/CDF  # Working Pdir for templates scripts
    DMONTOOLS=$SCRATCHDIR/DEV/DMONTOOLS_2.1                  # DMONTOOLS location
    OBS_DIR=$SHAREDSCRATCHDIR/DATA_OBS/                      # observations directory 
    BATCH=sbatch 

    ;;

    *)
    echo available machine is occigen  ; exit 1 ;;

esac

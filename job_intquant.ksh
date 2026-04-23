
set -xv
ulimit -s
ulimit -s unlimited


CONFIG=XXCONFXX
CASE=XXCASEXX
CASE2=XXCASE2XX
CONFCASE=${CONFIG}-${CASE}   ; CONFCASE2=${CONFIG}-${CASE2}
DATA_DIR="XXDATA_DIRXX"
GRID_DIR="XXGRID_DIRXX"
DATA_DIR2="XXDATA_DIRXX"
XIOS="XXXIOSFREQXX"

S_Y=XXSYEAXX 
E_Y=XXEYEAXX



. ../config_monarc.ksh

chkdir() { if [ ! -d $1 ] ; then mkdir -p $1 ; else echo $1 exists ; fi ; }

OUTFIGS=${INITDIR}/INTQUANT/FIGS-${CONFIG}/${CONFIG}-${CASE}
OUTNCDF=${INITDIR}/INTQUANT/NCDF-${CONFIG}/${CONFIG}-${CASE}
chkdir ${OUTFIGS}
chkdir ${OUTNCDF}

cd ${WPDIR}/MONARC/INTQUANT

cp ${INITDIR}/INTQUANT/CREG_intquant_func.py .
cp ${DMONTOOLS}/checkfile.py .


getfiles=XXGETFILEXX
if [ $getfiles -eq 1 ] ; then 
	# main program
	# Get all necessary data input files
        chkdir ${CONFIG}/${CONFCASE}-MEAN/DATA
	cd ${CONFIG}/${CONFCASE}-MEAN
	cy=$S_Y
	# loop over years
	while [ $cy -le ${E_Y} ]  ; do 
	       chkdir $cy      ;  cd $cy 
	       chkdir ${XIOS}  ;  cd ${XIOS}
	       # Monthly files 
	       file_base=${CONFCASE}_y${cy}m*
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_gridT.nc .
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_icemod.nc .
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${CONFCASE}_y${cy}.${XIOS}_SurOceCurl.nc .
               cd ../
               chkdir 1m   ;  cd 1m
	       # Yearly files 
	       file_base=${CONFCASE}_y${cy}.1m
	       ln -sf  ${DATA_DIR}/1m/$cy/${file_base}_SBC_scalar.nc
	       ln -sf  ${DATA_DIR}/1m/$cy/${file_base}_OCE_scalar.nc
               
	       let cy=$cy+1
	       cd ../..
	done
	
	cd ../..


	# Get In-situ data
	chkdir DATA  ;  cd DATA
	ln -sf ${OBS_DIR}/ICE/PIOMAS_icevol_maskedBeringSea_interp${CONFIG}_1-12_1979-2024.nc .
	ln -sf ${OBS_DIR}/ICE/NSIDC-G02202-V4_ice_area_and_extent_TiSe_y1978-11-2022-12_maskBeringSea_fullPoleGap.nc .
        ln -sf ${OBS_DIR}/OCEAN/FRAM_inflow.mat .
        ln -sf ${OBS_DIR}/ICE/IABP_ice_drift_BG_1979-2016.mat .
        ln -sf ${OBS_DIR}/OCEAN/BGFWC_OI_2023.nc .
        ln -sf ${OBS_DIR}/OCEAN/ArcticEkmanPumping_MonthlyMean.nc .
	cd ../

	# link the bathymetry file
	ln -sf ${GRID_DIR}/${BATH_FILE} Bathymetry.nc 

	# Get all necessary grid files
	chkdir ${CONFIG}/GRID
	cd ${CONFIG}/GRID
	if [ ! -f ${CONFCASE}_mesh_hgr.nc ] ; then ln -sf $GRID_DIR/${GRDH_FILE}  ${CONFCASE}_mesh_hgr.nc ; fi 
	if [ ! -f ${CONFCASE}_mesh_zgr.nc ] ; then ln -sf $GRID_DIR/${GRDZ_FILE}  ${CONFCASE}_mesh_zgr.nc ; fi
	if [ ! -f ${CONFCASE}_mask.nc ] ; then ln -sf $GRID_DIR/${MASK_FILE} ${CONFCASE}_mask.nc ; fi
	if [ ! -f ${CONFCASE}_coordinates.nc ] ; then ln -sf $GRID_DIR/${COOR_FILE} ${CONFCASE}_coordinates.nc ; fi

	cd ../../
fi

pwd 

./CREG_intquant_y${S_Y}${E_Y}.py

ls -lrt 

mv ${CONFCASE}*y${S_Y}*.png CREG*LGTS*.png ${OUTFIGS}/.

mv ./NETCDF/${CONFCASE}*.nc ${OUTNCDF}/.




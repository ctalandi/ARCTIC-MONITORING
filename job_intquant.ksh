
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
cp ${DMONTOOLS}/MONARC/checkfile.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf_plots.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf_colormaps.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf.py .


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
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_flxT.nc .
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
	if [ ${CONFIG} == 'CREG025.L75' ] ; then 
		ln -sf ${OBS_DIR}/ICE/PIOMAS_icevol_maskedBeringSea_interpCREG025.L75_1-12_1979-2018.nc .
	elif [ ${CONFIG} == 'CREG12.L75' ] ; then
		ln -sf ${OBS_DIR}/ICE/PIOMAS_icevol_maskedBeringSea_interpCREG12.L75_1-12_1979-2018.nc .
        fi

	ln -sf ${OBS_DIR}/ICE/NSIDC_ice_area_and_extent_maskBeringSea_fullPoleGap.nc .
        ln -sf ${OBS_DIR}/OCEAN/FRAM_inflow.mat .
        ln -sf ${OBS_DIR}/ICE/ice_drift_BG_1979-2011.mat .
        ln -sf ${OBS_DIR}/OCEAN/BeaufortGyreFWC-Obs-Proshutinsky_GRL2018_y2003-2017.nc .

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

mv ${CONFCASE}*y${S_Y}*.pdf CREG*LGTS*.pdf ${OUTFIGS}/.

mv ./NETCDF/${CONFCASE}*.nc ${OUTNCDF}/.




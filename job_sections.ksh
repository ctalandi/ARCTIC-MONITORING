
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

OUTFIGS=${INITDIR}/SECTIONS/FIGS-${CONFIG}/${CONFIG}-${CASE}
OUTNCDF=${INITDIR}/SECTIONS/NCDF-${CONFIG}/${CONFIG}-${CASE}
chkdir ${OUTFIGS}
chkdir ${OUTNCDF}

cd ${WPDIR}/MONARC/SECTIONS

cp ${INITDIR}/SECTIONS/CREG_sections_func.py .

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
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_gridU.nc .
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_gridV.nc .
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_gridW.nc .
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_icemod.nc .
               cd ../
               chkdir 1m  ;  cd 1m
	       # Yearly files 
	       file_base=${CONFCASE}_y${cy}.1m
	       ln -sf  ${DATA_DIR}/1m/$cy/${file_base}_SBC_scalar.nc
	       ln -sf  ${DATA_DIR}/1m/$cy/${file_base}_OCE_scalar.nc
               
	       let cy=$cy+1
	       cd ../..

	       # Link to the initial state file
	       ln -sf  ${GRID_DIR}/${INIT_FILE} ${CONFCASE}_init_gridT.nc
	done
	
	cd ../..

	# link the bathymetry file
	ln -sf ${GRID_DIR}/${BATH_FILE} Bathymetry.nc 
        ln -sf ${OBS_DIR}/OCEAN/FRAM_inflow.mat .

	# Get all necessary grid files
	chkdir ${CONFIG}/GRID
	cd ${CONFIG}/GRID
	if [ ! -f ${CONFCASE}_mesh_hgr.nc ] ; then ln -sf $GRID_DIR/${GRDH_FILE}  ${CONFCASE}_mesh_hgr.nc ; fi 
	if [ ! -f ${CONFCASE}_mesh_zgr.nc ] ; then ln -sf $GRID_DIR/${GRDZ_FILE}  ${CONFCASE}_mesh_zgr.nc ; fi
	if [ ! -f ${CONFCASE}_mask.nc ] ; then ln -sf $GRID_DIR/${MASK_FILE} ${CONFCASE}_mask.nc ; fi
	
	cd ../../
fi

pwd 

./CREG_sections_y${S_Y}${E_Y}.py

ls -lrt 

mv *LGTS*.pdf ${INITDIR}/INTQUANT/FIGS-${CONFIG}/${CONFIG}-${CASE}/.
mv ${CONFCASE}*y${S_Y}*.pdf ${OUTFIGS}/.

mv ./NETCDF/${CONFCASE}*.nc ${OUTNCDF}/.




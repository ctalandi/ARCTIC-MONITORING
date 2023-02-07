
set -xv
ulimit -s
ulimit -s unlimited


CONFIG=XXCONFXX
CASE=XXCASEXX
CASE2=XXCASE2XX
CONFCASE=${CONFIG}-${CASE}   ; CONFCASE2=${CONFIG}-${CASE2}
DATA_DIR="XXDATA_DIRXX"
GRID_DIR="XXGRID_DIRXX"
DATA_DIR2="XXDATA_DIR2XX"
XIOS="XXXIOSFREQXX"

S_Y=XXSYEAXX 
E_Y=XXEYEAXX

. ../config_monarc.ksh

chkdir() { if [ ! -d $1 ] ; then mkdir -p $1 ; else echo $1 exists ; fi ; }

OUTFIGS=${INITDIR}/MOORINGS/FIGS-${CONFIG}/${CONFIG}-${CASE}
OUTNCDF=${INITDIR}/MOORINGS/NCDF-${CONFIG}/${CONFIG}-${CASE}
DATADIR=${WPDIR}/MONARC/MOORINGS/DATA/${CONFIG}
chkdir ${OUTFIGS}
chkdir ${OUTNCDF}
chkdir ${DATADIR}

cd ${WPDIR}/MONARC/MOORINGS
cp ${INITDIR}/MOORINGS/CREG_moorings_func.py .
cp ${INITDIR}/MOORINGS/CREG_moorings_cont.py .
cp ${DMONTOOLS}/MONARC/checkfile.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf_plots.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf_colormaps.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf.py .

getfiles=XXGETFILEXX     
if [ $getfiles -eq 1 ] ; then 
	# main program
	# Get all necessary data input files
	chkdir ${CONFIG}/${CONFCASE}-MEAN
	cd ${CONFIG}/${CONFCASE}-MEAN
	cy=$S_Y
	# loop over years
	while [ $cy -le ${E_Y} ]  ; do 
	       chkdir $cy  ;  cd $cy 
	       # Monthly files 
	       file_base=${CONFCASE}_y${cy}m*
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_gridT.nc .
	       # Yearly files 
	       file_base=${CONFCASE}_y${cy}.${XIOS}
	       ln -sf  ${DATA_DIR}/${XIOS}/$cy/${file_base}_gridT.nc .
	       let cy=$cy+1
	       cd ..
	done

	# Link to the initial state file
	ln -sf  ${GRID_DIR}/${INIT_FILE} ${CONFCASE}_init_gridT.nc
	cd ../../

        chkdir DATA  ;  cd DATA
        ln -sf ${OBS_DIR}/OCEAN/TimeMean_profiles_moorings.mat .
	cd ../
	
	# Get all necessary grid files
	chkdir ${CONFIG}/GRID
	cd ${CONFIG}/GRID
	if [ ! -f ${CONFCASE}_mesh_hgr.nc ] ; then ln -sf $GRID_DIR/${GRDH_FILE}  ${CONFCASE}_mesh_hgr.nc ; fi 
	if [ ! -f ${CONFCASE}_mesh_zgr.nc ] ; then ln -sf $GRID_DIR/${GRDZ_FILE}  ${CONFCASE}_mesh_zgr.nc ; fi
	if [ ! -f ${CONFCASE}_mask.nc ] ; then ln -sf $GRID_DIR/${MASK_FILE} ${CONFCASE}_mask.nc ; fi
	cd ../../

fi

./${CASE}_MOOR_y${S_Y}${E_Y}.py
    
ls -lrt 

mv ${CONFCASE}*${S_Y}*.pdf ${CONFCASE}*LASTy.pdf ${OUTFIGS}/.

mv ./NETCDF/${CONFCASE}*.nc ${OUTNCDF}/.






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

OUTFIGS=${INITDIR}/MAPS/FIGS-${CONFIG}/${CONFIG}-${CASE}
OUTNCDF=${INITDIR}/MAPS/NCDF-${CONFIG}/${CONFIG}-${CASE}
chkdir ${OUTFIGS}
chkdir ${OUTNCDF}

cd ${WPDIR}/MONARC/MAPS

cp ${INITDIR}/MAPS/CREG_maps_cont.py .
cp ${INITDIR}/MAPS/CREG_maps_func.py .
cp ${DMONTOOLS}/MONARC/checkfile.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf_plots.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf_colormaps.py .
cp ${DMONTOOLS}/MONARC/PyRaf/PyRaf.py .

getfiles=XXGETFILEXX
if [ $getfiles -eq 1 ] ; then 
	# main program
	# Get all necessary data input files
        for ZEXP in `echo ${CASE} ${CASE2}`   ;  do 
             ZCONFCASE=${CONFIG}-${ZEXP}
             ZDATA_DIR=${DATA_DIR}
	     if [ ${ZEXP} == ${CASE2} ]  ;then 
                  ZDATA_DIR=${DATA_DIR2}
             fi
	     chkdir ${CONFIG}/${ZCONFCASE}-MEAN
	     cd ${CONFIG}/${ZCONFCASE}-MEAN
	     cy=$S_Y
	     # loop over years
	     while [ $cy -le ${E_Y} ]  ; do 
	            chkdir $cy     ;  cd $cy 
	            chkdir ${XIOS} ;  cd ${XIOS}

	            # Monthly files 
	            file_base=${ZCONFCASE}_y${cy}m*
	            ln -sf  ${ZDATA_DIR}/${XIOS}/$cy/${file_base}_gridT.nc .
	            ln -sf  ${ZDATA_DIR}/${XIOS}/$cy/${file_base}_icemod.nc .

	            # Yearly files 
	            file_base=${ZCONFCASE}_y${cy}.${XIOS}
	            ln -sf  ${ZDATA_DIR}/${XIOS}/$cy/${file_base}_gridT.nc .
	            ln -sf  ${ZDATA_DIR}/${XIOS}/$cy/${file_base}_icemod.nc .
	            ln -sf  ${ZDATA_DIR}/${XIOS}/$cy/${file_base}_PSI.nc .
	            ln -sf  ${ZDATA_DIR}/${XIOS}/$cy/${file_base}_EKE.nc .
                    cd ../
                    chkdir 1m  ;  cd 1m
	            # Yearly files 
	            file_base=${ZCONFCASE}_y${cy}.1m
	            ln -sf  ${ZDATA_DIR}/1m/$cy/${file_base}_SBC_scalar.nc
                    
	            let cy=$cy+1
	            cd ../..
	     done

	     # Link to the initial state file
	     ln -sf  ${GRID_DIR}/${INIT_FILE} ${CONFCASE}_init_gridT.nc
	
	     cd ../..
	done
	
	# link to the observations files
        chkdir DATA  ;  cd DATA
        ln -sf ${OBS_DIR}/OCEAN/phc3.0_annual.nc .
        ln -sf ${OBS_DIR}/OCEAN/EKE_DOT_based_2003-2014.nc .
        ln -sf ${OBS_DIR}/OCEAN/MLD_MIMOC_based_monthlyClim_rhocrit0.01.nc .
        ln -sf ${OBS_DIR}/OCEAN/MIMOC_ML_v2.2_PT_S_MLP_Clim.nc .
        ln -sf ${OBS_DIR}/OCEAN/BeaufortGyreFWC-Obs-Proshutinsky_GRL2018_y2003-2017.nc .

        ln -sf ${OBS_DIR}/ICE/NSIDC-0051_92585_monthly.nc
	if  [ ${CONFIG} == 'CREG025.L75'   ]  ; then 
        	ln -sf ${OBS_DIR}/ICE/PIOMAS_icethic_interpCREG025.L75_1-12_1979-2018.nc
	elif [ ${CONFIG} == 'CREG12.L75'   ]  ; then
        	ln -sf ${OBS_DIR}/ICE/PIOMAS_icethic_interpCREG12.L75_1-12_1979-2018.nc
	fi
        cd ../

	# link the bathymetry file
	ln -sf ${GRID_DIR}/${BATH_FILE} Bathymetry.nc 

	# Get all necessary grid files
	chkdir ${CONFIG}/GRID
	cd ${CONFIG}/GRID
	if [ ! -f ${CONFCASE}_mesh_hgr.nc ] ; then ln -sf $GRID_DIR/${GRDH_FILE}  ${CONFCASE}_mesh_hgr.nc ; fi 
	if [ ! -f ${CONFCASE}_mesh_zgr.nc ] ; then ln -sf $GRID_DIR/${GRDZ_FILE}  ${CONFCASE}_mesh_zgr.nc ; fi
	if [ ! -f ${CONFCASE}_mask.nc ] ; then ln -sf $GRID_DIR/${MASK_FILE} ${CONFCASE}_mask.nc ; fi
	if [ ${CONFIG} == 'CREG12.L75' ] ; then ln -sf ${SCRATCHDIR}/CREG025.L75/GRID/CREG025.L75_byte_mask.nc CREG025.L75_mask.nc ; fi

	cd ../../
fi

pwd 

if [ $MACHINE == 'occigen' ] ;   then
   module purge
   module list

   export GEOS_DIR=/opt/software/occigen/libraries/geos/3.6.2
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GEOS_DIR/lib
   module load intel/17.1
   module load python/2.7.12
fi

if [ ${PLT_MAP} -eq 1 ]  ;  then
   ./CREG_maps_geog.py 
   mv MONARC_ARC-GEOLOC.pdf ${INITDIR}/.
fi

for ZMAPS in `echo ${MAPS} ` ; do
    ./CREG_maps_y${S_Y}${E_Y}_${ZMAPS}.py 
    mv ${CONFCASE}*${ZMAPS}*${S_Y}*.pdf ${CONFCASE}*${ZMAPS}*${S_Y}*.png ${OUTFIGS}/.
done

mv ./NETCDF/${CONFCASE}*.nc ${OUTNCDF}/.

ls -lrt 





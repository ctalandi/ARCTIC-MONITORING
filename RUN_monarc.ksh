#!/bin/bash
# RUN_monarc.ksh version for SLURM on OCCIGEN
if [ $# = 0 ] ; then
   echo USAGE: RUN_monarc.ksh  syear eyear
   exit 1
fi
#set -x

# get configuration settings
. ./config_monarc.ksh

# set additional variables WPDIR set in config_moy
mkdir -p ${WPDIR}/MONARC/MAPS
mkdir -p ${WPDIR}/MONARC/MOORINGS
mkdir -p ${WPDIR}/MONARC/SECTIONS
mkdir -p ${WPDIR}/MONARC/INTQUANT

# copy required scripts there:
cp -f config_monarc.ksh ${WPDIR}/MONARC/

# get year to work with
STA_YEAR=$1
END_YEAR=$2

#---------------------------------------------------------------------------------------------------------------------------------
#--  MAPS  -----------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------
if [ ${DOMAPS} -eq 1 ]  ; then 

	cd ./MAPS/

	DEPEND="#"
	C_YEAR=${STA_YEAR}
	while [ ${C_YEAR} -le ${END_YEAR}  ]  ;  do
	
	      # Build Python script for each plot
	      for ZMAPS in `echo ${MAPS}`  ; do
                       ZAW_TMAX=False   ; ZFWC_MAPS=False ; ZICE_MAPS=False  ; ZMLD_MAPS=False  ; ZDYN_MAPS=False ; ZTSD_MAPS=False ; ZATL_MAPS=False ; ZMOC_MAPS=False  ; ZNCDF=False  ;   ZTSM_MAPS=False
	                if [ ${ZMAPS} = 'AWT' ] ; then  
                              ZAW_TMAX=True
	                elif [ ${ZMAPS} = 'FWC' ]   ; then 
                              ZFWC_MAPS=True 
	                elif [ ${ZMAPS} = 'ICE' ]   ; then 
                              ZICE_MAPS=True  
	                elif [ ${ZMAPS} = 'MLD' ]   ; then 
                              ZMLD_MAPS=True 
	                elif [ ${ZMAPS} = 'DYN' ]   ; then 
                              ZDYN_MAPS=True
	                elif [ ${ZMAPS} = 'TSD' ]   ; then 
                              ZTSD_MAPS=True
	                elif [ ${ZMAPS} = 'ATL' ]   ; then 
                              ZATL_MAPS=True
	                elif [ ${ZMAPS} = 'MOC' ]   ; then 
                              ZMOC_MAPS=True
	                elif [ ${ZMAPS} = 'MTS' ]   ; then 
                              ZTSM_MAPS=True
                        fi

	                if [ ${OUTNCDF} = 1 ] ; then  
                              ZNCDF=True
			fi

	              sed -e "s/XXSYEAXX/${C_YEAR}/" \
	                  -e "s/XXEYEAXX/${C_YEAR}/" \
	                  -e "s/XXCONFXX/${CONFIG}/" \
	                  -e "s/XXCASEXX/${CASE}/" \
	                  -e "s/XXCASE2XX/${CASE2}/" \
	                  -e "s/XXGRIDXX/${GTYPE}/" \
	                  -e "s;XXDATA_DIRXX;${DATA_DIR};" \
	                  -e "s;XXGRID_DIRXX;${GRID_DIR};" \
	                  -e "s;XXDATA_DIR2XX;${DATA_DIR2};" \
	                  -e "s;XXOBS_DIRXX;${OBS_DIR};" \
                          -e "s/XXXIOSFREQXX/'${XIOS}'/" \
	                  -e "s/XXAW_TMAXXX/${ZAW_TMAX}/" \
	                  -e "s/XXFWC_MAPSXX/${ZFWC_MAPS}/" \
	                  -e "s/XXICE_MAPSXX/${ZICE_MAPS}/" \
	                  -e "s/XXMLD_MAPSXX/${ZMLD_MAPS}/" \
	                  -e "s/XXDYN_MAPSXX/${ZDYN_MAPS}/" \
                          -e "s/XXTSD_MAPSXX/${ZTSD_MAPS}/" \
                          -e "s/XXATL_MAPSXX/${ZATL_MAPS}/" \
                          -e "s/XXMOC_MAPSXX/${ZMOC_MAPS}/" \
                          -e "s/XXMTS_MAPSXX/${ZTSM_MAPS}/" \
                          -e "s/XXNCDFOUTXX/${ZNCDF}/" CREG_maps.py > ${WPDIR}/MONARC/MAPS/CREG_maps_y${C_YEAR}${C_YEAR}_${ZMAPS}.py
	              chmod 750 ${WPDIR}/MONARC/MAPS/CREG_maps_y${C_YEAR}${C_YEAR}_${ZMAPS}.py
	      done

	      if [ ! -d ${WPDIR}/MONARC/MAPS/NETCDF ] ; then mkdir ${WPDIR}/MONARC/MAPS/NETCDF  ; fi
	
	      # Build the job which will be launched (one per grid)
	      sed -e "s/XXSYEAXX/${C_YEAR}/" \
	          -e "s/XXEYEAXX/${C_YEAR}/" \
	          -e "s/XXCONFXX/${CONFIG}/" \
	          -e "s/XXCASEXX/${CASE}/" \
	          -e "s/XXCASE2XX/${CASE2}/" \
	          -e "s;XXDATA_DIRXX;${DATA_DIR};" \
	          -e "s;XXGRID_DIRXX;${GRID_DIR};" \
	          -e "s;XXDATA_DIR2XX;${DATA_DIR2};" \
	          -e "s;XXOBS_DIRXX;${OBS_DIR};" \
	          -e "s/XXGETFILEXX/${GETFILE}/" \
                  -e "s/XXXIOSFREQXX/'${XIOS}'/" \
	          -e "s/XXMAPSXX/${MAPS}/"   job_maps.ksh > ${WPDIR}/MONARC/MAPS/tmp_job_maps

	      sed -e "s/XXDEPENDENCYXX/${DEPEND}/" \
	          -e "s/XXPTYPEXX/MAPS/" \
	          -e "s/XXSPTYPEXX/MAPS/" \
                  -e "s/XXSYEAXX/${C_YEAR}/" \
                  -e "s/XXEYEAXX/${C_YEAR}/" \
                  -e "s/XXCASEXX/${CASE}/" \
	          -e "s;XXINITDIRXX;${INITDIR};" \
	          -e "s/XXTIMEXX/45/" ../JHEAD_${MACHINE}.bash > ${WPDIR}/MONARC/MAPS/tmp_head_maps

              cat ${WPDIR}/MONARC/MAPS/tmp_head_maps ${WPDIR}/MONARC/MAPS/tmp_job_maps > ${WPDIR}/MONARC/MAPS/tmp_job_${CASE}_maps_y${C_YEAR}${C_YEAR}
	
	      ${BATCH} ${WPDIR}/MONARC/MAPS/tmp_job_${CASE}_maps_y${C_YEAR}${C_YEAR}
	
	
	   let C_YEAR=$C_YEAR+1
	
	done

	if [ ${PLT_MAP} -eq 1 ]  ;  then
		sed -e "s/XXCONFXX/${CONFIG}/" \
                    -e "s/XXCASEXX/${CASE}/" CREG_maps_geog.py > ${WPDIR}/MONARC/MAPS/CREG_maps_geog.py
                chmod 750 ${WPDIR}/MONARC/MAPS/CREG_maps_geog.py
	fi

	# Perform a climatological mean
	if [ $MKCLIM -eq 1 ] && [ ${STA_YEAR} != ${END_YEAR} ]   ;  then 
		
		# Build Python script for each plot
		for ZMAPS in `echo ${MAPS}`  ; do
        	         ZAW_TMAX=False   ; ZFWC_MAPS=False ; ZICE_MAPS=False  ; ZMLD_MAPS=False  ; ZDYN_MAPS=False    ; ZTSD_MAPS=False ; ZATL_MAPS=False ; ZMOC_MAPS=False
		          if [ ${ZMAPS} = 'AWT' ] ; then  
        	                ZAW_TMAX=True
		          elif [ ${ZMAPS} = 'FWC' ]   ; then 
        	                ZFWC_MAPS=True 
		          elif [ ${ZMAPS} = 'ICE' ]   ; then 
        	                ZICE_MAPS=True  
		          elif [ ${ZMAPS} = 'MLD' ]   ; then 
        	                ZMLD_MAPS=True 
		          elif [ ${ZMAPS} = 'DYN' ]   ; then 
        	                ZDYN_MAPS=True
	                  elif [ ${ZMAPS} = 'TSD' ]   ; then 
                                ZTSD_MAPS=True
		          elif [ ${ZMAPS} = 'ATL' ]   ; then 
        	                ZATL_MAPS=True
	                elif [ ${ZMAPS} = 'MOC' ]   ; then 
                                ZMOC_MAPS=True
        	          fi

		        sed -e "s/XXSYEAXX/${STA_YEAR}/" \
		            -e "s/XXEYEAXX/${END_YEAR}/" \
		            -e "s/XXCONFXX/${CONFIG}/" \
		            -e "s/XXCASEXX/${CASE}/" \
		            -e "s/XXCASE2XX/${CASE2}/" \
		            -e "s/XXGRIDXX/${GTYPE}/" \
		            -e "s;XXDATA_DIRXX;${DATA_DIR};" \
		            -e "s;XXGRID_DIRXX;${GRID_DIR};" \
		            -e "s;XXDATA_DIR2XX;${DATA_DIR2};" \
		            -e "s;XXOBS_DIRXX;${OBS_DIR};" \
        	            -e "s/XXXIOSFREQXX/'${XIOS}'/" \
		            -e "s/XXAW_TMAXXX/${ZAW_TMAX}/" \
		            -e "s/XXFWC_MAPSXX/${ZFWC_MAPS}/" \
		            -e "s/XXICE_MAPSXX/${ZICE_MAPS}/" \
		            -e "s/XXMLD_MAPSXX/${ZMLD_MAPS}/" \
		            -e "s/XXDYN_MAPSXX/${ZDYN_MAPS}/" \
        		    -e "s/XXTSD_MAPSXX/${ZTSD_MAPS}/" \
		            -e "s/XXATL_MAPSXX/${ZATL_MAPS}/" \
                            -e "s/XXMOC_MAPSXX/${ZMOC_MAPS}/" CREG_maps.py > ${WPDIR}/MONARC/MAPS/CREG_maps_y${STA_YEAR}${END_YEAR}_${ZMAPS}.py
		        chmod 750 ${WPDIR}/MONARC/MAPS/CREG_maps_y${STA_YEAR}${END_YEAR}_${ZMAPS}.py
		done
		
		# Build the job which will be launched (one per grid)
		sed -e "s/XXSYEAXX/${STA_YEAR}/" \
		    -e "s/XXEYEAXX/${END_YEAR}/" \
		    -e "s/XXCONFXX/${CONFIG}/" \
		    -e "s/XXCASEXX/${CASE}/" \
		    -e "s/XXCASE2XX/${CASE2}/" \
		    -e "s;XXDATA_DIRXX;${DATA_DIR};" \
		    -e "s;XXGRID_DIRXX;${GRID_DIR};" \
		    -e "s;XXDATA_DIR2XX;${DATA_DIR2};" \
		    -e "s;XXOBS_DIRXX;${OBS_DIR};" \
		    -e "s/XXGETFILEXX/${GETFILE}/" \
        	    -e "s/XXXIOSFREQXX/'${XIOS}'/" \
		    -e "s/XXMAPSXX/${MAPS}/"   job_maps.ksh > ${WPDIR}/MONARC/MAPS/tmp_job_maps

		sed -e "s/XXDEPENDENCYXX/${DEPEND}/" \
		    -e "s/XXPTYPEXX/MAPS/" \
		    -e "s/XXSPTYPEXX/MAPS/" \
        	    -e "s/XXSYEAXX/${STA_YEAR}/" \
        	    -e "s/XXEYEAXX/${END_YEAR}/" \
        	    -e "s/XXCASEXX/${CASE}/" \
		    -e "s;XXINITDIRXX;${INITDIR};" \
		    -e "s/XXTIMEXX/45/" ../JHEAD_${MACHINE}.bash > ${WPDIR}/MONARC/MAPS/tmp_head_maps

        	cat ${WPDIR}/MONARC/MAPS/tmp_head_maps ${WPDIR}/MONARC/MAPS/tmp_job_maps > ${WPDIR}/MONARC/MAPS/tmp_job_${CASE}_maps_y${STA_YEAR}${END_YEAR}
		
		${BATCH} ${WPDIR}/MONARC/MAPS/tmp_job_${CASE}_maps_y${STA_YEAR}${END_YEAR}
	fi
	
	cd ../
fi

#---------------------------------------------------------------------------------------------------------------------------------
#--  MOORINGS  -------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------
if [ ${DOMOORINGS} -eq 1 ]  ; then 

	cd ./MOORINGS/
	
	LGTS_S=${REF_YEAR}
	LGTS_E=${REF_YEAR}
	C_YEAR=${STA_YEAR}
	DEPEND="#"
	iyear=1
        ZNCDF=False
	if [ ${OUTNCDF} = 1 ] ; then  
              ZNCDF=True
	fi

	while [ ${C_YEAR} -le ${END_YEAR}  ]  ;  do
	
	        if [ ${C_YEAR} -eq ${END_YEAR}  ] ; then LGTS_E=${END_YEAR} ; fi

	        sed -e "s/XXSYEAXX/${C_YEAR}/" \
	            -e "s/XXEYEAXX/${C_YEAR}/" \
	            -e "s/XXCONFXX/${CONFIG}/" \
	            -e "s/XXCASEXX/${CASE}/" \
	            -e "s/XXCASE2XX/${CASE2}/" \
                    -e "s/XXLGTSSXX/${LGTS_S}/" \
                    -e "s/XXLGTSEXX/${LGTS_E}/" \
                    -e "s/XXXIOSFREQXX/'${XIOS}'/" \
	            -e "s/XXALLBOXXX/${BOX}/" \
	            -e "s/XXVARXX/votemper/"  \
	            -e "s/XXNCDFOUTXX/${ZNCDF}/" CREG_moorings.py > ${WPDIR}/MONARC/MOORINGS/${CASE}_MOOR_y${C_YEAR}${C_YEAR}.py
	        chmod 750 ${WPDIR}/MONARC/MOORINGS/${CASE}_MOOR_y${C_YEAR}${C_YEAR}.py
	
	        if [ ${C_YEAR} -eq ${END_YEAR}  ] && [ $(( $END_YEAR-$STA_YEAR )) -gt 0 ] ; then 
		     if [ ${MACHINE} == 'occigen' ]  ; then  DEPEND="#SBATCH --dependency=afterany:${jobids}" ; fi
		     if [ ${MACHINE} == 'datarmor' ] ; then  DEPEND="#PBS -W depend=afterany:${jobids}"       ; fi
		fi

	        if [ ! -d ${WPDIR}/MONARC/MOORINGS/NETCDF ] ; then mkdir ${WPDIR}/MONARC/MOORINGS/NETCDF ; fi

	        sed -e "s/XXSYEAXX/${C_YEAR}/" \
	            -e "s/XXEYEAXX/${C_YEAR}/" \
	            -e "s/XXCONFXX/${CONFIG}/" \
	            -e "s/XXCASEXX/${CASE}/" \
	            -e "s/XXCASE2XX/${CASE2}/" \
	            -e "s;XXDATA_DIRXX;${DATA_DIR};" \
	            -e "s;XXGRID_DIRXX;${GRID_DIR};" \
	            -e "s;XXDATA_DIR2XX;${DATA_DIR2};" \
	            -e "s;XXINITDIRXX;${INITDIR};" \
	            -e "s/XXGETFILEXX/${GETFILE}/" \
                    -e "s/XXXIOSFREQXX/'${XIOS}'/" job_moorings.ksh > ${WPDIR}/MONARC/MOORINGS/tmp_job_moor

	        sed -e "s/XXDEPENDENCYXX/${DEPEND}/" \
	            -e "s/XXPTYPEXX/MOORINGS/" \
	            -e "s/XXSPTYPEXX/MOOR/" \
                    -e "s/XXSYEAXX/${C_YEAR}/" \
                    -e "s/XXEYEAXX/${C_YEAR}/" \
                    -e "s/XXCASEXX/${CASE}/" \
	            -e "s;XXINITDIRXX;${INITDIR};" \
	            -e "s/XXTIMEXX/55/" ../JHEAD_${MACHINE}.bash > ${WPDIR}/MONARC/MOORINGS/tmp_head_moor

                cat ${WPDIR}/MONARC/MOORINGS/tmp_head_moor ${WPDIR}/MONARC/MOORINGS/tmp_job_moor > ${WPDIR}/MONARC/MOORINGS/job_${CASE}_MOOR_y${C_YEAR}${C_YEAR}.tmp
	
                tmp_jobid=$( ${BATCH} ${WPDIR}/MONARC/MOORINGS/job_${CASE}_MOOR_y${C_YEAR}${C_YEAR}.tmp )

                if [ ${MACHINE} == 'occigen' ]   ;  then jobid=$( echo $tmp_jobid | awk '{print $4}'  ) ; else jobid=$tmp_jobid   ; fi
	        
	        if [ $iyear -eq 1 ] ;  then jobids="${jobid}" ; else  jobids="${jobids}:${jobid}" ; fi

	        let iyear=$iyear+1
		let C_YEAR=$C_YEAR+1

	done

	
	cd ../
fi

#---------------------------------------------------------------------------------------------------------------------------------
#--  SECTIONS  -------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------
if [ ${DOSECTIONS} -eq 1 ]  ; then 

	cd ./SECTIONS/
	
	LGTS_S=${REF_YEAR}
	LGTS_E=${REF_YEAR}
	C_YEAR=${STA_YEAR}
	DEPEND="#"
	iyear=1
        ZNCDF=False
	if [ ${OUTNCDF} = 1 ] ; then  
              ZNCDF=True
	fi

	while [ ${C_YEAR} -le ${END_YEAR}  ]  ;  do

	      if [ ${C_YEAR} -eq ${END_YEAR}  ] ; then LGTS_E=${END_YEAR} ; fi
	
	      # Build Python script for each plot
	      sed -e "s/XXSYEAXX/${C_YEAR}/" \
	          -e "s/XXEYEAXX/${C_YEAR}/" \
	          -e "s/XXCONFXX/${CONFIG}/" \
	          -e "s/XXCASEXX/${CASE}/" \
	          -e "s/XXCASE2XX/${CASE2}/" \
                  -e "s/XXXIOSFREQXX/'${XIOS}'/" \
	          -e "s/XXLGTSSXX/${LGTS_S}/" \
	          -e "s/XXLGTSEXX/${LGTS_E}/" \
	          -e "s/XXDIAGSSECXX/${MKDIAGSSEC}/" \
	          -e "s/XXSECPLOTXX/${MKSECPLOT}/" \
	          -e "s/XXNCDFOUTXX/${ZNCDF}/" CREG_sections.py > ${WPDIR}/MONARC/SECTIONS/CREG_sections_y${C_YEAR}${C_YEAR}.py
	      chmod 750 ${WPDIR}/MONARC/SECTIONS/CREG_sections_y${C_YEAR}${C_YEAR}.py
	
	      # Build the job which will be launched (one per grid)
              if [ ${C_YEAR} -eq ${END_YEAR}  ] && [ $(( $END_YEAR-$STA_YEAR )) -gt 0 ] ; then 
                   if [ ${MACHINE} == 'occigen' ]  ; then  DEPEND="#SBATCH --dependency=afterany:${jobids}" ; fi
                   if [ ${MACHINE} == 'datarmor' ] ; then  DEPEND="#PBS -W depend=afterany:${jobids}"        ; fi
              fi
              
	      if [ ! -d ${WPDIR}/MONARC/SECTIONS/NETCDF ] ; then mkdir ${WPDIR}/MONARC/SECTIONS/NETCDF ; fi

	      sed -e "s/XXSYEAXX/${C_YEAR}/" \
	          -e "s/XXEYEAXX/${C_YEAR}/" \
	          -e "s/XXCONFXX/${CONFIG}/" \
	          -e "s/XXCASEXX/${CASE}/" \
	          -e "s/XXCASE2XX/${CASE2}/" \
	          -e "s;XXDATA_DIRXX;${DATA_DIR};" \
	          -e "s;XXGRID_DIRXX;${GRID_DIR};" \
	          -e "s;XXDATA_DIR2XX;${DATA_DIR2};" \
                  -e "s/XXXIOSFREQXX/'${XIOS}'/" \
	          -e "s;XXINITDIRXX;${INITDIR};" \
	          -e "s/XXGETFILEXX/${GETFILE}/"  job_sections.ksh > ${WPDIR}/MONARC/SECTIONS/tmp_job_sections

	       sed -e "s/XXDEPENDENCYXX/${DEPEND}/" \
	           -e "s/XXSPTYPEXX/SECTI/" \
	           -e "s/XXPTYPEXX/SECTIONS/" \
                   -e "s/XXSYEAXX/${C_YEAR}/" \
                   -e "s/XXEYEAXX/${C_YEAR}/" \
                   -e "s/XXCASEXX/${CASE}/" \
	           -e "s;XXINITDIRXX;${INITDIR};" \
	           -e "s/XXTIMEXX/40/" ../JHEAD_${MACHINE}.bash > ${WPDIR}/MONARC/SECTIONS/tmp_head_sections

               cat ${WPDIR}/MONARC/SECTIONS/tmp_head_sections ${WPDIR}/MONARC/SECTIONS/tmp_job_sections > ${WPDIR}/MONARC/SECTIONS/tmp_job_${CASE}_sections_y${C_YEAR}${C_YEAR}
	
              tmp_jobid=$( ${BATCH} ${WPDIR}/MONARC/SECTIONS/tmp_job_${CASE}_sections_y${C_YEAR}${C_YEAR} )

              if [ ${MACHINE} == 'occigen' ]   ;  then jobid=$( echo $tmp_jobid | awk '{print $4}'  ) ; else jobid=$tmp_jobid   ; fi

	      if [ $iyear -eq 1 ] ;  then jobids="${jobid}" ; else  jobids="${jobids}:${jobid}" ; fi
	      let iyear=$iyear+1

	      let C_YEAR=$C_YEAR+1
	
	done
	
	cd ../
fi
	
#---------------------------------------------------------------------------------------------------------------------------------
#--  INTQUANT  -------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------
if [ ${DOINTQUANT} -eq 1 ]  ; then 

	cd ./INTQUANT/
	
	LGTS_S=${REF_YEAR}
	LGTS_E=${REF_YEAR}
	C_YEAR=${STA_YEAR}
	DEPEND="#"
	iyear=1
        ZNCDF=False
	if [ ${OUTNCDF} = 1 ] ; then  
              ZNCDF=True
	fi

	while [ ${C_YEAR} -le ${END_YEAR}  ]  ;  do

	      if [ ${C_YEAR} -eq ${END_YEAR}  ] ; then LGTS_E=${END_YEAR} ; fi
	
	      sed -e "s/XXSYEAXX/${C_YEAR}/" \
	          -e "s/XXEYEAXX/${C_YEAR}/" \
	          -e "s/XXCONFXX/${CONFIG}/" \
	          -e "s/XXCASEXX/${CASE}/" \
	          -e "s/XXCASE2XX/${CASE2}/" \
                  -e "s/XXXIOSFREQXX/'${XIOS}'/" \
	          -e "s/XXLGTSSXX/${LGTS_S}/" \
	          -e "s/XXLGTSEXX/${LGTS_E}/" \
	          -e "s/XXNCDFOUTXX/${ZNCDF}/" CREG_intquant.py > ${WPDIR}/MONARC/INTQUANT/CREG_intquant_y${C_YEAR}${C_YEAR}.py
	      chmod 750 ${WPDIR}/MONARC/INTQUANT/CREG_intquant_y${C_YEAR}${C_YEAR}.py
	
	      # Build the job which will be launched (one per grid)
              if [ ${C_YEAR} -eq ${END_YEAR}  ] && [ $(( $END_YEAR-$STA_YEAR )) -gt 0 ] ; then 
                   if [ ${MACHINE} == 'occigen' ]  ; then  DEPEND="#SBATCH --dependency=afterany:${jobids}" ; fi
                   if [ ${MACHINE} == 'datarmor' ] ; then  DEPEND="#PBS -W depend=afterany:${jobids}"        ; fi
              fi

	      if [ ! -d ${WPDIR}/MONARC/INTQUANT/NETCDF ] ; then mkdir ${WPDIR}/MONARC/INTQUANT/NETCDF ; fi
              
	      sed -e "s/XXSYEAXX/${C_YEAR}/" \
	          -e "s/XXEYEAXX/${C_YEAR}/" \
	          -e "s/XXCONFXX/${CONFIG}/" \
	          -e "s/XXCASEXX/${CASE}/" \
	          -e "s/XXCASE2XX/${CASE2}/" \
	          -e "s;XXDATA_DIRXX;${DATA_DIR};" \
	          -e "s;XXGRID_DIRXX;${GRID_DIR};" \
	          -e "s;XXDATA_DIR2XX;${DATA_DIR2};" \
                  -e "s/XXXIOSFREQXX/'${XIOS}'/" \
	          -e "s;XXINITDIRXX;${INITDIR};" \
	          -e "s/XXGETFILEXX/${GETFILE}/"  job_intquant.ksh > ${WPDIR}/MONARC/INTQUANT/tmp_job_intquant

	       sed -e "s/XXDEPENDENCYXX/${DEPEND}/" \
	           -e "s/XXSPTYPEXX/INTQT/" \
	           -e "s/XXPTYPEXX/INTQUANT/" \
                   -e "s/XXSYEAXX/${C_YEAR}/" \
                   -e "s/XXEYEAXX/${C_YEAR}/" \
                   -e "s/XXCASEXX/${CASE}/" \
	           -e "s;XXINITDIRXX;${INITDIR};" \
	           -e "s/XXTIMEXX/40/" ../JHEAD_${MACHINE}.bash > ${WPDIR}/MONARC/INTQUANT/tmp_head_intquant

               cat ${WPDIR}/MONARC/INTQUANT/tmp_head_intquant ${WPDIR}/MONARC/INTQUANT/tmp_job_intquant > ${WPDIR}/MONARC/INTQUANT/tmp_job_${CASE}_intquant_y${C_YEAR}${C_YEAR}
	
              tmp_jobid=$( ${BATCH} ${WPDIR}/MONARC/INTQUANT/tmp_job_${CASE}_intquant_y${C_YEAR}${C_YEAR} )

              if [ ${MACHINE} == 'occigen' ]   ;  then jobid=$( echo $tmp_jobid | awk '{print $4}'  ) ; else jobid=$tmp_jobid   ; fi

	      if [ $iyear -eq 1 ] ;  then jobids="${jobid}" ; else  jobids="${jobids}:${jobid}" ; fi
	      let iyear=$iyear+1

	      let C_YEAR=$C_YEAR+1
	
	done
	
	cd ../
fi
	

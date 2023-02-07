import os 

def chkfile(chk_file,zstop=False,zscript=None):

	file_pres = os.path.exists(chk_file)
	if not file_pres : 
		print " 			>>>>>> WARNING: Check the path or the file name syntax or presence  " 
		print "           		>>>>>>          ", chk_file
		print 
		if zstop :
			print " 			#######################################################################"
			print " 			#######################################################################"
			print " 			>>>>>> THE FILE ", chk_file, " IS MISSING "
			print " 			>>>>>> PYTHON SCRIPT:", zscript, " IS STOPPED "
			print " 			#######################################################################"
			print " 			#######################################################################"

	return file_pres

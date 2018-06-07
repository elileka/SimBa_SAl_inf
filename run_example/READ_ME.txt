############################################### Note 1 ##################################################
# This directory contains what you need for a small run example.                                        #
# The paths are relative to the working dir. It's possible to replace them with full paths.             #
# Place the run_example directory in your working directory.                                            #
# If you build a Visual Studio solution, use this link to figure out where the working dir is:          #
# https://stackoverflow.com/questions/3423538/how-can-i-find-out-the-value-of-projectdir                #
#########################################################################################################

############################################### Note 2 ##################################################
# This example uses "full" optimization, meaning a dp procedure will be computed with each input table. #
# The example relies on a list file pointing at two table files. The structure of the list file is:     #
# <TABLE_FILE_OUTPUT_OF_SIMBA_SAL_SIM>	<r_val>	<mu_val>	<gamma_val>	<max_indel_length>	<t>         #
# The fileds are tab-separated. Each field contains the value of that parameter in the tables files.    #
# If you produce tables using SimBa_SAl_sim, you should list this way.                                  #
# A Perl script to list tables is provided to assist with listing tables.                               #
#########################################################################################################

############################################### Note 3 ##################################################
# To get to know more of the SimBa_SAl_inf parameters, run the program without arguments.               #
# For the theory and algorithmic details we kindly refer you to the manuscript:                         #
# Levy Karin, Ashkenazy, Hein, and Pupko, 2018.                                                         #
# and to the code documentation file.                                                                   #
#########################################################################################################

##################################### The SimBa_SAl_inf command #########################################
# In Visual Studio, put this in the program arguments:                                                  #
input_fasta_file=.\run_example\HOMSTRAD_true_MIP.fas file_with_paths_to_chop_tables=.\run_example\list_of_tables.txt type_seq=AA path_to_chebi_JTT_coef_file=.\run_example\chebi_coef.txt type_file=aligned out_files_prefix=.\run_example\EXAMPLE_OUT_ opt_method=full
# In UNIX:                                                                                              #
<PATH_TO_SIMBA_SAL_INF_PROGRAM> input_fasta_file=.\run_example\HOMSTRAD_true_MIP.fas file_with_paths_to_chop_tables=.\run_example\list_of_tables.txt type_seq=AA path_to_chebi_JTT_coef_file=.\run_example\chebi_coef.txt type_file=aligned out_files_prefix=.\run_example\EXAMPLE_OUT_ opt_method=full
#########################################################################################################

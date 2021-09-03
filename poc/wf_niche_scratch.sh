            #######################################
            #--     Individual Niches Scratch   --#
            #--         Scott Yanco, Phd        --#
            #--       scott.yanco@yale.edu      --#
            #######################################

# Proof of concept (POC) workflow script for processing individual, time-dynamic
# Grinellian niches.  Script pulls data from movebank.org into a SQL database,
# processes output and gets niche estimation.

# TODO: modify /src/poc/ filepaths after workflow finalization

#-----------------------------------------------------------------#

##-- BUILD DATABASE --##

# Build and load SQL db with desired data from movebank.org
# Study's to process can be controled in the study.csv control file
# located at /$PD/analysis/movedb/ctfs/
~/projects/niche_scratch/src/poc/create_db.sh

##-- PROCESS MOVEMENT DATA --##

# TODO:  Go back and add segemntation and filtering steps;
# For now, proceed to annotation and niche estimation


    	##-- ANNOTATE LOCATIONS --##

    #- Create csv for annotation

# Set working directory
wd=/home/syanco/projects/niche_scratch/

#set analysis path for output
out=~/projects/niche_scratch/analysis

# Change to wd
cd $wd

# Create csv from db
~/projects/niche_scratch/src/poc/create_MOL.r $out/dat_MOL.csv -t

    #- Upload to MOL - this happens outside workflow (for now)

# TODO: Check with MOL team to see if this can be done programmatically 
# DONE:  It cannot apparently - thus for now the plan is to try the on-the-fly annotator 
#           in STOAT - we'll see if it's compuationally reasonable...

    #- Annotate Tracks - ON-THE-FLY ANNOTATOR
~/projects/niche_scratch/src/poc/annotate_tracks.r $out/dat_MOL.csv $out/dat_anno.csv -t



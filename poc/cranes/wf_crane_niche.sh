#######################################
#--     Individual Niches Scratch   --#
#--           1000 Cranes           --#
#--         Scott Yanco, Phd        --#
#--       scott.yanco@yale.edu      --#
#######################################

# Proof of concept (POC) workflow script for processing individual, time-dynamic
# Grinellian niches.  

# TODO: modify /src/poc/cranes filepaths after workflow finalization

#-----------------------------------------------------------------#

#set working directory
wd=~/projects/niche_scratch

#move to WD
cd $wd


##-- Join annotations --##

# Activate conda env
conda activate db

# Run annotation join script - results returned to db
# NOTE: control file must be updated to select which variables are joined
Rscript $wd/src/poc/cranes/join_annos.r $wd/data/anno_move.db $wd/data/anno_move.db $wd/analysis/cranes/ctfs/anno_vars.csv 


##-- Calculate time-dynamic niche metrics --##

  #-  Niche Size
#activate conda env
conda activate MVNH_db

# Calculate  weekly niche size estimates using `MVNH`
Rscript $wd/src/poc/cranes/calc_niche_sizes.r $wd/data/anno_move.db $wd/analysis/cranes/niche_sizes.csv $wd/analysis/cranes/ctfs/anno_vars.csv anno_join_2021-08-23

  #- Nich Dissimilarity
# Calculate weekly niche dissimilarity estimates using `MVNH`
# NOTE: this output needs to be saved as '.Rdata' not '.csv' b/c some cells are 
  # actually vectors of output from the `MVNH` package
Rscript $wd/src/poc/cranes/calc_niche_dissim.r $wd/data/anno_move.db $wd/analysis/cranes/niche_dissims.Rdata $wd/analysis/cranes/ctfs/anno_vars.csv anno_join_2021-08-23


##-- Calculate NSD --##

#activate conda env
# TODO: consider using this same env in the niche steps, 
  # this is same env as MVNH_db, but includes `amt`...
conda activate amt_db

# Calculate weekly NSD for each individual
Rscript $wd/src/poc/cranes/calc_nsd.r $wd/data/anno_move.db $wd/analysis/cranes/nsd.csv

##-- Plots --##

#activate conda env
conda activate plots

#produce plots across individuals
Rscript $wd/src/poc/cranes/plot_niches.r $wd/analysis/cranes $wd/analysis/cranes

#stitch pdfs together
cd analysis/cranes/niche_plots
pdftk *.pdf cat output crane_niches.pdf
cd $wd
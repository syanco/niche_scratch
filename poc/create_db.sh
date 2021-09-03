#!/usr/bin/env bash

#-- Create SQL DB and download Movebank Data --#
#-- Scott Yanco; scott.yanco@yale.edu --#

# Assumes user has downloaded mosey_db from https://github.com/benscarlson/mosey_db
# Also assumes user is using Breezy file structure as in: https://github.com/benscarlson/breezy

#project directory
PD=/home/syanco/projects/niche_scratch

#create dir for db
mkdir $PD/analysis/movedb

#set working dir to the db dir inside /analysis
wd=$PD/analysis/movedb

#set var for path to mosey_db scripts
export MOVEDB_SRC=/home/syanco/tools/mosey_db

#move to db dir
cd $wd

#create sub-dir for data
mkdir -p $wd/data

#create the db
cat $MOVEDB_SRC/db/create_db.sql | sqlite3 $wd/data/move.db

#define csvdir (location for raw and cleaned gps output) given as arg to load_studies.sh
csvdir=$wd/data

# Load desired studies into db - studies are defined in the study control file ('wd/ctfs/study.csv')
# requires auth.yml with Movebank UN and PW - auth.yml is excluded from this repo for security purposes
source $MOVEDB_SRC/db/load_studies.sh $csvdir
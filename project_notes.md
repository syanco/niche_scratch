## TODO
* Join annotations from the database
* Fix the niche estimator script to loop through individuals
  * Lots to sort here:
    * What's the baseline for dissimilarity calculations?
    * Common set of niceh variables?
* Put conda info in the README
  * export relevant .yml files
  * include instructions for creating conda env from yml


**OLD TODO**
* Finish annotate_tracks.r script
  *Confirm ability to pass variables from ctfs to annotator.
  *rstoat throwing error with particular products
* Next script is to feed annotated data into MVNH niche estimator

## Activity log

|Date|Activity|
|:-|:------------|
|2021-07-20|Pilot WF completed up to annotator call|
|2021-08-06|Setting aside old wf dev b/c we have a fully-annotated crane dataset, proceeding with new wf building from that|

## Notes
* Cranes pre-anotated by STOAT team at Map of Life.

**OLD NOTES**
* The anotator relies on breaking data up into chunks of 1000 - this could probably be parallelized easily for large datasets
  * Consider writing a version of the annotate_tracks.r script for the HPC
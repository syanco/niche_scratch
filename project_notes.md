## TODO
* Plotting script to combine NSD, velocity, niche size and dissim in plots for every individual
* Need to come up with "final" and well-justified set of niche variables  
  * *and* be able to annotate with those vars...  
    *May require modes to `join_annos.r` script
* Think about some sort of spatial size measure on same weekly time scales - do niche dynamics conform to spatial dynamics?  
*Find a way to aggregate/summarize relationships across individuals so I can make inter-species comparisons  
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
|2021-08-16|Crane workflow built through flexible joining of annotations and time-dynamic niche size estimation.|
|2021-08-18|Working scripts to produce NSD, velocity, ncihe size, niche dissim!|

## Notes
* Cranes pre-anotated by STOAT team at Map of Life.  
* Questions:  
  * What's the baseline for dissimilarity calculations?
  * Common set of niche variables?


**OLD NOTES**
* The anotator relies on breaking data up into chunks of 1000 - this could probably be parallelized easily for large datasets
  * Consider writing a version of the annotate_tracks.r script for the HPC
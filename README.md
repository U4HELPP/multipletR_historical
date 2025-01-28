# multipletR: Multiplet Classification for Single Cell RNASeq PDX Samples

This project is focused on developing a new method for multiplet classification of PDX single cell sequencing samples using the 10X GEM counts.

## Description

__The Challenge:__ 
PDX sequencing data contains both human and mouse genomic content that must be seperated prior to downstream bioinformatics analysis. In single cell data, each cell barcode should ideally relate to a single cell that is either of human or mouse origin.  However, sometimes a single GEM droplet contains cells from both human and mouse, resulting in a "multiplet". Methods are needed to detect these multiplets and remove them from the analysis as they could contaminate downstream analyses.  

__Current Process:__ 
The 10X pipeline automatically classifies each barcode/cell as Human, Mouse, or a Multiplet and stores those calls in a file named "GEM_classification.csv". The algorithm is based on setting a threshold using the 10th percentile values, and is described on 10X's website [HERE](https://kb.10xgenomics.com/hc/en-us/articles/115003517183-How-does-cellranger-count-identify-multiplets).

__The Problem:__ 
The current 10X GEM classification method assumes there is a 1:1 ratio of human and mouse cells in the sample.  Based on the data being generated from the U4HELPP project this is very far from reality. In the current PDX data, many samples are seeing a large amount of clearly human cells being classified as multiplets.

__The Solution:__
We are developing new methods that can utilize the 10X GEM count data and do a better job seperating multiplets from human or mouse cells.

## Key Resources and Links

__Google Drive:__

* [Agenda and Meeting Notes](https://docs.google.com/document/d/1vdZ5HJRoQ1ginSqaSBghWv4SpT6nFEB82Rzk9n_Y9yA/edit?tab=t.0)

* Results are stored and shared on Google Drive in the following directory: [0_CHarrell/2023_U4HELPP_LivingPDXProgram/Data/scRNASeq/MESSY-Multiplet_Evaluation_for_ScSeq_barnYard_samples](https://drive.google.com/drive/u/0/folders/18DVetp2fzvMKtdTw3CEFbqsp3y_YA_C1)

__Apollo:__

* Additional data, scripts, and results are also on Apollo (/lustre/home/harrell_lab/scRNASeq/exploratory_analyses/multiplet_exploration).  This includes actual results from runs using the semi-log space thresholding method that were used in data analysis for the U4HELPP project while the development of this method is in progress.

__10X Genomics Information:__

* [Current automated species classification method.](https://kb.10xgenomics.com/hc/en-us/articles/115003517183-How-does-cellranger-count-identify-multiplets)

* [Manual Loupe classification method if automated method fails.](https://kb.10xgenomics.com/hc/en-us/articles/28552478030093-Species-Cell-Assignments-on-Loupe-for-Barnyard-samples)

* [Benchmark 10X Datasets](https://www.10xgenomics.com/datasets?configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bplatform%5D%5B0%5D=Chromium%20Single%20Cell&refinementList%5Bspecies%5D%5B0%5D=Human%2C%20Mouse)



## Getting Started

### Dependencies

* R (version dependency?)

### Git Team Workflow

All team members should adhere to the following protocol to keep the main branch clean and functional.

* Create your own development branch to work on.
* Update your development branch as needed until ready to merge with main.
* Create a pull request to merge with main branch and resolve any conflicts.
* Merge with main branch and tag the commit with a minor version number.
* Notify team to update their local repositories.
* Create a new dev branch, or continue editing on your existing branch.
* Tag major releases as needed, or as new functionality is added.

### Installing (TBD)

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program (TBD)

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```


## Authors

* Amy Olex
* Mikhail Dozmorov
* Brydon Wall

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the GNU GPLv3 License - see the LICENSE.md file for details

## Acknowledgments

# metstats
Library for standardized analysis and visualization of metabolomics data

Calico Life Sciences, LLC

[Demo analysis of Metabolomics Workbench ST003519 published in quarto pub](https://delfarah.quarto.pub/metstats-demo---metabolomics-workbench-st003519-dd5a/)

# Analysis workflow

1. Convert mass spec raw files to mzML using [ProteWizard MSConvert](https://proteowizard.sourceforge.io/tools/msconvert.html) software.
e.g. Metabolomics Workbench [ST003519](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&DataMode=AllData&StudyID=ST003519&StudyType=MS&ResultType=1#DataTabs) for demo, raw files are already mzML.

2. Run dataset through [MAVEN peakdetector](https://github.com/eugenemel/maven/releases) to generate ```mzrollDB``` file for peak identification and quantification.

3. Open mzrollDB in [MAVEN 2 software](https://www.mdpi.com/2218-1989/12/8/684) and annotate peaks.

4. Create a project folder and put the following files in the folder:
- Annotated ```mzrollDB``` file
- qmd file ```metstats/vignettes/ST003519_demo.Qmd``` for the demo

5. In terminal go to the project folder and run the following to set up a new Quarto website project:
```
quarto create project my-website
```
- Choose a folder name for the project/website.

6. Render the qmd to generate html webpage in localhost.

7. To publish html to Quarto pub, run the following in terminal:
```
quarto publish quarto-pub ST003519_demo.qmd
```
- a [Quarto pub](https://quartopub.com/) account is required for publishing.

# Network Analysis

To run Cytoscape and generate metabolic networks:

1. Download latest version of [Cytoscape](https://cytoscape.org/download.html)

2. Install [RCy3 2.5.1](https://github.com/delfarahalireza/RCy3) Cytoscape-R API
- Newer versions of RCy3 are problematic.
- Recommended method of installation:

```r
# Clone the RCy3 repository via terminal:
git clone git@github.com:delfarahalireza/RCy3.git

# Install the package from source in R:
install.packages("cloned directory of RCy3", repos = NULL, type = "source")

# Restart R session and load RCy3
library(RCy3)

# Check that version of RCy3 is 2.5.1 after loading library(Rcy3)
packageVersion("RCy3")

```

3. Open ```metstats/inst/extdata/Networks.cys``` file

4. In R, run the ```network plots```, ```pathways color```, ```pathways shape```, and ```pathways``` chunks similar to ```metstats/vignettes/ST003519_demo.Qmd```

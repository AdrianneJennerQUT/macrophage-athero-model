This code accompanies the manuscript *"Investigating atherosclerotic plaque development from macrophage lineages using an agent-based model"* by Adrianne L Jenner<sup>1</sup>, Angela M. Reynolds<sup>2</sup>, Keith L. Chambers<sup>3</sup>, Andrew J. Murphy<sup>4</sup>, Helen Byrne<sup>3,5</sup>, Mary R. Myerscough<sup>6</sup>

<sup>1</sup>School of Mathematical Sciences, Queensland University of Technology, Brisbane, QLD, Australia,
<sup>2</sup>Virginia Commonwealth University, Richmond, USA
<sup>3</sup>Mathematical Institute, The University of Oxford, Oxford, UK 
<sup>4</sup Baker Heart and Diabetes Institute, Melbourne, Australia
<sup>5</sup>Ludwig Institute for Cancer Research, The University of Oxford, Oxford, UK
<sup>6</sup> School of Mathematics and Statistics, The University of Sydney, Sydney, Australia

![ControlGrowth](https://user-images.githubusercontent.com/85978071/152262842-950fbb47-dcd1-4da0-9c19-1c6bae2230b0.png)

In this repository are two ABM models built in PhysiCell for atherosclerotic plaque developement. The firs is in "simple_macropahge_model" and the second is in "athero". 
These two models only differ by their lipid model. The first has a simple homogeneous IC for lipid. The second has that lipid is secreted by ghost cells, where
secretion is a function of the distance from the top boundary, i.e. cells secrete lipid the further they are from the blood. 

To run the models in PhysiCell, download the "simple_macrophage_model" folder and "athero" folder and save in the "sample_projects" folder. 
Then download "Makefile-default" and save in the "sample_projects folder" overwriting the current version in this folder

To run the first model, open your command line and change directory into the PhysiCell folder, then run: 

make reset,
make simple-macrophage-model-sample,
make

if on windows, then type 
.\simple_macrophage_model_sample.exe

if on mac, then type
./simple_macrophage_model_sample

To run the second model, open your command line and change directory into the PhysiCell folder, then run: 

make reset,
make athero-sample,
make

if on windows, then type 
.\athero.exe

if on mac, then type
./athero

Download the matlab script to plot the video of the microenvironment and save in the "matlab" folder in the main PhysiCell branch.

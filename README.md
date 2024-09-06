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

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
		
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
		
	initialize_cell_definitions_from_pugixml(); 
		
	Cell_Definition* pMacrophage = find_cell_definition( "macrophage" ); 
	
	pMacrophage->functions.update_phenotype = macrophage_function; 
	
	pMacrophage->phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	build_cell_definitions_maps(); 

	setup_signal_behavior_dictionaries(); 	

	setup_cell_rules(); 

	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	
	// initialise one macrophage at the top of the domain
	Cell* pC;
	Cell_Definition* pCD = find_cell_definition( "macrophage" );
	
	// set initial x and y position for cell
	double x = 0;
	double y = 140;
	pC = create_cell( *pCD ); 
	pC->assign_position( x,y, 0.0 );
		
	return; 
}

std::vector<std::string> coloring_function( Cell* pCell )
{
	// default mac colour
	// output[0] = nucleus colour
	// output[1] = boundary of nucleus colour
	// output[2] = cytoplasm
	// output[3] = exterior boundary
	std::vector<std::string> output = { "blue" , "black" , "blue", "black" }; 
	
	// check if cell is dead
	if( pCell->phenotype.death.dead == true )
	{
		// if dead, change nucleus and cytoplasm to be red and dark red
		 output[0] = "red";  
		 output[2] = "darkred"; 
		 
		 return output; 
	}
		
	return output; 
}

std::vector<Cell*> get_possible_neighbors( Cell* pCell )
{
	std::vector<Cell*> neighbors = {}; 

	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end =
		pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ neighbors.push_back( *neighbor ); }

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{ neighbors.push_back( *neighbor ); }
	}
	
	return neighbors; 
}

void macrophage_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	// bookkeeping 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	static int lipid_index = microenvironment.find_density_index( "lipid" ); 
		
	double u_max = parameters.doubles("u_max"); // max lipid uptake rate
	double lipid_half = parameters.doubles("lipid_half"); // lipid uptake half-effect 	
	double lipid_internal = pCell->phenotype.molecular.internalized_total_substrates[lipid_index];
	
	// check my internal lipid and if it's above a certain threshold, stop uptaking lipid and become apoptotic
	if( lipid_internal> parameters.doubles("lipid_threshold"))
	{
		std::cout<<"I died from eating too much"<<std::endl;
		pCell->start_death( apoptosis_model_index );
		pCell->functions.update_phenotype = NULL; 
	}	
		
	// check how much lipid I've take up and based on that, change my speed
	double max_speed = parameters.doubles("macrophage_max_speed");
	double speed_threshold = parameters.doubles("macrophage_speed_threshold");
	
	
	if( lipid_internal<speed_threshold )
	{
		double new_speed = max_speed-max_speed*lipid_internal/speed_threshold;
		set_single_behavior( pCell, "migration speed" , new_speed ); 
	}
	else
	{
		//set_single_behavior( pCell, "migration speed" ,0);
	}
	
	// check for contact with a cell and then see if I eat that cell
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);
	
	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		// if it is not me 
		if( pTestCell != pCell)
		{
			// calculate distance to the cell 
			std::vector<double> displacement = pTestCell->position;
			displacement -= pCell->position;
			double distance = norm( displacement ); 
			
			double max_distance = pCell->phenotype.geometry.radius + 
				pTestCell->phenotype.geometry.radius; 
			max_distance *= 1.05; 
			
			double rand_sample = UniformRandom();
			
			// if the cell is dead and close enough and probability is satisfied
			if( distance < max_distance && pTestCell->phenotype.death.dead == true &&
				rand_sample<parameters.doubles("prob_engulf"))
			{
				std::cout << "\t\tnom nom nom" << std::endl; 
				pCell->ingest_cell( pTestCell ); 
			}
			
		
		}
	}
	
	
	return; 
}

void macrophage_arrival( double dt )
{
	// sampling a random number from U(0,1)
	double rand_sample = UniformRandom();
	
	// if arrival probability satisfied then 1 mac arrives
	if( rand_sample<parameters.doubles("macrophage_arrival_rate")*dt)
	{
		std::cout<<"Adding 1 macrophage"<<std::endl;
		
		static Cell_Definition* pCD = find_cell_definition( "macrophage" );
		Cell* pC = create_cell( *pCD ); 
		double length_x = microenvironment.mesh.bounding_box[3] - 
			microenvironment.mesh.bounding_box[0]; 
		double length_y = microenvironment.mesh.bounding_box[4] - 
			microenvironment.mesh.bounding_box[1]; 

		double Ymin = parameters.doubles("min_y_entry");
		double Yrange = parameters.doubles("y_entry_range");
			
		// create macrophage at random position in top of domain 
		std::vector<double> position = {0,0,0}; 
		position[0] = microenvironment.mesh.bounding_box[0] + UniformRandom() * length_x; 
		position[1] = Ymin + UniformRandom()*Yrange; 

		pC->assign_position( position );
	}
	return;	
}

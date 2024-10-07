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
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
		
	initialize_cell_definitions_from_pugixml(); 

	Cell_Definition* pMacrophage_C1 = find_cell_definition( "macrophage C1" ); 
	Cell_Definition* pMacrophage_C2 = find_cell_definition( "macrophage C2" ); 
	
	static int lipid_index = microenvironment.find_density_index( "lipid" ); 		 
	
	
	pMacrophage_C1->functions.update_phenotype = macrophage_function; 
	pMacrophage_C1->phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	pMacrophage_C2->functions.update_phenotype = macrophage_function; 
	pMacrophage_C2->phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	pMacrophage_C1->phenotype.molecular.fraction_transferred_when_ingested[ lipid_index ] = 1;
	pMacrophage_C2->phenotype.molecular.fraction_transferred_when_ingested[ lipid_index ] = 1;	
	
	// add endogenous lipid			
	pMacrophage_C1->phenotype.molecular.internalized_total_substrates[lipid_index] = parameters.doubles("endogenous_lipid");
	pMacrophage_C2->phenotype.molecular.internalized_total_substrates[lipid_index] = parameters.doubles("endogenous_lipid");	
						
	Cell_Definition* pGhost = find_cell_definition( "ghost" ); 
	pGhost->functions.update_phenotype = ghost_secretion_model; 
    //pGhost->functions.contact_function = NULL; 	
	pGhost->phenotype.mechanics.cell_cell_adhesion_strength = 0;
	pGhost->phenotype.mechanics.cell_cell_repulsion_strength = 0;
	
			
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
	static int lipid_index = microenvironment.find_density_index( "lipid" ); 
	
	//initialise ghost cells throughout the domain
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	
	Cell* pC;
	
	Cell_Definition* pCD_ghost = find_cell_definition( "ghost"); 
	double number_x_dir = (Xmax-Xmin)/16;
	double number_y_dir = (Ymax-Ymin)/16;
	
	for( int k=0; k < number_x_dir; k++ )
	{
		double x = Xmin + k*16;
		
		for( int j = 0; j < number_y_dir; j++)
		{
			double y = Ymin + j*16;
			
			pC = create_cell( *pCD_ghost ); 
			pC->assign_position( x,y, 0.0 );
			pC->is_movable = false; 
			//pC->is_active = false;
			
			//setting unique secretion and saturation density targets based on distance to vessel wall
			pC->phenotype.secretion.secretion_rates[lipid_index] = parameters.doubles("base_secretion_rate_ghosts")*(Ymax-y)/Ymax;		
			pC->phenotype.secretion.saturation_densities[lipid_index] = parameters.doubles("base_secretion_target")*(Ymax-y)/Ymax;
		}
	}
	
	// initialise one macrophage at the top of the domain
	static Cell_Definition* pCD_mac_C1 = find_cell_definition( "macrophage C1"); 
	static Cell_Definition* pCD_mac_C2 = find_cell_definition( "macrophage C2"); 
	
	// set initial x and y position for cell
	double x = 650;
	double y = 280;
	
	
	pC = create_cell( *pCD_mac_C1 ); 
	pC->assign_position( x,y, 0.0 );
		
	// set initial x and y position for cell
	x = 250;
	y = 250;
	pC = create_cell( *pCD_mac_C2 ); 
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
	std::vector<std::string> output = { "slateblue" , "black" , "slateblue", "black" }; 
	
	
	static Cell_Definition* pCD_mac_C1 = find_cell_definition( "macrophage C1"); 
	static Cell_Definition* pCD_mac_C2 = find_cell_definition( "macrophage C2"); 
	
	static Cell_Definition* pCD_ghost = find_cell_definition( "ghost"); 
	
	// check if cell is dead
	if( pCell->phenotype.death.dead == true )
	{
		// if dead, change nucleus and cytoplasm to be red and dark red
		 output[0] = "goldenrod";  
		 output[1] = "black";
		 output[2] = "saddlebrown"; 
		 output[3] = "black";
		 
		 return output; 
	}
	else if( pCell->type == pCD_mac_C1->type)
	{
		 output[0] = "deepskyblue"; 
		 output[1] = "black";		 
		 output[2] = "deepskyblue"; 
		 output[3] = "black";
	}
	else if( pCell->type == pCD_ghost->type)
	{
		 output[0] = "white"; 
		 output[1] = "white";
		 output[2] = "white";
		 output[3] = "white"; 	
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
void ghost_secretion_model( Cell* pCell, Phenotype& phenotype, double dt )
{
		
	return;
	
}
void macrophage_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	// bookkeeping 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	static int lipid_index = microenvironment.find_density_index( "lipid" ); 
		
	static Cell_Definition* pCD_mac_C1 = find_cell_definition( "macrophage C1"); 
	
	double u_max = parameters.doubles("u_max"); // max lipid uptake rate
	double lipid_half = parameters.doubles("lipid_half"); // lipid uptake half-effect 	
	double lipid_internal = pCell->phenotype.molecular.internalized_total_substrates[lipid_index];
	
	//std::cout<<lipid_internal<<std::endl; // prints out internal lipid for each agent
		
	// add individual cell efflux of lipid
	if( pCell->type == pCD_mac_C1->type)
	{
		//!! don't secrete my endogenous lipid
		double s_base = parameters.doubles("base_secretion_rate");
		pCell->phenotype.molecular.internalized_total_substrates[lipid_index] = lipid_internal - s_base*lipid_internal*dt;
		lipid_internal = pCell->phenotype.molecular.internalized_total_substrates[lipid_index];
	
	}
	
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
				
				//check internal amount in cell to be eaten + amount in eating cell is not more than threshold
				double lipid_internal_apop = pTestCell->phenotype.molecular.internalized_total_substrates[lipid_index];
				if(lipid_internal_apop+lipid_internal<parameters.doubles("lipid_threshold"))
				{
					// can eat cell
					pCell->ingest_cell( pTestCell ); 
				}
				else
				{
					// can only eat some of cell
					nibble_cell( pTestCell, pCell); 
				}			
			}
			
		
		}
	}
	
	#pragma omp critical
	// check if I am near the bottom boundaru
	if (pCell->position[1]<30)
	{
		std::cout<<"I am free"<<std::endl;
		delete_cell(pCell);// as cell emigrates from plaque
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
		
		rand_sample = UniformRandom();
		if(rand_sample <parameters.doubles("prob_lineage_class")) // add C1
		{pCD = find_cell_definition( "macrophage C1" );}
		else // add C2
		{pCD = find_cell_definition( "macrophage C2" );}
		
		Cell* pC = create_cell( *pCD ); 
		double length_x = microenvironment.mesh.bounding_box[3] - 
			microenvironment.mesh.bounding_box[0]; 
		double length_y = microenvironment.mesh.bounding_box[4] - 
			microenvironment.mesh.bounding_box[1]; 

		double Ymin = parameters.doubles("min_y_entry");
		double Yrange = microenvironment.mesh.bounding_box[4]-Ymin;//parameters.doubles("y_entry_range");
			
		// create macrophage at random position in top of domain 
		std::vector<double> position = {0,0,0}; 
		position[0] = microenvironment.mesh.bounding_box[0] + UniformRandom() * length_x; 
		position[1] = Ymin + UniformRandom()*Yrange; 

		pC->assign_position( position );
	}
	return;	
}

void nibble_cell( Cell* pCell_to_eat, Cell* pCell )
{
	
	// don't ingest self 
	if( pCell_to_eat == pCell )
	{ return; } 
	
	// don't ingest a cell that's already ingested 
	if( pCell_to_eat->phenotype.volume.total < 1e-15 )
	{ return; } 
		
	// make this thread safe 
	#pragma omp critical
	{
		/*
		if( pCell_to_eat->phenotype.death.dead == true )
		{ std::cout << this->type_name << " (" << this << ")" << " eats dead " << pCell_to_eat->type_name << " (" << pCell_to_eat 
			<< ") of size " << pCell_to_eat->phenotype.volume.total << std::endl; }
		else
		{ std::cout << this->type_name << " (" << this << ")" << " eats live " << pCell_to_eat->type_name << " (" << pCell_to_eat 
			<< ") of size " << pCell_to_eat->phenotype.volume.total << std::endl; }
		*/


		// mark cell to eat as dead 
		pCell_to_eat->phenotype.death.dead = true; 
		
		// mark cell that is eating as dead
		pCell->phenotype.death.dead = true; 
			
		// set secretion and uptake to zero of cell to be eaten
		pCell_to_eat->phenotype.secretion.set_all_secretion_to_zero( );  
		pCell_to_eat->phenotype.secretion.set_all_uptake_to_zero( ); 
		
		// set secretion and uptake to zero of cell eating
		pCell->phenotype.secretion.set_all_secretion_to_zero( );  
		pCell->phenotype.secretion.set_all_uptake_to_zero( ); 
		
		// deactivate all custom function 
		pCell_to_eat->functions.custom_cell_rule = NULL; 
		pCell_to_eat->functions.update_phenotype = NULL; 
		pCell_to_eat->functions.contact_function = NULL; 
		
		// deactivate all custom function 
		pCell->functions.custom_cell_rule = NULL; 
		pCell->functions.update_phenotype = NULL; 
		pCell->functions.contact_function = NULL; 
		
		// should set volume fuction to NULL too! 
		pCell_to_eat->functions.volume_update_function = NULL; 
		pCell->functions.volume_update_function = NULL; 

		// set cell as movable and non-secreting 
		pCell_to_eat->is_movable = true; 
		pCell_to_eat->is_active = false; 
		pCell->is_movable = true; 
		pCell->is_active = false; 

		// absorb a fraction of the volume(s)
		static int lipid_index = microenvironment.find_density_index( "lipid" ); 
		double lipid_internal_dead_cell = pCell_to_eat->phenotype.molecular.internalized_total_substrates[lipid_index];
		double lipid_internal_eating_cell = pCell->phenotype.molecular.internalized_total_substrates[lipid_index];
			
		double fraction_to_absorb = (parameters.doubles("lipid_threshold")-lipid_internal_eating_cell)/lipid_internal_dead_cell;
				
		// absorb fluid volume (fraction into the cytoplasm) 
		pCell->phenotype.volume.cytoplasmic_fluid += fraction_to_absorb*pCell_to_eat->phenotype.volume.fluid; 
		pCell_to_eat->phenotype.volume.cytoplasmic_fluid = (1-fraction_to_absorb)*pCell_to_eat->phenotype.volume.cytoplasmic_fluid; 
		
		// absorb nuclear and cyto solid volume (fraction into the cytoplasm) 
		pCell->phenotype.volume.cytoplasmic_solid += fraction_to_absorb*pCell_to_eat->phenotype.volume.cytoplasmic_solid; 
		pCell_to_eat->phenotype.volume.cytoplasmic_solid = (1-fraction_to_absorb)*pCell_to_eat->phenotype.volume.cytoplasmic_solid; 
		
		pCell->phenotype.volume.cytoplasmic_solid += fraction_to_absorb*pCell_to_eat->phenotype.volume.nuclear_solid; 
		pCell_to_eat->phenotype.volume.cytoplasmic_solid = (1-fraction_to_absorb)*pCell_to_eat->phenotype.volume.nuclear_solid; 
		
		// consistency calculations 		
		pCell->phenotype.volume.fluid = pCell->phenotype.volume.nuclear_fluid + 
			pCell->phenotype.volume.cytoplasmic_fluid; 
		pCell_to_eat->phenotype.volume.fluid = pCell_to_eat->phenotype.volume.nuclear_fluid + 
			pCell_to_eat->phenotype.volume.cytoplasmic_fluid; 
		
		pCell->phenotype.volume.solid = pCell->phenotype.volume.cytoplasmic_solid + 
			pCell->phenotype.volume.nuclear_solid; 
		pCell_to_eat->phenotype.volume.solid = pCell_to_eat->phenotype.volume.cytoplasmic_solid + 
			pCell_to_eat->phenotype.volume.nuclear_solid; 
			
		// no change to nuclear volume (initially) 
		//pCell_to_eat->phenotype.volume.nuclear = 0.0; 
		//pCell_to_eat->phenotype.volume.nuclear_fluid = 0.0; 
		
		pCell->phenotype.volume.cytoplasmic = pCell->phenotype.volume.cytoplasmic_solid + 
			pCell->phenotype.volume.cytoplasmic_fluid; 
		pCell_to_eat->phenotype.volume.cytoplasmic = pCell_to_eat->phenotype.volume.cytoplasmic_solid + 
			pCell_to_eat->phenotype.volume.cytoplasmic_fluid; 
			
		pCell->phenotype.volume.total = pCell->phenotype.volume.nuclear + 
			pCell->phenotype.volume.cytoplasmic; 
		pCell_to_eat->phenotype.volume.total = pCell_to_eat->phenotype.volume.nuclear + 
			pCell_to_eat->phenotype.volume.cytoplasmic; 

		pCell->phenotype.volume.fluid_fraction = pCell->phenotype.volume.fluid / 
			(  pCell->phenotype.volume.total + 1e-16 ); 
		pCell_to_eat->phenotype.volume.fluid_fraction = pCell_to_eat->phenotype.volume.fluid / 
			(  pCell_to_eat->phenotype.volume.total + 1e-16 ); 
		
		pCell->phenotype.volume.cytoplasmic_to_nuclear_ratio = pCell->phenotype.volume.cytoplasmic_solid / 
			( pCell->phenotype.volume.nuclear_solid + 1e-16 );
		pCell_to_eat->phenotype.volume.cytoplasmic_to_nuclear_ratio = pCell_to_eat->phenotype.volume.cytoplasmic_solid / 
			( pCell_to_eat->phenotype.volume.nuclear_solid + 1e-16 );		
		
		// update corresponding BioFVM parameters (self-consistency) 
		pCell->set_total_volume( pCell->phenotype.volume.total ); 
		pCell_to_eat->set_total_volume( pCell_to_eat->phenotype.volume.total ); 
		
		// absorb the internalized substrates 
		
		// multiply by the fraction that is supposed to be ingested (for each substrate) 
		//*(pCell_to_eat->internalized_substrates) *= 
		//	*(pCell_to_eat->fraction_transferred_when_ingested); // 
		
		//*(pCell->internalized_substrates) += *(fraction_to_absorb*pCell_to_eat->internalized_substrates); 
		//*(pCell_to_eat->internalized_substrates) -= *(fraction_to_absorb*pCell_to_eat->internalized_substrates); 
		
		pCell->phenotype.molecular.internalized_total_substrates[lipid_index] += fraction_to_absorb*pCell_to_eat->phenotype.molecular.internalized_total_substrates[lipid_index]; 
		pCell_to_eat->phenotype.molecular.internalized_total_substrates[lipid_index] -= fraction_to_absorb*pCell_to_eat->phenotype.molecular.internalized_total_substrates[lipid_index]; 
		
		//static int n_substrates = pCell->internalized_substrates->size(); 
		//pCell_to_eat->internalized_substrates->assign( n_substrates , 0.0 ); 	
		
		// trigger removal from the simulation 
		// pCell_to_eat->die(); // I don't think this is safe if it's in an OpenMP loop 
		
		// flag it for removal 
		// pCell_to_eat->flag_for_removal(); 

		// remove all adhesions 
		// pCell_to_eat->remove_all_attached_cells();
		
	}

	// things that have their own thread safety 
	//pCell_to_eat->flag_for_removal();
	//pCell_to_eat->remove_all_attached_cells();
	//pCell_to_eat->remove_all_spring_attachments();
	
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	pCell->start_death( apoptosis_model_index );
	
	
	return; 
}
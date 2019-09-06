/* This file is part of the EpiCell library.
*
* The most recent release of EpiCell can be downloaded at
* <http://epicells.unige.ch/>
*
* EpiCell is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* EpiCell is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with EpiCell.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
	\file
	The main function simulates a rectangular tissue of hexagonal cells.
	After the creation of the tissue, it goes through an initial relaxation.
	External mechanical constraints are set. They are applied by elastic bottom, left and right "borders".
	One selected cell is submitted to an external force.
	The simulation stops when it reaches a given number of iterations.
	This simulation is mainly used to study tissue deformation and cell migration when subjected to an external force. \n
	Examples of execution arguments: \n
	(Ex.1) ./cellMovement 1 1 1e9 100e-12 -0.05 0.04 20 100 2e-7
*/
#include "epicell.h"
#include "epicell.hh"
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <vector>
//to create a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <ctime>

#define MAXIT 10000000
using namespace epc;
using namespace std;

typedef double T;
int main(int argc, char *argv[]){
	// Check number of arguments.
	if(argc != 10){
		std::cout << std::endl
			<< "Incorrect number of parameters." << std::endl
			<< "9 Parameters: nRep myseed K A0 lambdaNorm gammaNorm ncols nrows" << std::endl
			<< std::endl
			<< "nRep - number of replicates" << std::endl
			<< "myseed - seed for random number generator" << std::endl
			<< "K - elasticity of the cell area" << std::endl
			<< "A0 -  cell preferred area" << std::endl
			<< "lambdaNorm - normalized line tension of the edges" << std::endl
			<< "gammaNorm - normalized contractility of the cell perimeter" << std::endl
			<< "ncols - number of columns of rectangular epithelium" << std::endl
			<< "nrows - number of columns of rectangular epithelium" << std::endl
			<< "force - force to be applied on the cell" << std::endl
			<< std::endl;
		exit(1);
	}

	// Assign arguments to model parameters.
	const unsigned long nRep = atoi(argv[1]);
	unsigned long myseed = atoi(argv[2]);
	double K = atof(argv[3]);
	double A0 = atof(argv[4]);
	double lambdaNorm = atof(argv[5]);//give lambdaNorm as parameter
	double gammaNorm = atof(argv[6]);//give gammaNorm as parameter
	unsigned long ncols = atoi(argv[7]);
	unsigned long nrows = atoi(argv[8]);
	double force = atof(argv[9]); // ex. 2e-7

	// External constraints.
	double elasticCoef = 0; //5; //5e-4;

	double sdLambda = double(0.0);
	double sdGamma = double(0.0);

	double r0 = sqrt(2*A0/(3*sqrt(3)));//0.5;
	double dt = 0.1; //0.02; //0.1
	double damping = 1.0;
	double gamma = gammaNorm*(K*A0);
	double lambda = lambdaNorm*(K*pow(A0,1.5));

	double boundaryLambda = 0.0;

	bool activatedFunction = true; //true; false;
	int wantedNumberIterations = 50000;

	// TODO: change this into a vector
	int density_1 = 0;
	int density_2 = 0;
	int density_3 = 0;
	int density_4 = 0;
	int density_5 = 0;
	int density_6 = 0;
	int density_7 = 0;
	int density_8 = 0;
	int density_9 = 0;
	int density_10 = 0;

	//bool scylla = false;

	// Creating the path to the folder containing the results of the simulation.
	// --------------------------------------------------------------------------
	std::string workdir = "./tmp";
	// Check if the upper folder exists
	struct stat st = {0};
	if (stat(workdir.c_str(), &st) == -1) {
	    mkdir(workdir.c_str(), 0700);
	}

	std::ostringstream ossGammaLambdaN;

	ossGammaLambdaN << "gammaN_" << gammaNorm << "_lambdaN_" << lambdaNorm << "_F_" << force;
	workdir.append(std::string("/").append(ossGammaLambdaN.str()));
	//Check if the upper folder exists
	if (stat(workdir.c_str(), &st) == -1) {
	    mkdir(workdir.c_str(), 0700);
	}

	// Date and time.
	time_t now = time(0);
	tm *ltm = localtime(&now);
	std::ostringstream oss_time;
	oss_time << "_" << 1900+ltm->tm_year << "_" << 1+ltm->tm_mon << "_"<< ltm->tm_mday << "_" << ltm->tm_hour << "h" << ltm->tm_min << "m" << ltm->tm_sec ;

	// Create a default cell type.
	std::shared_ptr<CellType> type0 = std::shared_ptr<CellType>(new CellType());
	// Create a forced cell type: just to identify it on the visualization.
	std::shared_ptr<CellType> type1 = std::shared_ptr<CellType>(new CellType());

	unsigned long myRep = 1;
	while(myRep <= nRep) {

		std::ostringstream ossSeed;
		ossSeed << myseed;
		std::string ossArgStr = ossSeed.str();
		for (int i=4;i<8;i++){
			ossArgStr.append(std::string("_").append(std::string(argv[i])));
		}

		std::string workdir2 = std::string(workdir).append(std::string("/").append(ossArgStr)).append(oss_time.str());

		// Check if the upper folder exists
		if (stat(workdir2.c_str(), &st) == -1) {
		    mkdir(workdir2.c_str(), 0700);
		}

		// Create a text file containing all the parameter values used in this execution.
		std::string exec_cmd_folder = workdir2;
		std::ofstream exec_cmd_stream(exec_cmd_folder.append("/exec_cmd.txt").c_str());
		exec_cmd_stream << "\nSimulation" << std::endl;
		exec_cmd_stream << " seed = " << myseed << std::endl;
		exec_cmd_stream << "\nInitial Cells" << std::endl;
		exec_cmd_stream << " K = " << K << ", A0 = " << A0 << std::endl;
		exec_cmd_stream << " gammaNorm = " << gammaNorm << " ==> gamma = " <<  gamma << std::endl;
		exec_cmd_stream << " lambdaNorm = " << lambdaNorm << " ==> lambda = " << lambda << std::endl;
		exec_cmd_stream << " boundaryLambda = " << boundaryLambda << std::endl;
		exec_cmd_stream.close();

		std::cout << "nRep= "<<nRep<<std::endl;

		unsigned long iT = 0;
		auto dyn = std::shared_ptr<Dynamics<2>>(new Verlet<2>(damping));

		RNG::seed(myseed);
		NRNG::seed(myseed);

        // Create epithelium tissue.
        // -------------------------
		Epithelium<2> epithelium(HexagonalCellsRectangularCluster<2>(nrows,ncols,r0,K,gamma,lambda,sdGamma,sdLambda,dyn, false,iT));
		// Add default cell type type0 to the cell cluster.
		epithelium.getCellCluster().getCellTypesManager().addCellType(type0);
		epithelium.getCellCluster().getCellTypesManager().setLambdaBetweenCellTypes(lambda,type0,type0);
		// Assign type0 to all cells.
		for (auto c : epithelium.getCellCluster().getCells()){
			c->setType(*type0);
		}
		// Set the lambda of boundary edges. It is automatically assigned to all boundary edges.
		epithelium.getCellCluster().setBoundaryLambda(boundaryLambda);

		// Set boundary vertices.
		epithelium.getCellCluster().setBoundaryVertices();
		// std::cout << "Number of boundary vertices = " << epithelium.getCellCluster().getBoundaryVertices().size() << std::endl;

		//Convergence (threshold, maximum-window-data)
		double conv_relax_thresh = 1e-4;//double(mathConstants::ForceEpsilon);
		statistics::Convergence conv_relax(conv_relax_thresh, 1000/dt);
		double energy = dataAnalysis::computeNetEnergy(epithelium.getCellCluster());
		conv_relax.takeValue(energy);//using the energy of the tissue as a parameter to check the convergence


		// Initial tissue relaxation.
		//---------------------------
		std::cout << "Initial Relaxation" << std::endl;
		while((conv_relax.hasConverged() == false)&&(iT < MAXIT)){//the model will run until there is no significative change in the regularity of the tissue (by significative we mean that the stdev/mean of the retularity is not higher than an epsilon
	           	if (iT % (int)(100/dt) == 0){
        		    	std::cout << "Relaxation: iT = " << iT << ", converged = " << conv_relax.hasConverged()<< "(" << conv_relax.getNormalizedStd () << " ?< "<< conv_relax_thresh << ")" <<std::endl;
			}
			epithelium.performTimeStep(dt,false,false);
			iT++;
			energy = dataAnalysis::computeNetEnergy(epithelium.getCellCluster());
			conv_relax.takeValue(energy);//using the energy of the tissue as a parameter to check the convergence
			// VTK.
			if (iT % (int)(1000/dt) == 0){
                VTKwriter<2> vtk(epithelium.getCellCluster(), workdir2+"/initRelax"+std::to_string(iT)+".vtk");
                vtk.addScalarToCell("id",dataFunctionals2D::getId);
                vtk.addScalarToCell("degree",dataFunctionals2D::computeDegree);
                vtk.addScalarToCell("onMitosis",dataFunctionals2D::isOnMitosis);
                vtk.addScalarToCell("type",dataFunctionals2D::getType);
                vtk.addScalarToCell("age",dataFunctionals2D::computeAge);
                vtk.addScalarToCell("area",dataFunctionals2D::computeArea);
                vtk.addVectorToVertex("force",dataFunctionals2D::computeForce);
            }
		}//end-while of Relaxation
        std::cout << "Relaxation: iT = " << iT << ", converged = " << conv_relax.hasConverged()<< "(" << conv_relax.getNormalizedStd () << " ?< "<< conv_relax_thresh << ")" <<std::endl;

		// Compute average relaxed area A_relax. It can be used as a reference area for mitosis triggering, when the latter is used in simulation. Ex. Cells with area < 90% A_relax can not start mitosis.
        double A_relax = dataAnalysis::computeAverageCellArea(epithelium.getCellCluster().getCells());
        std::cout << "A_relax = " << A_relax << std::endl;

        VTKwriter<2> vtk(epithelium.getCellCluster(), workdir2+"/initRelax"+std::to_string(iT)+".vtk");
        vtk.addScalarToCell("id",dataFunctionals2D::getId);
        vtk.addScalarToCell("degree",dataFunctionals2D::computeDegree);
        vtk.addScalarToCell("onMitosis",dataFunctionals2D::isOnMitosis);
        vtk.addScalarToCell("type",dataFunctionals2D::getType);
        vtk.addScalarToCell("age",dataFunctionals2D::computeAge);
        vtk.addScalarToCell("area",dataFunctionals2D::computeArea);
        vtk.addVectorToVertex("force",dataFunctionals2D::computeForce);

        // Get from the relaxed tissue, the min and max values of two directions.
        double xmin = std::numeric_limits<double>::max();
				double xmax = std::numeric_limits<double>::min();
				double ymin = std::numeric_limits<double>::max();
				double ymax = std::numeric_limits<double>::min();

        for (auto v: epithelium.getCellCluster().getBoundaryVertices()){
            if (v->getPosition()[0] < xmin){
                xmin = v->getPosition()[0];
            }
            if (v->getPosition()[0] > xmax){
                xmax = v->getPosition()[0];
            }
            if (v->getPosition()[1] < ymin){
                ymin = v->getPosition()[1];
            }
            if (v->getPosition()[1] > ymax){
                ymax = v->getPosition()[1];
            }
        }

				std::cout << "xmin = " << xmin << '\n';
				std::cout << "xmax = " << xmax << '\n';
				std::cout << "diff = " << xmax-xmin << '\n';
				// std::cout << "h = " << (xmax-xmin)/10 << '\n';

		// Add elastic borders on the left, right, bottom and top.
		// -------------------------------------------------------
		ElasticLeftBorderEnvironment<2> elasticLeftBorder(elasticCoef,xmin);
		ElasticRightBorderEnvironment<2> elasticRightBorder(elasticCoef,xmax);
		ElasticBottomBorderEnvironment<2> elasticBottomBorder(elasticCoef,ymin);
		ElasticTopBorderEnvironment<2> elasticTopBorder(elasticCoef,ymax);

		GradientEnvironment<2> gradientEnvironment(activatedFunction, xmin, xmax, ymin, ymax);

		// Start applying external force on selected cell.
		// -----------------------------------------------
		T1Swap<2> t1(damping);

		std::vector<std::shared_ptr<Vertex<2>>> boundaryVertices;
		std::vector<std::shared_ptr<Vertex<2>>> allVertices;
		std::vector<std::shared_ptr<Vertex<2>>> innerVertices;


		// Identify the cell that will be constrained by a constant force.
		std::vector<std::shared_ptr<Cell<2>>> forcedCells;
		int forceCellIndex = 0.10*nrows*(ncols-0.5)+0.1*ncols; //0.10*nrows*(ncols-0.5)+0.5*ncols; //0.10*nrows*(ncols-0.5)+0.1*ncols; //0; //5*nrows + 5;
		// std::cout << "forceCellId = " << forceCellIndex << std::endl;
		std::shared_ptr<Cell<2>> forcedC = epithelium.getCellCluster().getCells()[forceCellIndex];
		// std::cout << "forced cell : " << *forcedC << '\n';
		forcedCells.push_back(forcedC);



		// // Force the cell cluster made up of the central cell and its neighbors.
		// for (auto n: epithelium.getCellCluster().getNeighborCells(forcedC)){
		// 	forcedCells.push_back(n);
		// }

		if (activatedFunction) {
			force = 0;
		}

		Array<double,2> cellForce(force,0.0); //cellForce(force,0.0); // cellForce(1.0e-7,1.0e-7);
		epithelium.addSpecialCellsAndForceFunction(
			forcedCells, std::bind(dataFunctionals2D::constForce,std::placeholders::_1,cellForce));

		// Add forced cell type to the cell cluster.
		epithelium.getCellCluster().getCellTypesManager().addCellType(type1);
		epithelium.getCellCluster().getCellTypesManager().setLambdaBetweenCellTypes(lambda,type1,type0);
		epithelium.getCellCluster().getCellTypesManager().setLambdaBetweenCellTypes(lambda,type1,type1);
		// Assign type1 to forced cells.
		for (auto c : forcedCells){
			c->setType(*type1);
		}

		// Create a text file containing all the parameter values used in this execution.
		std::string results_folder = workdir2;
		std::ofstream results_stream(results_folder.append("/results.txt").c_str());
		results_stream << "force\ttime_sec\tnbT1\tdegree\tposition_x\tposition_y" << std::endl;

		/*
		*/
		std::string density_folder = workdir2;
		std::ofstream density_stream(density_folder.append("/density.txt").c_str());
		density_stream << "time_sec\tnbT1\tarea_1\tarea_2\tarea_3\tarea_4\tarea_5\tarea_6\tarea_7\tarea_8\tarea_9\tarea_10" << std::endl;
		/*
		*/


		int nbT1 = 0;
		iT = 0;



		while(iT < MAXIT){// The simulation run until it reaches MAXIT iterations.
			// Important to consider the updated list of boundary vertices. It can change with transitions.
			boundaryVertices = epithelium.getCellCluster().getBoundaryVertices();
			allVertices = epithelium.getCellCluster().getVertices();
			innerVertices = epithelium.getCellCluster().getVertices();

			for (size_t i = 0; i < boundaryVertices.size(); i++) {
				for (size_t j = 0; j < innerVertices.size(); j++) {
					if (boundaryVertices[i] == innerVertices[j]) {
						innerVertices.erase(innerVertices.begin() + j);
					}
				}
			}

			epithelium.removeAllForcedVerticesAndForceFunction();

			epithelium.addSpecialVerticesAndForceFunction(
			boundaryVertices, std::bind(dataFunctionals2D::elasticForce,std::placeholders::_1,elasticLeftBorder) );
			epithelium.addSpecialVerticesAndForceFunction(
			boundaryVertices, std::bind(dataFunctionals2D::elasticForce,std::placeholders::_1,elasticRightBorder) );
			epithelium.addSpecialVerticesAndForceFunction(
			boundaryVertices, std::bind(dataFunctionals2D::elasticForce,std::placeholders::_1,elasticBottomBorder) );
			epithelium.addSpecialVerticesAndForceFunction(
			boundaryVertices, std::bind(dataFunctionals2D::elasticForce,std::placeholders::_1,elasticTopBorder) );

			epithelium.addSpecialVerticesAndForceFunction(
			innerVertices, std::bind(dataFunctionals2D::gradientForce,std::placeholders::_1, gradientEnvironment) );


			// Add here the fuction that adds forces on all vertices (or except the boundary ones).

			epithelium.performTimeStep(dt, false, false);

			if (iT % 1 == 0){ //if (iT % (int)(10/dt) == 0){
				bool done_t1 = true;
				while(done_t1){
					done_t1 = t1.check(epithelium.getCellCluster(),0.1*r0);
					if (done_t1){
						nbT1 ++;
						//std::cout << iT*dt << "\t" << nbT1 << "\t" << forcedCells[0]->computeDegree() << std::endl;
						results_stream << iT*dt << "\t" << nbT1 << "\t" << forcedCells[0]->computeDegree() << "\t" << forcedCells[0]->computeCentroid() << std::endl;
					}
				}
			}

			// VTK.
			// 1800/dt is number of .vtk files, or number od frames in vidéos.
			if (iT % (int)(100/dt) == 0){
	        	VTKwriter<2> vtk(epithelium.getCellCluster(), workdir2+"/tissue"+std::to_string(iT)+".vtk");
	        	vtk.addScalarToCell("id",dataFunctionals2D::getId);
	        	vtk.addScalarToCell("degree",dataFunctionals2D::computeDegree);
	        	vtk.addScalarToCell("onMitosis",dataFunctionals2D::isOnMitosis);
	        	vtk.addScalarToCell("type",dataFunctionals2D::getType);
 	        	vtk.addScalarToCell("age",dataFunctionals2D::computeAge);
	        	vtk.addScalarToCell("area",dataFunctionals2D::computeArea);
	        	vtk.addVectorToVertex("force",dataFunctionals2D::computeForce);
        	}

        	// Store time, nbT1, cell degree and position... each minute (60 sec)
        	if (iT % (int)(60/dt) == 0){
        		Array<double,2> position = forcedCells[0]->computeCentroid();
	        	results_stream << force << "\t" << iT*dt << "\t" << nbT1 << "\t" << forcedCells[0]->computeDegree() << "\t" << position[0] << "\t" << position[1] << std::endl;

						// "\ttime_sec\tnbT1\tarea_1\tarea_2\tarea_3\tarea_4\tarea_5\tarea_6\tarea_7\tarea_8\tarea_9\tarea_10"


						// Code to compute the density of cells along portion of the tissue
						if (activatedFunction) {

							double h = (xmax-xmin)/10;
							// Might be replaced by a vector
							density_1 = 0;
							density_2 = 0;
							density_3 = 0;
							density_4 = 0;
							density_5 = 0;
							density_6 = 0;
							density_7 = 0;
							density_8 = 0;
							density_9 = 0;
							density_10 = 0;

							for (auto cell: epithelium.getCellCluster().getCells()){
								Array<double,2> center = cell->computeCentroid();
									if (center[0] < xmin+h){
											density_1 = density_1 + 1;
									} else if (center[0] < xmin+2*h){
											density_2 = density_2 + 1;
									}	else if (center[0] < xmin+3*h){
											density_3 = density_3 + 1;
									}	else if (center[0] < xmin+4*h){
											density_4 = density_4 + 1;
									}	else if (center[0] < xmin+5*h){
											density_5 = density_5 + 1;
									}	else if (center[0] < xmin+6*h){
											density_6 = density_6 + 1;
									}	else if (center[0] < xmin+7*h){
											density_7 = density_7 + 1;
									}	else if (center[0] < xmin+8*h){
											density_8 = density_8 + 1;
									}	else if (center[0] < xmin+9*h){
											density_9 = density_9 + 1;
									}	else if (center[0] < xmin+10*h){
											density_10 = density_10 + 1;
									}
							}
						}
						density_stream << iT*dt << "\t\t" << nbT1 << "\t" << density_1  << "\t" << density_2 << "\t" << density_3 << "\t" << density_4 << "\t" << density_5 << "\t" << density_6 << "\t" << density_7 << "\t" << density_8 << "\t" << density_9 << "\t" << density_10 << std::endl;
        	}

					if (iT % (int)(10000) == 0){
						std::cout << "frame :" << iT/10000 << '\n';
					}

			iT++;

		}//end-while

		results_stream.close();
		density_stream.close();
		myseed++;
		myRep= myRep+1;

	}//end-while-myRep
	return 0;
}

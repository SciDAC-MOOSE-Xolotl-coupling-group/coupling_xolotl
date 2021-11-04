//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include <cmath>
#include "XolotlNetworkProblem.h"
#include "SystemBase.h"

using std::max;

registerMooseObject("coupling_xolotlApp", XolotlNetworkProblem);

template<>
InputParameters validParams<XolotlNetworkProblem>() {
	InputParameters params = validParams<ExternalProblem>();

	// Parameter for the Xolotl file name
	params.addRequiredParam < FileName
			> ("network_xolotl_filename", "Name with the path for the Xolotl input file with the full network");
	params.addRequiredParam < std::vector<FileName>
			> ("subnetwork_xolotl_filenames", "Name with the path for the Xolotl input files");
	return params;
}

XolotlNetworkProblem::XolotlNetworkProblem(const InputParameters &params) :
		ExternalProblem(params), _network_xolotl_filename(
				getParam < FileName > ("network_xolotl_filename")), _subnetwork_xolotl_filenames(
						getParam < std::vector<FileName>
								> ("subnetwork_xolotl_filenames")), _current_time(
				declareRestartableData < Real > ("current_time", 0.0)), _current_dt(
				declareRestartableData < Real > ("current_dt", 0.0)), _previous_time(
				declareRestartableData < Real > ("previous_time", 0.0)), _conc_vector(
				declareRestartableData
						< std::vector<
								std::vector<
										std::vector<
												std::vector<
														std::pair<
																xolotl::IdType,
																Real> > > > >
						> ("conc_vector")) {
	// Create the whole network interface
	int argc = 2;
	const char *argv[argc + 1];
	std::string fakeAppName = "mainXolotl";
	argv[0] = fakeAppName.c_str();
	argv[1] = _network_xolotl_filename.c_str();

	_networkInterface = std::make_shared<XolotlInterface>();
	_networkInterface->initializeXolotl(argc, argv,
			(_app.getCommunicator())->get());

	_subInterfaces.clear();
	// Loop on the number of parameter files
	for (auto name : _subnetwork_xolotl_filenames) {
		_subInterfaces.push_back(std::make_shared<XolotlInterface>());

		std::string fakeAppName = "subXolotl";
		argv[0] = fakeAppName.c_str();
		argv[1] = name.c_str();

		(_subInterfaces.back())->initializeXolotl(argc, argv, (_app.getCommunicator())->get());
	}

	// Exchange information about the sub networks
	std::vector < std::vector<std::vector<std::uint32_t>> > allBounds;
	// Loop on the sub interfaces
	for (auto inter : _subInterfaces) {
		// Get the bounds
		auto bounds = inter->getAllClusterBounds();

		// Add them to the main vector
		allBounds.push_back(bounds);
		_subDOFs.push_back(bounds.size());
	}

	// Pass it the the network instance
	_networkInterface->initializeClusterMaps(allBounds);

	// Take care of the fluxes
	auto fluxVector = _networkInterface->getImplantedFlux();
	for (auto i = 0; i < _subInterfaces.size(); i++) {
		_subInterfaces[i]->setImplantedFlux(fluxVector[i]);
	}
}

void XolotlNetworkProblem::externalSolve() {
	// Check that the next time is larger than the current one
	if (time() > _current_time) {
		std::vector < std::vector<double> > conc;
		// Loop on the sub interfaces to get all the concentrations
		for (auto i = 0; i < _subInterfaces.size(); i++) {
			auto sparseConc = _subInterfaces[i]->getConcVector();
			std::vector<double> subConc(_subDOFs[i], 0.0);
			for (auto pair : sparseConc[0][0][0]) {
				if (pair.first < _subDOFs[i])
					subConc[pair.first] = pair.second;
			}
			conc.push_back(subConc);
		}

		// Print the result
		_networkInterface->outputData(_current_time, conc);

		// Compute the new rates
		auto constantRates = _networkInterface->computeConstantRates(conc);

		// Pass them
		for (auto i = 0; i < _subInterfaces.size(); i++) {
			_subInterfaces[i]->setConstantRates(constantRates[i]);
			// Set the time we want to reach
			_subInterfaces[i]->setTimes(time(), dt());
			// Run the solver
			_subInterfaces[i]->solveXolotl();
		}
		// Save the current time
		_current_time = time();
	}
}

bool XolotlNetworkProblem::converged() {
	return true;
}

void XolotlNetworkProblem::saveState() {
	// Update the values from Xolotl
	_conc_vector = _networkInterface->getConcVector();
	_current_dt = _networkInterface->getCurrentDt();
	_previous_time = _networkInterface->getPreviousTime();

	xolotl::IdType i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
	_networkInterface->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
}

void XolotlNetworkProblem::setState() {
	// Set them in Xolotl
	_networkInterface->setConcVector(_conc_vector);
	_networkInterface->setCurrentTimes(_current_time, _current_dt);
	_networkInterface->setPreviousTime(_previous_time);
}

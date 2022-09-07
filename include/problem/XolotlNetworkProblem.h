//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XOLOTLNETWORKPROBLEM_H
#define XOLOTLNETWORKPROBLEM_H

#include "ExternalProblem.h"
#include "coupling_xolotlApp.h"

class XolotlNetworkProblem;

template<>
InputParameters validParams<XolotlNetworkProblem>();

/**
 * This is an interface to call an external solver
 */
class XolotlNetworkProblem: public ExternalProblem {
public:
	XolotlNetworkProblem(const InputParameters &params);

	~XolotlNetworkProblem() {
		std::vector<double> temperatures;
		std::vector<double> depths;
		_subInterfaces[0]->getNetworkTemperature(temperatures, depths);
		// Loop on the grid points
		std::vector < std::vector<std::vector<double> > > fullConc;
		// 0D
		if (temperatures.size() < 2) {
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
			fullConc.push_back(conc);
		}
		// 1D
		else {
			for (auto j = 0; j < temperatures.size() - 2; j++) {
				// Loop on the grid points
				for (auto j = 0; j < temperatures.size() - 2; j++) {
					std::vector < std::vector<double> > conc;
					// Loop on the sub interfaces to get all the concentrations
					for (auto i = 0; i < _subInterfaces.size(); i++) {
						auto sparseConc = _subInterfaces[i]->getConcVector();
						std::vector<double> subConc(_subDOFs[i], 0.0);
						for (auto pair : sparseConc[0][0][j]) {
							if (pair.first < _subDOFs[i]) {
								subConc[pair.first] = pair.second;
							}
						}
						conc.push_back(subConc);
					}

					fullConc.push_back(conc);
				}
			}
		}

		// Print the result
		_networkInterface->outputData(_current_time, fullConc,
				std::max((int) temperatures.size() - 2, 1));
	}

	virtual void externalSolve() override;
	virtual bool converged() override;
	virtual void syncSolutions(Direction /*direction*/) override {
		return;
	}

	// Methods for restart
	void saveState();
	void setState();

private:
	/// The path to the input file for Xolotl
	FileName _network_xolotl_filename;
	std::vector<FileName> _subnetwork_xolotl_filenames;
	std::shared_ptr<XolotlInterface> _networkInterface;
	std::vector<std::shared_ptr<XolotlInterface> > _subInterfaces;
	std::vector<xolotl::IdType> _subDOFs;
	Real &_current_time;
	Real _max_dt;
	xolotl::IdType _localXS;
	xolotl::IdType _localXM;

	// Variables for restart
	Real &_current_dt;
	Real &_previous_time;
	std::vector<
			std::vector<
					std::vector<std::vector<std::pair<xolotl::IdType, Real> > > > > &_conc_vector;

}
;

#endif /* XOLOTLNETWORKPROBLEM_H */

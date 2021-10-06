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
		_networkInterface->finalizeXolotl();
	}

	virtual void externalSolve() override;
	virtual bool converged() override;
	virtual void syncSolutions(Direction /*direction*/) override {return;}

	// Methods for restart
	void saveState();
	void setState();

private:
	/// The path to the input file for Xolotl
	FileName _network_xolotl_filename;
	std::shared_ptr<XolotlInterface> _networkInterface;
	std::vector<std::shared_ptr<XolotlInterface> > _subInterfaces;
	std::vector<xolotl::IdType> _subDOFs;
	Real &_current_time;

	// Variables for restart
	Real &_current_dt;
	Real &_previous_time;
	std::vector<std::vector<std::vector<std::vector<std::pair<xolotl::IdType, Real> > > > > &_conc_vector;

};

#endif /* XOLOTLNETWORKPROBLEM_H */

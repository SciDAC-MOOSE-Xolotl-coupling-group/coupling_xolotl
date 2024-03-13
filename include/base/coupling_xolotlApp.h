//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

#include <memory>

#include "MooseApp.h"
#include <xolotl/interface/XolotlInterface.h>

using XolotlInterface = xolotl::interface::XolotlInterface;

class coupling_xolotlApp: public MooseApp {
public:
	static InputParameters validParams();

	coupling_xolotlApp(InputParameters parameters);
	virtual ~coupling_xolotlApp();

	void createInterface(FileName paramName);

	std::shared_ptr<XolotlInterface> getInterface() {
		return _interface;
	}
	TS& getXolotlTS() {
		return _interface->getTS();
	}
	static void registerApps();
	static void registerAll(Factory &f, ActionFactory &af, Syntax &s);

	// For restart capabilities
	std::shared_ptr<Backup> backup();
	virtual void preBackup() override;
	virtual void postRestore(bool for_restart = false) override;

private:
	std::shared_ptr<XolotlInterface> _interface;
	bool _is_xolotl_app;
};

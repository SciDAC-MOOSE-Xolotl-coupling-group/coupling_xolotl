# default length unit: nm
# default time unit: s
# default mass unit: ?

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
[]

[Problem]
 type = XolotlNetworkProblem
 network_xolotl_filename = 'param_full_network.txt'
 subnetwork_xolotl_filenames = 'param_B.txt param_V.txt param_I.txt'
 max_dt = 1
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e-10
    growth_factor = 1.1
    cutback_factor = 0.1
  [../]
#  [./TimeStepper]
#    type = ConstantDT
#    dt = 1.0
#  [../]
  start_time = 0
  end_time = 1.0e3
  num_steps = 10
[]

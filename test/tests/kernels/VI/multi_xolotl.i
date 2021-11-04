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
 subnetwork_xolotl_filenames = 'param_V.txt param_I.txt'
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = ConstantDT
    dt = 0.01
  [../]
  start_time = 0
  end_time = 1.0
  num_steps = 5  
[]

# default length unit: nm
# default time unit: s
# default mass unit: ?
[Mesh]
  type = XolotlReflectedMesh
  dim = 1
  subnetwork_xolotl_filenames = 'param_V.txt'
[]

[Problem]
 type = XolotlNetworkProblem
 network_xolotl_filename = 'param_full_network.txt'
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = ConstantDT
    dt = 0.01
  [../]
  start_time = 0
  end_time = 1.0
  num_steps = 1  
[]

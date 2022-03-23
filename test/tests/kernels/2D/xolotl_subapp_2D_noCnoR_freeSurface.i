# default length unit: nm
# default time unit: s
# default mass unit: ?
[Mesh]
  type = XolotlReflectedMesh
  dim = 2
  XolotlInput_path_name = './param_2D_noCnoR_freeSurface.txt'
[]

[AuxVariables]
  [./Auxv]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxGB]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxMono]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxFrac]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxPid]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [pid_aux]
    type = ProcessorIDAux
    variable = AuxPid
    execute_on = 'INITIAL'
  []
[]

[Problem]
 type = XolotlProblem
 sync_rate = Auxv
 sync_GB = AuxGB
 sync_mono = AuxMono
 sync_frac = AuxFrac
 free_surface = no
[]

[Variables]
  [./d]
  [../]
[]
[ICs]
  [./Init_Aux_const]
    type = ConstantIC
    variable = Auxv
    value = 0
  [../]
  [./Init_Aux_gb]
    type = ConstantIC
    variable = AuxGB
    value = 1
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = ConstantDT
    dt = 5000.0
  [../]
  start_time = 0
  end_time = 1.2e8
  # end_time = 20000.0
[]

[Outputs]
  exodus = true
[]

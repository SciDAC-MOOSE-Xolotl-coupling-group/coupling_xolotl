# default length unit: nm
# default time unit: s
# default mass unit: ?
[Mesh]
  type = XolotlReflectedMesh
  dim = 2
  XolotlInput_path_name = './param_2D_W.json'
[]

[AuxVariables]
  [./AuxHRate]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxGB]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxVRate]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxHConc]
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
 sync_H_rate = AuxHRate
 sync_V_rate = AuxVRate
 sync_H = AuxHConc
 sync_GB = AuxGB
 free_surface = no
[]

[Variables]
  [./d]
  [../]
[]
[ICs]
  [./Init_AuxH_const]
    type = ConstantIC
    variable = AuxHRate
    value = 0
  [../]
  [./Init_AuxV_const]
    type = ConstantIC
    variable = AuxVRate
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
    dt = 1.0e-4
  [../]
  start_time = 0
  end_time = 1.2e8
  # end_time = 20000.0
[]

[Outputs]
  exodus = true
[]

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 32
    ny = 32
    xmin = 0
    xmax = 20000
    ymin = 0
    ymax = 20000
    partition = square
  []
  parallel_type = distributed
[]

[GlobalParams]
  op_num = 11
  var_name_base = etam
  grain_num = 32
  int_width = 480
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt
    # coloring_algorithm = bt # We must use bt to force the UserObject to assign one grain to each op
    rand_seed = 5
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
  [../]
  [./term]
    type = Terminator
    expression = 'grain_tracker < 15'
  [../]
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[AuxVariables]
  [./bubble]
    order = FIRST
    family = LAGRANGE
  [../]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./time]
  [../]
  [./unique_grains]
    order = FIRST
    family = MONOMIAL
  [../]
  [./var_indices]
    order = FIRST
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./time]
    type = FunctionAux
    variable = time
    function = 't'
  [../]

  #additional_for_grain_tracker
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
  [./bnds]
    type = ConstantIC
    variable = bnds
    value = 1
  [../]
  [./bubble_IC_no_bubble]
    variable = bubble
    type = ConstantIC
    value = 0
  [../]
[]

[BCs]
  [./Periodic]
    [./x_y_direction_periodic]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    T = 500 # K
    wGB = 480 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./num_DOF]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./elapsed]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./num_nonlinear_residuals]
    type = NumNonlinearIterations
  [../]
  [./num_linear_residuals]
    type = NumLinearIterations
  [../]
  [./numgrain]
    type = ElementExtremeValue
    variable = unique_grains
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -pc_factor_levels'
  petsc_options_value = 'asm      lu              1             2'
  nl_max_its = 15
  l_max_its = 15
  l_tol = 1.0e-3
  # nl_max_its = 20
  # l_max_its = 40
  # l_tol = 1e-04
  nl_rel_tol = 1.0e-8
  num_steps = 500
  end_time = 1e8
  nl_abs_tol = 1e-9
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 5
    growth_factor = 1.2
    cutback_factor = 0.8
    # dt = 0.5
    dt = 1
    # adapt_log = true
  [../]
[]

[Outputs]
  [./nemesis]
    type = Nemesis
  [../]
  perf_graph = true
  csv = true
[]


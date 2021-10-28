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
  grain_num = 14
  int_width = 480
[]

[UserObjects]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
  [../]
  [./term]
    type = Terminator
    expression = 'grain_tracker < 11'
  [../]
  [./poly]
    type = Polycrystal
    bnds_var = 'bnds'
    num_bub = 5
    radius_bub = 480
    rand_seed = 10
    execute_on = timestep_begin
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
  # For use of GrainTracker
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
  [./bubble_IC]
    variable = bubble
    type = PolycrystalGBBubbleIC
    invalue = 1.0
    outvalue = 0.0
    poly_ic_uo = poly
    bubspac = 0
    numbub = 5
    radius = 400
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
    c = bubble
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
  [./bubble_fraction_new]
    type = ElementAverageValue
    variable = bubble
    # outputs = csv
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
  num_steps = 1000000
  nl_rel_tol = 1.0e-8
  # end_time = 1e8
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

[MultiApps]
  [./sub]
    type = FullSolveMultiApp
    positions = '0. 0. 0.'
    input_files = step1_parallel.i
    execute_on = 'initial'
    output_in_position = true
  [../]
[]

[Transfers]
  [./from_sub0]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam0
    variable = etam0
    multi_app = sub
  [../]
  [./from_sub1]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam1
    variable = etam1
    multi_app = sub
  [../]
  [./from_sub2]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam2
    variable = etam2
    multi_app = sub
  [../]
  [./from_sub3]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam3
    variable = etam3
    multi_app = sub
  [../]
  [./from_sub4]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam4
    variable = etam4
    multi_app = sub
  [../]
  [./from_sub5]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam5
    variable = etam5
    multi_app = sub
  [../]
  [./from_sub6]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam6
    variable = etam6
    multi_app = sub
  [../]
  [./from_sub7]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam7
    variable = etam7
    multi_app = sub
  [../]
  [./from_sub8]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam8
    variable = etam8
    multi_app = sub
  [../]
  [./from_sub9]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam9
    variable = etam9
    multi_app = sub
  [../]
  [./from_sub10]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = etam10
    variable = etam10
    multi_app = sub
  [../]
  [./from_sub11]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = bnds
    variable = bnds
    multi_app = sub
  [../]
[]

[Outputs]
  [./nemesis]
    type = Nemesis
  [../]
  perf_graph = true
  csv = true
[]

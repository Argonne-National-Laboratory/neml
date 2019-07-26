# Simple 3D test

[GlobalParams]
  displacements = "disp_x disp_y disp_z"
[]

[Mesh]
  type = FileMesh
  file = "small.exo"
[]

[NEMLMechanics]
  kinematics = large
  add_all_output = true
[]

[BCs]
  [./leftx]
    type = PresetBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./lefty]
    type = PresetBC
    boundary = back
    variable = disp_y
    value = 0.0
  [../]
  [./leftz]
    type = PresetBC
    boundary = bottom
    variable = disp_z
    value = 0.0
  [../]
  [./pull_z]
    type = FunctionPresetBC
    boundary = top
    variable = disp_z
    function = pfn
  [../]
[]

[Functions]
  [./pfn]
    type = PiecewiseLinear
    x = '0    10'
    y = '0.00 0.1'
  [../]
[]

[Materials]
  [./stress1]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_1"
    large_kinematics = true
    block = 1
  [../]
  [./stress2]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_2"
    large_kinematics = true
    block = 2
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  
  dt = 1
  end_time = 10.0
[]

[Outputs]
  exodus = true
[]

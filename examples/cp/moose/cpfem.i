[Mesh]
  type = FileMesh
  file = 'n10.exo'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
      [./disp_x]
            order = second
      [../]
      [./disp_y]
            order = second
      [../]
      [./disp_z]
            order = second
      [../]
[]

[Functions]
  [./pfn]
    type = PiecewiseLinear
    x = '0    100'
    y = '0.00 0.01'
  [../]
[]

[Kernels]
  [./sdx]
      type = StressDivergenceNEML
      variable = disp_x
      component = 0
      use_displaced_mesh = true
  [../]
  [./sdy]
      type = StressDivergenceNEML
      variable = disp_y
      component = 1
      use_displaced_mesh = true
  [../]
  [./sdz]
      type = StressDivergenceNEML
      variable = disp_z
      component = 2
      use_displaced_mesh = true
  [../]
[]

[BCs]
  [./left]
     type = PresetBC
     variable = disp_x
     boundary = left
     value = 0.0
  [../]

  [./bottom]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./back]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./front]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = pfn
  [../]
[]

[Materials]
  [./strain]
    type = ComputeNEMLStrain
    large_kinematics = true
  [../]
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
  [./stress3]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_3"
    large_kinematics = true
    block = 3
  [../]
  [./stress4]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_4"
    large_kinematics = true
    block = 4
  [../]
  [./stress5]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_5"
    large_kinematics = true
    block = 5
  [../]
  [./stress6]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_6"
    large_kinematics = true
    block = 6
  [../]
  [./stress7]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_7"
    large_kinematics = true
    block = 7
  [../]
  [./stress8]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_8"
    large_kinematics = true
    block = 8
  [../]
  [./stress9]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_9"
    large_kinematics = true
    block = 9
  [../]
  [./stress10]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "grain_10"
    large_kinematics = true
    block = 10
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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  nl_max_its = 20

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8

  end_time = 100.0
  dtmin = 0.1
  dt = 1.0
[]

[Outputs]
  exodus = true
[]

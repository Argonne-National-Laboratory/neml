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
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "elastic_model"
    large_kinematics = true
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

  petsc_options_iname = -pc_type
  petsc_options_value = lu

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8

  end_time = 100.0
  dtmin = 1.0
  dt = 1.0
[]

[Outputs]
  exodus = true
[]

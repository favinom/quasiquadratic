[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 40
 ny = 4
 xmin = 0.0
 xmax = 10.0
 ymin = 0.0
 ymax = 1.0
 elem_type = TRI6
[]

 [Functions]
 [./parsed_function]
    type = ParsedFunction
    value = '0.0'
 [../]
 []
 
# [ICs]
# [./u_ic]
# type = FunctionIC
# variable = 'disp_x'
# function = parsed_function
# [../]
# []

 
[Variables]
 [./disp_x]  order = SECOND  family=LAGRANGE
 
 [./InitialCondition]
 type = BoundingBoxIC
 x1 = -0.0001
 x2 =  0.0001
 y1 = -0.0001
 y2 =  0.0001
 inside =   0.0
 outside =  0.0
 [../]

 
 [../]
 [./disp_y]  order = SECOND  family=LAGRANGE [../]
 []
 
[Kernels]
[./QQElementKernel_x]
 type = NeohookeanQQ_slow
 component =0
 variable = disp_x
 disp_x = disp_x
 disp_y = disp_y
[../]

 [./QQElementKernel_y]
 type = NeohookeanQQ_slow
 component =1
 variable = disp_y
 disp_x = disp_x
 disp_y = disp_y
 [../]

 []
 
 [BCs]
#active = 'bottom_y'
 [./left_x]
 type = DirichletBC
 variable = disp_x
 boundary = left
 value = 0
 [../]

 [./left_y]
 type = DirichletBC
 variable = disp_y
 boundary = left
 value = 0
 [../]

 [./bottom_y]
 type = NeumannBC
 variable = disp_y
 boundary = bottom
 value = 0.0001
 [../]
 
 []
 
 [Preconditioning]
 [./SMP]
type = FDP
#type = PBP
#type = SMP
 full = true
 [../]
 []
 
[Executioner]
 type=Steady
 solve_type=NEWTON

line_search = 'none'
nl_abs_tol = 1e-8 
#petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
# petsc_options_value='   preonly   lu       NONZERO               mumps                '

 petsc_options_iname=' -ksp_type -pc_type' #  -mat_view '
 petsc_options_value='   preonly   lu    ' #::ascii_matlab       '
 
 [./Quadrature]
 type = SIMPSON
 [../]

 
 []

 [Problem]
 type = FEProblem
 solve = true
 kernel_coverage_check = true
 []
 
 
[Outputs]
 file_base = incompressible
 exodus = true
 print_linear_residuals = true
 print_perf_log = true
 []

[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 1
 ny = 1
 xmin = 0.0
 xmax = 1.0
 ymin = 0.0
 ymax = 1.0
 elem_type = TRI6
[]
 
[Functions]
 [./parsed_function]
 type = ParsedFunction
 value = 'x'
 [../]
[]

[Variables]
 [./disp_x]  order = SECOND  family=LAGRANGE [../]
 [./disp_y]  order = SECOND  family=LAGRANGE [../]
 []

 [ICs]
 [./u_ic]
 type = FunctionIC
 variable = 'disp_x'
 function = parsed_function
 [../]
 []
 
[Kernels]
[./QQElementKernel]
 type = QQElement
 variable = disp_x
 disp_x = disp_x
 disp_y = disp_y
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
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='   preonly   lu       NONZERO               mumps                '
 
 [./Quadrature]
 type = SIMPSON
 [../]

 
 []

 [Problem]
 type = FEProblem
 solve = true
 kernel_coverage_check = false
 []
 
 
[Outputs]
 file_base = incompressible
 exodus = true
 print_linear_residuals = true
 print_perf_log = true
 []

[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 2
 ny = 2
 xmin = 0.0
 xmax = 1.0
 ymin = 0.0
 ymax = 1.0
 elem_type = TRI6
[]

[Variables]
 [./disp_x]  order = SECOND  family=LAGRANGE [../]
 [./disp_y]  order = SECOND  family=LAGRANGE [../]
 []
 
[Kernels]
[./QQElementKernel]
 type = QQElement
 variable = disp_x
 disp_x = disp_x
 disp_y = disp_y
 [../]
 []
 
 [BCs]
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

 [./bottom_x]
 type = NeumannBC
 variable = disp_y
 boundary = bottom
 value = -0.000001
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

petsc_options_iname=' -ksp_type -pc_type  -mat_view '
 petsc_options_value='   preonly   lu    ::ascii_matlab       '


 
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

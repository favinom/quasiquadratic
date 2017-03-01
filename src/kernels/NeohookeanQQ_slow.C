#include "NeohookeanQQ_slow.h"

#include "MooseMesh.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

template <>
InputParameters validParams<NeohookeanQQ_slow>() {
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "");
  params.addRequiredCoupledVar("disp_y", "");
  //params.addCoupledVar        ("disp_z", "");

 params.addRequiredParam<unsigned>("component", "component");
    
  return params;
}

NeohookeanQQ_slow::NeohookeanQQ_slow(const InputParameters & params) :
    Kernel(params),
    _component(getParam<unsigned>("component")),
    _disp_x_var(coupled("disp_x")),
    _disp_y_var(coupled("disp_y")),
    //_disp_z_var(_mesh.dimension() == 3 ? coupled("disp_z") : 100000),
    _disp_x(coupledValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    //_disp_z(coupledValue("disp_z")),
    _grad_disp_x(coupledGradient("disp_x")),
    _grad_disp_y(coupledGradient("disp_y"))/*,
    _grad_disp_z(coupledGradient("disp_z"))*/
{

    if (_mesh.dimension() == 3)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);

    if (_mesh.dimension() == 2)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
    
    // here we will store the gradient of basis functions: 0,1,2 are the fine ones, 3 the coarse ones
    _qq_grad_test=new RealVectorValue * [4];
    for (int i=0 ; i < 4; ++i)
        _qq_grad_test[i]=new RealVectorValue [3];
    
    // here we allocate mechanical quantities per element
    U = new RealTensorValue [4];
    F = new RealTensorValue [4];
    C = new RealTensorValue [4];
    S = new RealTensorValue [4];

    // EQQ is per node coarse
    CQQ = new RealTensorValue [3];
    
    SQQ = new RealTensorValue [3];
    

    // eps_lin, per dimension, per triangle, per node
    _eps_lin = new RealTensorValue ** [2];
    for (int i=0; i<2; ++i)
        _eps_lin[i]=new RealTensorValue * [4];
 
    for (int i=0; i<2; ++i)
        for (int j=0; j<4; ++j)
            _eps_lin[i][j]=new RealTensorValue[6];
    
    // _eps_lin_QQ per dimension, per node fine (dof), per node coarse
    _eps_lin_QQ = new RealTensorValue ** [2];
    for (int i=0; i<2; ++i)
        _eps_lin_QQ[i]=new RealTensorValue * [6];
    
    for (int dim=0; dim<2; ++dim)
       for (int j=0; j<6; ++j)
           _eps_lin_QQ[dim][j]= new RealTensorValue [3];
    
    
    _local_to_global= new int *[4];
    for (int i=0; i<4;++i)
        _local_to_global[i]=new int [3];
    
    _local_to_global[0][0]=0;
    _local_to_global[0][1]=3;
    _local_to_global[0][2]=4;
    _local_to_global[1][0]=1;
    _local_to_global[1][1]=5;
    _local_to_global[1][2]=3;
    _local_to_global[2][0]=2;
    _local_to_global[2][1]=4;
    _local_to_global[2][2]=5;
    _local_to_global[3][0]=0;
    _local_to_global[3][1]=1;
    _local_to_global[3][2]=2;
    
    simpson_to_tri6= new int [6];
    
    simpson_to_tri6[0]=0;
    simpson_to_tri6[1]=1;
    simpson_to_tri6[2]=2;
    simpson_to_tri6[3]=3;
    simpson_to_tri6[4]=5;
    simpson_to_tri6[5]=4;
    
    _mu=1.0;
    _lambda=3.0;
    

}

void NeohookeanQQ_slow::computeResidual()
{
    if (_mesh.dimension()==2)
        computeResidual2D();
}

void NeohookeanQQ_slow::computeResidual2D()
{
    Real area=0.0;
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        area+=_JxW[_qp];

    DenseMatrix<Number> mass_matrix;
    mass_matrix.resize(3,3);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            if (i==j)
                mass_matrix(i,j)=area/6.0;
            else
                mass_matrix(i,j)=area/12.0;

    mass_matrix.resize(3,3);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            if (i==j)
                mass_matrix(i,j)=area/3.0;
            else
                mass_matrix(i,j)=0.0;
    
    
    for (int tri=0; tri<4; ++tri)
        computeGradient(_q_point[_local_to_global[tri][0]],
                        _q_point[_local_to_global[tri][1]],
                        _q_point[_local_to_global[tri][2]],_qq_grad_test[tri]);
    
    
    RealVectorValue Temp1,Temp2,Temp3;
    
    // Strain
    
    for (int tri=0; tri<4; ++tri) {

        Temp1 = _qq_grad_test[tri][0]*_disp_x[_local_to_global[tri][0]]+
                _qq_grad_test[tri][1]*_disp_x[_local_to_global[tri][1]]+
                _qq_grad_test[tri][2]*_disp_x[_local_to_global[tri][2]];
        Temp2 = _qq_grad_test[tri][0]*_disp_y[_local_to_global[tri][0]]+
                _qq_grad_test[tri][1]*_disp_y[_local_to_global[tri][1]]+
                _qq_grad_test[tri][2]*_disp_y[_local_to_global[tri][2]];
        Temp3 = RealVectorValue(0.0,0.0,0.0);
    
        U[tri]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));
    }
    
    for (int i=0; i<4; ++i)
    {
        F[i]=U[i]+_identity;
        C[i]=F[i].transpose()*F[i];
        
        // CALCOLO INVERSA
        C[i](2,2)=1.0;
        RealTensorValue invC=C[i].inverse();
        C[i](2,2)=0.0;
        
        invC(2,2)=0.0;
        S[i]=_identity-invC;
    }
    
    for (int nodo_coarse=0; nodo_coarse<3; ++nodo_coarse)
    {
        CQQ[nodo_coarse]=2.0*C[nodo_coarse]-C[3];
        CQQ[nodo_coarse](2,2)=1.0;
        
        RealTensorValue invCQQ;
        invCQQ=CQQ[nodo_coarse].inverse();
        invCQQ(2,2)=0.0;
        CQQ[nodo_coarse](2,2)=0.0;
        SQQ[nodo_coarse]=_identity-invCQQ;
//        SQQ[nodo_coarse]=2.0*S[nodo_coarse]-S[3];
    }
    // FINE ASSEMBLING STRAIN QQ

    for (int dim=0; dim<2; ++dim)
        for (int tri=0; tri<4; ++tri)
            assembleStrainLin(dim,tri,_local_to_global[tri]);
    
    for (int dim=0; dim<2; ++dim)
        for (int nodo=0; nodo<6; ++nodo)
            for (int nodo_coarse=0; nodo_coarse<3; ++nodo_coarse)
                _eps_lin_QQ[dim][nodo][nodo_coarse]=2.0*_eps_lin[dim][nodo_coarse][nodo]-_eps_lin[dim][3][nodo];
    
    
    DenseVector<Number> & f_x = _assembly.residualBlock(_disp_x_var);
    DenseVector<Number> & f_y = _assembly.residualBlock(_disp_y_var);
    DenseVector<Number> _f_local[2];
    
    for (int dim=0; dim<2; ++dim)
    {
        _f_local[dim].resize(f_x.size());
        _f_local[dim].zero();
    }
 
    
    DenseVector<Number> V(3);
    DenseVector<Number> VLin(3);
    DenseVector<Number> res(3);
    
    for (int dim=0; dim<2; ++dim)
        for (int nodo = 0; nodo < 6; ++nodo)
        {
            _f_local[dim](simpson_to_tri6[nodo])=0.0;
            
            
            for (int i_local=0; i_local<2; ++i_local)
                for (int j_local=0; j_local<2; ++j_local)
                {
                    
                    for (int nodo_coarse = 0; nodo_coarse < 3; ++nodo_coarse)  // nuova quadratura
                    {
                        V(nodo_coarse)=SQQ[nodo_coarse](i_local,j_local);
                        VLin(nodo_coarse)=_eps_lin_QQ[dim][nodo][nodo_coarse](i_local,j_local);
                    }
                    
                    mass_matrix.vector_mult(res,V);
                    _f_local[dim](simpson_to_tri6[nodo]) += res.dot(VLin);
                }
            
        }
    
    
    
    if(_component==0)
        f_x += _f_local[0];
    if(_component==1)
        f_y += _f_local[1];
    
//    std::cout<<_f_local[0]<<std::endl;
//    std::cout<<_f_local[1]<<std::endl;
//    
//    for (int i=0; i<6; ++i)
//        std::cout<<_q_point[i]<<" "<<_disp_x[i]<<" "<<_u[i]<<std::endl;
    
//    std::cout<<std::endl<<std::endl;

//    int i;
//    std::cin>>i;
    
    //exit(1);
    
}

//void QQElement::computeOffDiagJacobian(unsigned int jvar)
//{
//    if (jvar==_var.number())
//        computeJacobian2D();
//
//    
//    if(_has_diag_save_in)
//    {
//        mooseError("Error: diag in already saved in");
//    }
//    
//    
//}

/*void QQElement::computeJacobian2D()
{
 
//    // References to the right components of the stiffness matrix
//    DenseMatrix<Number> & A_R_xx = _assembly.jacobianBlock(_disp_real_x_var, _disp_real_x_var);
//    DenseMatrix<Number> & A_R_xy = _assembly.jacobianBlock(_disp_real_x_var, _disp_real_y_var);
//    DenseMatrix<Number> & B_R_xp = _assembly.jacobianBlock(_disp_real_x_var, _p_real_var);
//    
//    DenseMatrix<Number> & A_R_yx = _assembly.jacobianBlock(_disp_real_y_var, _disp_real_x_var);
//    DenseMatrix<Number> & A_R_yy = _assembly.jacobianBlock(_disp_real_y_var, _disp_real_y_var);
//    DenseMatrix<Number> & B_R_yp = _assembly.jacobianBlock(_disp_real_y_var, _p_real_var);
//    
//    DenseMatrix<Number> & A_I_xx = _assembly.jacobianBlock(_disp_imag_x_var, _disp_imag_x_var);
//    DenseMatrix<Number> & A_I_xy = _assembly.jacobianBlock(_disp_imag_x_var, _disp_imag_y_var);
//    DenseMatrix<Number> & B_I_xp = _assembly.jacobianBlock(_disp_imag_x_var, _p_imag_var);
//    
//    DenseMatrix<Number> & A_I_yx = _assembly.jacobianBlock(_disp_imag_y_var, _disp_imag_x_var);
//    DenseMatrix<Number> & A_I_yy = _assembly.jacobianBlock(_disp_imag_y_var, _disp_imag_y_var);
//    DenseMatrix<Number> & B_I_yp = _assembly.jacobianBlock(_disp_imag_y_var, _p_imag_var);
//    
//    DenseMatrix<Number> & C_RR = _assembly.jacobianBlock(_p_real_var, _p_real_var);
//    DenseMatrix<Number> & M_RI = _assembly.jacobianBlock(_p_real_var, _p_imag_var);
//    DenseMatrix<Number> & M_IR = _assembly.jacobianBlock(_p_imag_var, _p_real_var);
//    DenseMatrix<Number> & C_II = _assembly.jacobianBlock(_p_imag_var, _p_imag_var);
//    
//    DenseMatrix<Number> & B_pR_xI = _assembly.jacobianBlock(_p_real_var,_disp_imag_x_var);
//    DenseMatrix<Number> & B_pR_yI = _assembly.jacobianBlock(_p_real_var,_disp_imag_y_var);
//    
//    DenseMatrix<Number> & B_pI_xR = _assembly.jacobianBlock(_p_imag_var,_disp_real_x_var);
//    DenseMatrix<Number> & B_pI_yR = _assembly.jacobianBlock(_p_imag_var,_disp_real_y_var);
//    
//    DenseMatrix<Number> _A_local[2][2];
//    DenseMatrix<Number> _BT_local[2];
//    DenseMatrix<Number> _C_local;
//    DenseMatrix<Number> _M_local;
//    
//    for (int i=0; i<2; ++i)
//    {
//        for (int j=0; j<2; ++j)
//        {
//            _A_local[i][j].resize(A_R_xx.m(),A_R_xx.n());
//            _A_local[i][j].zero();
//        }
//        _BT_local[i].resize(B_R_xp.m(),B_R_xp.n());
//        _BT_local[i].zero();
//    }
//    _C_local.resize(C_RR.m(),C_RR.n());
//    _C_local.zero();
//    _M_local.resize(M_RI.m(),M_RI.n());
//    _M_local.zero();
//    
//    initTensorVariables();
//    
//    for (int idim=0; idim<2; ++idim)
//        for (int jdim=idim; jdim<2; ++jdim)
//            for (_i = 0; _i < _P2.size(); ++_i)
//                for (_j = 0; _j < _P2.size(); ++_j)
//                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                    {
//                        RealTensorValue & sigma=_sigma_lin[jdim][_j][_qp];
//                        RealTensorValue & V=_V[idim][_i][_qp];
//                        _A_local[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*sigma.contract(V);
//                    }
//    
//    for (int idim=0; idim<2; ++idim)
//        for (_i = 0; _i < _P2.size(); ++_i)
//            for (_j = 0; _j < _P1.size(); ++_j)
//                for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                {
//                    RealTensorValue & V=_V[idim][_i][_qp];
//                    _BT_local[idim](_i,_j)+=_JxW[_qp] * _coord[_qp]*_alpha[_qp]*_P1[_j][_qp]*V.tr();
//                }
//    
//    DenseMatrix<Number> _B_local[2];
//    _BT_local[0].get_transpose(_B_local[0]);
//    _BT_local[1].get_transpose(_B_local[1]);
//
//
//    
//    
//    for (int i=0; i<_A_local[1][0].m(); ++i)
//        for (int j=0; j<_A_local[1][0].m(); ++j)
//            _A_local[1][0](i,j)=_A_local[0][1](j,i);
//    
//    for (_i=0;_i<_P1.size(); ++_i)
//        for (_j=0;_j<_P1.size(); ++_j)
//            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                _C_local(_i,_j)+=_JxW[_qp] * _coord[_qp]*_diffusion[_qp]*_grad_P1[_i][_qp]*_grad_P1[_j][_qp];
//
//    for (_i=0;_i<_P1.size(); ++_i)
//        for (_j=0;_j<_P1.size(); ++_j)
//            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                _M_local(_i,_j)+=_JxW[_qp] * _coord[_qp]*_inv_m[_qp]*_P1[_i][_qp]*_P1[_j][_qp];
//
//    
//    
//    A_R_xx += _A_local[0][0];
//    A_R_xy += _A_local[0][1];
//    A_R_yx += _A_local[1][0];
//    A_R_yy += _A_local[1][1];
//    
//    A_I_xx += _A_local[0][0];
//    A_I_xy += _A_local[0][1];
//    A_I_yx += _A_local[1][0];
//    A_I_yy += _A_local[1][1];
//    
//    B_R_xp -= _BT_local[0];
//    B_R_yp -= _BT_local[1];
//    
//    B_I_xp -= _BT_local[0];
//    B_I_yp -= _BT_local[1];
//    
//    C_RR   -= _C_local;
//    C_II   -= _C_local;
//
//    M_RI   += _M_local;
//    M_IR   -= _M_local;
//    
//    B_pR_xI+=_B_local[0];
//    B_pR_yI+=_B_local[1];
//    B_pI_xR-=_B_local[0];
//    B_pI_yR-=_B_local[1];
    
}*/

/*void QQElement::initTensorVariables()
{
//        for (int d=0; d<_mesh.dimension(); ++d)
//            for ( _j = 0; _j < _P2.size() ; ++_j)
//                for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                {
//                    _V[d][_j][_qp]=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
//                    for (int j=0; j<_mesh.dimension(); ++j)
//                        _V[d][_j][_qp](d,j)=_grad_P2[_j][_qp](j);
//                    
//                    RealTensorValue & V=_V[d][_j][_qp];
//                    _eps_lin[d][_j][_qp]=0.5*( V+V.transpose() );
//                    
//                    RealTensorValue & eps=_eps_lin[d][_j][_qp];
//                    _sigma_lin[d][_j][_qp]=2.0 * _mu[_qp]*eps + _lambda[_qp]*eps.tr()*_identity;
//                }
//    
}*/

void NeohookeanQQ_slow::computeGradient(RealVectorValue x0, RealVectorValue x1, RealVectorValue x2, RealVectorValue * Gradient)
{
    RealTensorValue coordinates;
    coordinates(0,0)=x0(0);
    coordinates(0,1)=x0(1);
    coordinates(0,2)=1.0;
    coordinates(1,0)=x1(0);
    coordinates(1,1)=x1(1);
    coordinates(1,2)=1.0;
    coordinates(2,0)=x2(0);
    coordinates(2,1)=x2(1);
    coordinates(2,2)=1.0;
    
    coordinates=coordinates.inverse();//RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    Gradient[0]=RealVectorValue(coordinates(0,0),coordinates(1,0),0.0);
    Gradient[1]=RealVectorValue(coordinates(0,1),coordinates(1,1),0.0);
    Gradient[2]=RealVectorValue(coordinates(0,2),coordinates(1,2),0.0);
}

void NeohookeanQQ_slow::assembleStrainLin(int dim, int triangle, int * _local_to_global)
{
    // set to zero strain_lin
    for (int i=0; i <6; ++i)   // i is the quadrature node
        _eps_lin[dim][triangle][i]=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    for (int i=0; i<3; ++i) // iteration over the nodes
    {
        _V=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        for (int k=0; k<3; ++k)
            _V(dim,k)=_qq_grad_test[triangle][i](k);
        
        _eps_lin[dim][triangle][_local_to_global[i]]=0.5*(F[triangle].transpose()*_V+_V.transpose()*F[triangle]);
    }

}




























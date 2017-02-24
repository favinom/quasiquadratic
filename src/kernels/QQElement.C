#include "QQElement.h"

#include "MooseMesh.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

template <>
InputParameters validParams<QQElement>() {
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "");
  params.addRequiredCoupledVar("disp_y", "");
  //params.addCoupledVar        ("disp_z", "");

  return params;
}

QQElement::QQElement(const InputParameters & params) :
    Kernel(params),
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
    
    
    _qq_grad_test=new RealVectorValue * [4];
    for (int i=0 ; i < 4; ++i)
        _qq_grad_test[i]=new RealVectorValue [3];
    
    U = new RealTensorValue [4];
    F = new RealTensorValue [4];
    C = new RealTensorValue [4];
    E = new RealTensorValue [4];
    EQQ = new RealTensorValue [6];
    
    
//    _V          = new RealTensorValue **[_mesh.dimension()];
//    _eps_lin    = new RealTensorValue **[_mesh.dimension()];
//    _sigma_lin  = new RealTensorValue **[_mesh.dimension()];
//    
//    for (int d=0; d<_mesh.dimension(); ++d)
//    {
//        _V         [d] = new RealTensorValue *[64];
//        _eps_lin   [d] = new RealTensorValue *[64];
//        _sigma_lin [d] = new RealTensorValue *[64];
//    }
//    
//    for (int d=0; d<_mesh.dimension(); ++d)
//        for (int j = 0; j < 64 ; ++j)
//        {
//            _V         [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
//            _eps_lin   [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
//            _sigma_lin [d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
//        }
    
    
}

void QQElement::computeResidual()
{
    std::cout<<"chiamato\n";
    if (_mesh.dimension()==2)
        computeResidual2D();
}

void QQElement::computeResidual2D()
{
    //std::cout<<"numero di punti "<<_qrule->n_points()<<std::endl;
     for (_qp = 0; _qp < _qrule->n_points(); _qp++)
         std::cout<<_coord[_qp]<<std::endl;
    
    computeGradient(_q_point[0],_q_point[3],_q_point[4],_qq_grad_test[0]);
    computeGradient(_q_point[1],_q_point[5],_q_point[3],_qq_grad_test[1]);
    computeGradient(_q_point[2],_q_point[4],_q_point[5],_qq_grad_test[2]);
    computeGradient(_q_point[0],_q_point[1],_q_point[2],_qq_grad_test[3]);

    RealVectorValue Temp1,Temp2,Temp3;
    
    // Strain triangolino vicino al nodo 0
    
    Temp1=_qq_grad_test[0][0]*_disp_x[0]+_qq_grad_test[0][1]*_disp_x[3]+_qq_grad_test[0][2]*_disp_x[4];
    Temp2=_qq_grad_test[0][0]*_disp_y[0]+_qq_grad_test[0][1]*_disp_y[3]+_qq_grad_test[0][2]*_disp_y[3];
    Temp3=RealVectorValue(0.0,0.0,0.0);
    
    U[0]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));

    // Strain triangolino vicino al nodo 1

    Temp1=_qq_grad_test[1][0]*_disp_x[1]+_qq_grad_test[1][1]*_disp_x[5]+_qq_grad_test[1][2]*_disp_x[3];
    Temp2=_qq_grad_test[1][0]*_disp_y[1]+_qq_grad_test[1][1]*_disp_y[5]+_qq_grad_test[1][2]*_disp_y[3];
    Temp3=RealVectorValue(0.0,0.0,0.0);
    
    U[1]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));

    // Strain triangolino vicino al nodo 2
    
    Temp1=_qq_grad_test[2][0]*_disp_x[2]+_qq_grad_test[2][1]*_disp_x[4]+_qq_grad_test[2][2]*_disp_x[5];
    Temp2=_qq_grad_test[2][0]*_disp_y[2]+_qq_grad_test[2][1]*_disp_y[4]+_qq_grad_test[2][2]*_disp_y[5];
    Temp3=RealVectorValue(0.0,0.0,0.0);
    
    U[2]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));

    // Strain triangolone
    
    Temp1=_qq_grad_test[3][0]*_disp_x[0]+_qq_grad_test[3][1]*_disp_x[1]+_qq_grad_test[3][2]*_disp_x[2];
    Temp2=_qq_grad_test[3][0]*_disp_y[0]+_qq_grad_test[3][1]*_disp_y[1]+_qq_grad_test[3][2]*_disp_y[2];
    Temp3=RealVectorValue(0.0,0.0,0.0);
    
    U[3]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));
    
    for (int i=0; i<4; ++i)
    {
        F[i]=U[i]+_identity;
        C[i]=F[i].transpose()*F[i];
        E[i]=0.5*(C[i]-_identity);
    }
    
    EQQ[0]=2.0*E[0]-E[3];
    EQQ[1]=2.0*E[1]-E[3];
    EQQ[2]=2.0*E[2]-E[3];
    EQQ[3]=0.5*(E[0]+E[1]);
    EQQ[4]=0.5*(E[2]+E[0]);
    EQQ[5]=0.5*(E[1]+E[2]);

    
    // FINE ASSEMBLING STRAIN QQ
    
    for (int i=0; i<6;++i)
    std::cout<<EQQ[i]<<std::endl;
    
    exit(1);

//    std::cout<<_qq_grad_test[2][0]<<std::endl;
//    std::cout<<_qq_grad_test[2][1]<<std::endl;
//    std::cout<<_qq_grad_test[2][2]<<std::endl;
    
    return;
    
    
//    DenseVector<Number> & f_R_x = _assembly.residualBlock(_disp_real_x_var);
//    DenseVector<Number> & f_R_y = _assembly.residualBlock(_disp_real_y_var);
//    DenseVector<Number> & f_I_x = _assembly.residualBlock(_disp_imag_x_var);
//    DenseVector<Number> & f_I_y = _assembly.residualBlock(_disp_imag_y_var);
//    
//    DenseVector<Number> & f_R_p = _assembly.residualBlock(_p_real_var);
//    DenseVector<Number> & f_I_p = _assembly.residualBlock(_p_imag_var);
//    
//    DenseVector<Number> _f_R_local[2];
//    DenseVector<Number> _f_I_local[2];
//    
//    DenseVector<Number> _f_R_p_local;
//    DenseVector<Number> _f_I_p_local;
//    
//    for (int i=0; i<2; ++i)
//    {
//        _f_R_local[i].resize(f_R_x.size());
//        _f_R_local[i].zero();
//        _f_I_local[i].resize(f_R_x.size());
//        _f_I_local[i].zero();
//    }
//    _f_R_p_local.resize(f_R_p.size());
//    _f_R_p_local.zero();
//    _f_I_p_local.resize(f_I_p.size());
//    _f_I_p_local.zero();
//    
//    
//    
//    for (int d=0; d<2; ++d)
//        for (_i = 0; _i < _P2.size(); _i++)
//            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//            {
//                RealTensorValue P_R=_sigma_real[_qp]-_alpha[_qp]*_pres_real[_qp]*_identity;
//                RealTensorValue P_I=_sigma_imag[_qp]-_alpha[_qp]*_pres_imag[_qp]*_identity;
//                
//                _f_R_local[d](_i) += _JxW[_qp] * _coord[_qp] *
//                (
//                 P_R(d , 0) * _grad_P2[_i][_qp] (0)+
//                 P_R(d , 1) * _grad_P2[_i][_qp] (1)
//                 );
//                
//                _f_I_local[d](_i) += _JxW[_qp] * _coord[_qp] *
//                (
//                 P_I(d , 0) * _grad_P2[_i][_qp] (0)+
//                 P_I(d , 1) * _grad_P2[_i][_qp] (1)
//                 );
//
//            }
//    
//    for (_i=0; _i<_P1.size(); ++_i)
//        for (_qp=0; _qp<_qrule->n_points(); ++_qp)
//        {
//            _f_R_p_local(_i)+=_JxW[_qp] * _coord[_qp] *(_alpha[_qp]*_tr_eps_imag[_qp]+_inv_m[_qp]*_pres_imag[_qp])*_P1[_i][_qp];
//            _f_R_p_local(_i)-=_JxW[_qp] * _coord[_qp] *_diffusion[_qp]*_grad_pres_real[_qp]*_grad_P1[_i][_qp];
//            
//            _f_I_p_local(_i)+=-1.0*_JxW[_qp] * _coord[_qp] *(_alpha[_qp]*_tr_eps_real[_qp]+_inv_m[_qp]*_pres_real[_qp])*_P1[_i][_qp];
//            _f_I_p_local(_i)-=_JxW[_qp] * _coord[_qp] *_diffusion[_qp]*_grad_pres_imag[_qp]*_grad_P1[_i][_qp];
//            
//        }
//    
//    f_R_x += _f_R_local[0];
//    f_R_y += _f_R_local[1];
//    f_I_x += _f_I_local[0];
//    f_I_y += _f_I_local[1];
//    
//    f_R_p += _f_R_p_local;
//    f_I_p += _f_I_p_local;
}

void QQElement::computeOffDiagJacobian(unsigned int jvar)
{
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
}

void QQElement::computeJacobian2D()
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
    
}

void QQElement::initTensorVariables()
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
}

void QQElement::computeGradient(RealVectorValue x0, RealVectorValue x1, RealVectorValue x2, RealVectorValue * Gradient)
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

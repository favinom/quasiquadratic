#include "GuccioneStrain_slow.h"

#include "MooseMesh.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

template <>
InputParameters validParams<GuccioneStrain_slow>() {
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "");
  params.addRequiredCoupledVar("disp_y", "");
  //params.addCoupledVar        ("disp_z", "");

 params.addRequiredParam<unsigned>("component", "component");
    
  return params;
}

GuccioneStrain_slow::GuccioneStrain_slow(const InputParameters & params) :
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
    
    _mu=0.1;
    _lambda=100.0;
    strain_projection=1;

}


void GuccioneStrain_slow::computeResidual()
{
    if (_mesh.dimension()==2)
        computeResidual2D();
}

void GuccioneStrain_slow::computeResidual2D()
{
    Real area=0.0;
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        area+=_JxW[_qp];

//    DenseMatrix<Number> mass_matrix;
//    mass_matrix.resize(3,3);
//    for (int i=0; i<3; ++i)
//        for (int j=0; j<3; ++j)
//            if (i==j)
//                mass_matrix(i,j)=area/6.0;
//            else
//                mass_matrix(i,j)=area/12.0;
    
    
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

    }
    
    for (int nodo_coarse=0; nodo_coarse<3; ++nodo_coarse)
    {
        CQQ[nodo_coarse]=2.0*C[nodo_coarse]-C[3];
    }

    RealTensorValue SQQ_qp[3];
    if (strain_projection)
    {

        RealTensorValue CQQ_qp[3];
        CQQ_qp[0]=2.0/3.0*CQQ[0]+1.0/6.0*CQQ[1]+1.0/6.0*CQQ[2];
        CQQ_qp[1]=1.0/6.0*CQQ[0]+2.0/3.0*CQQ[1]+1.0/6.0*CQQ[2];
        CQQ_qp[2]=1.0/6.0*CQQ[0]+1.0/6.0*CQQ[1]+2.0/3.0*CQQ[2];
        
        for (_qp=0; _qp<3; ++_qp)
        {
            SQQ_qp[_qp]=assembleStress(CQQ_qp[_qp]);
        }
    }
    
    
    
    // FINE ASSEMBLING STRAIN QQ

    for (int dim=0; dim<2; ++dim)
        for (int tri=0; tri<4; ++tri)
            assembleStrainLin(dim,tri,_local_to_global[tri]);
    
    for (int dim=0; dim<2; ++dim)
        for (int nodo=0; nodo<6; ++nodo)
            for (int nodo_coarse=0; nodo_coarse<3; ++nodo_coarse)
                _eps_lin_QQ[dim][nodo][nodo_coarse]=2.0*_eps_lin[dim][nodo_coarse][nodo]-_eps_lin[dim][3][nodo];
    
    RealTensorValue _eps_lin_QQ_qp[2][6][3];
    if (strain_projection)
    {
        for (int dim=0; dim<2; ++dim)
        {
            for (int nodo=0; nodo<6; ++nodo)
            {
                _eps_lin_QQ_qp[dim][nodo][0]=2.0/3.0*_eps_lin_QQ[dim][nodo][0]+1.0/6.0*_eps_lin_QQ[dim][nodo][1]+1.0/6.0*_eps_lin_QQ[dim][nodo][2];
                _eps_lin_QQ_qp[dim][nodo][1]=1.0/6.0*_eps_lin_QQ[dim][nodo][0]+2.0/3.0*_eps_lin_QQ[dim][nodo][1]+1.0/6.0*_eps_lin_QQ[dim][nodo][2];
                _eps_lin_QQ_qp[dim][nodo][2]=1.0/6.0*_eps_lin_QQ[dim][nodo][0]+1.0/6.0*_eps_lin_QQ[dim][nodo][1]+2.0/3.0*_eps_lin_QQ[dim][nodo][2];
            }
        }
        
    }
    
    DenseVector<Number> & f_x = _assembly.residualBlock(_disp_x_var);
    DenseVector<Number> & f_y = _assembly.residualBlock(_disp_y_var);
    DenseVector<Number> _f_local[2];
    
    for (int dim=0; dim<2; ++dim)
    {
        _f_local[dim].resize(f_x.size());
        _f_local[dim].zero();
    }
    

        for (int dim=0; dim<2; ++dim)
            for (_j=0; _j<6; ++_j)
                for (_qp =0; _qp<3; ++_qp)
                    _f_local[dim](simpson_to_tri6[_j])+=1.0/3.0*area*SQQ_qp[_qp].contract(_eps_lin_QQ_qp[dim][_j][_qp]);
        
    
    
    
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


void GuccioneStrain_slow::computeGradient(RealVectorValue x0, RealVectorValue x1, RealVectorValue x2, RealVectorValue * Gradient)
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

void GuccioneStrain_slow::assembleStrainLin(int dim, int triangle, int * _local_to_global)
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


RealTensorValue GuccioneStrain_slow::assembleStress(RealTensorValue const & C)
{
    RealTensorValue E = 0.5 * (C - _identity);
    
    Real _stiffening =  std::exp(_mu*E.contract(E)+_lambda/2.0 *E.tr()*E.tr());
    
    RealTensorValue S = _stiffening * ( 2.0 * _mu* E +_lambda*E.tr()*_identity );
    
    return S;

}





















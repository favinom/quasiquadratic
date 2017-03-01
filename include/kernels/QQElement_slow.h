#ifndef QQElement_slow_H
#define QQElement_slow_H

#include "Kernel.h"

class QQElement_slow;

template<>
InputParameters validParams<QQElement_slow>();


class QQElement_slow : public Kernel
{
public:
    QQElement_slow( InputParameters const & params);
protected:
    virtual Real computeQpResidual(){return 0.0;};
    //virtual Real computeQpJacobian(){return 0.0;};
    //virtual Real computeQpOffDiagJacobian(unsigned int jvar){std::cout<<jvar<<std::endl; return 0.0;};

    virtual void computeResidual();
    //virtual void computeJacobian(){};
    //virtual void computeOffDiagJacobian(unsigned int jvar);
    //virtual void computeOffDiagJacobianScalar(unsigned int jvar){std::cout<<jvar<<std::endl;};
    
    //void initTensorVariables();
    void computeResidual2D();
    //void computeJacobian2D();
    
    virtual void assembleStrainLin(int dim, int triangle, int * _local_to_global);
    
    void computeGradient(RealVectorValue x0, RealVectorValue x1, RealVectorValue x2, RealVectorValue * Gradient);

    int _component;
    
    unsigned int _disp_x_var;
    unsigned int _disp_y_var;
//    unsigned int _disp_z_var;

    VariableValue    const & _disp_x;
    VariableValue    const & _disp_y;
//    VariableValue    const & _disp_z;
    VariableGradient const & _grad_disp_x;
    VariableGradient const & _grad_disp_y;
//    VariableGradient const & _grad_disp_z;

    RealTensorValue _identity;
    RealTensorValue _V;
    RealTensorValue ***_eps_lin;
    
    RealTensorValue ***_eps_lin_QQ;
    
    RealTensorValue * EQQ_qp;
    RealTensorValue *** Elin_qp;
    
    
//    RealTensorValue ***_sigma_lin;
    
    RealVectorValue **_qq_grad_test;
    RealTensorValue *U;
    RealTensorValue *F;
    RealTensorValue *C;
    RealTensorValue *E;
    RealTensorValue *EQQ;

    int **_local_to_global;
    
    int * simpson_to_tri6;
    
    Real _mu,_lambda;
};

#endif 

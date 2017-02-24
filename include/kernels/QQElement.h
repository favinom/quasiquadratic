#ifndef QQElement_H
#define QQElement_H

#include "Kernel.h"

class QQElement;

template<>
InputParameters validParams<QQElement>();


class QQElement : public Kernel
{
public:
    QQElement( InputParameters const & params);
protected:
    virtual Real computeQpResidual(){return 0.0;};
    virtual Real computeQpJacobian(){return 0.0;};
    virtual Real computeQpOffDiagJacobian(unsigned int jvar){std::cout<<jvar<<std::endl; return 0.0;};

    virtual void computeResidual();
    virtual void computeJacobian(){};
    virtual void computeOffDiagJacobian(unsigned int jvar);
    virtual void computeOffDiagJacobianScalar(unsigned int jvar){std::cout<<jvar<<std::endl;};
    
    void initTensorVariables();
    void computeResidual2D();
    void computeJacobian2D();
    
    void computeGradient(RealVectorValue x0, RealVectorValue x1, RealVectorValue x2, RealVectorValue * Gradient);

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
    RealTensorValue ***_V;
    RealTensorValue ***_eps_lin;
    RealTensorValue ***_sigma_lin;
    
    RealVectorValue **_qq_grad_test;
    RealTensorValue *U;
    RealTensorValue *F;
    RealTensorValue *C;
    RealTensorValue *E;
    RealTensorValue *EQQ;
    
    Real _mu,_lambda;
};

#endif 

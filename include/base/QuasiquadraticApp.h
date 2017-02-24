#ifndef QUASIQUADRATICAPP_H
#define QUASIQUADRATICAPP_H

#include "MooseApp.h"

class QuasiquadraticApp;

template<>
InputParameters validParams<QuasiquadraticApp>();

class QuasiquadraticApp : public MooseApp
{
public:
  QuasiquadraticApp(InputParameters parameters);
  virtual ~QuasiquadraticApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* QUASIQUADRATICAPP_H */

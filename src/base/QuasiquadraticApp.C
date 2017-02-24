#include "QuasiquadraticApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template<>
InputParameters validParams<QuasiquadraticApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

QuasiquadraticApp::QuasiquadraticApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  QuasiquadraticApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  QuasiquadraticApp::associateSyntax(_syntax, _action_factory);
}

QuasiquadraticApp::~QuasiquadraticApp()
{
}

// External entry point for dynamic application loading
extern "C" void QuasiquadraticApp__registerApps() { QuasiquadraticApp::registerApps(); }
void
QuasiquadraticApp::registerApps()
{
  registerApp(QuasiquadraticApp);
}

// External entry point for dynamic object registration
extern "C" void QuasiquadraticApp__registerObjects(Factory & factory) { QuasiquadraticApp::registerObjects(factory); }
void
QuasiquadraticApp::registerObjects(Factory & factory)
{
}

// External entry point for dynamic syntax association
extern "C" void QuasiquadraticApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { QuasiquadraticApp::associateSyntax(syntax, action_factory); }
void
QuasiquadraticApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

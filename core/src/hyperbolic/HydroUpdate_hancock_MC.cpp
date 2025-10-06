#include "states/State_hydro_MC.h"

#include "hyperbolic/policy/HyperbolicPolicy_Hydro_MC.h"

#include "hyperbolic/scheme/Hyperbolic_hancock_MC.h"

namespace dyablo{

class HydroUpdate_hancock_MC
  : public Hyperbolic_hancock_MC<HyperbolicPolicy_Hydro_MC>
{
public:
  using Hyperbolic_hancock_MC<HyperbolicPolicy_Hydro_MC>::Hyperbolic_hancock_MC;
};

} //namespace dyablo

FACTORY_REGISTER( dyablo::HyperbolicUpdateFactory, 
                  dyablo::HydroUpdate_hancock_MC, 
                  "HydroUpdate_hancock_MC")
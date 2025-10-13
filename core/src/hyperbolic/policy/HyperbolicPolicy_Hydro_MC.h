#pragma once

#include "states/State_Ops.h"
#include "HyperbolicPolicy_base.h"
#include "HyperbolicPolicy_Slope.h"
#include "HyperbolicPolicy_BoundaryConditions.h"

namespace dyablo{

/**
 * @brief Structure holding conservative hydrodynamics variables
 **/ 
struct HyperbolicPolicy_ConsHydroMCState {
  enum VarIndex : dyablo::VarIndex
  {
    Irho,
    Ie_tot,
    Irho_vx,
    Irho_vy,
    Irho_vz,
    Iflux_x_l,
    Iflux_x_r,
    Iflux_y_l,
    Iflux_y_r,
    Iflux_z_l,
    Iflux_z_r,
    Irho_old
  }; 

  static FieldManager getFieldManager()
  {
    return FieldManager( {VarIndex::Irho, VarIndex::Ie_tot, VarIndex::Irho_vx, VarIndex::Irho_vy, VarIndex::Irho_vz, 
                          VarIndex::Iflux_x_l, VarIndex::Iflux_x_r, VarIndex::Iflux_y_l, VarIndex::Iflux_y_r, VarIndex::Iflux_z_l, VarIndex::Iflux_z_r, VarIndex::Irho_old } );
  } 

  real_t rho = 0;
  real_t e_tot = 0;
  real_t rho_u = 0;
  real_t rho_v = 0;
  real_t rho_w = 0;
  
  real_t flux_x_l = 0;
  real_t flux_x_r = 0;
  real_t flux_y_l = 0;
  real_t flux_y_r = 0;
  real_t flux_z_l = 0;
  real_t flux_z_r = 0;
  real_t rho_old  = 0;
};

DECLARE_STATE_TYPE( HyperbolicPolicy_ConsHydroMCState, 12 );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 0, rho );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 1, e_tot );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 2, rho_u );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 3, rho_v );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 4, rho_w );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 5, flux_x_l );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 6, flux_x_r );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 7, flux_y_l );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 8, flux_y_r );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 9, flux_z_l );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 10, flux_z_r );
DECLARE_STATE_GET( HyperbolicPolicy_ConsHydroMCState, 11, rho_old );

/**
 * @brief Structure holding primitive hydrodynamics variables
 */
struct HyperbolicPolicy_PrimHydroMCState {
  enum VarIndex : dyablo::VarIndex
  {
    Irho,
    Ip,
    Iu,
    Iv,
    Iw,
    Iflux_x_l,
    Iflux_x_r,
    Iflux_y_l,
    Iflux_y_r,
    Iflux_z_l,
    Iflux_z_r,
    Irho_old
  };

  static FieldManager getFieldManager()
  {
    return FieldManager( {VarIndex::Irho, VarIndex::Ip, VarIndex::Iu, VarIndex::Iv, VarIndex::Iw, 
                          VarIndex::Iflux_x_l, VarIndex::Iflux_x_r, VarIndex::Iflux_y_l, VarIndex::Iflux_y_r, VarIndex::Iflux_z_l, VarIndex::Iflux_z_r, VarIndex::Irho_old} );
  }

  real_t rho = 0;
  real_t p = 0;
  real_t u = 0;
  real_t v = 0;
  real_t w = 0;
  real_t flux_x_l = 0;
  real_t flux_x_r = 0;
  real_t flux_y_l = 0;
  real_t flux_y_r = 0;
  real_t flux_z_l = 0;
  real_t flux_z_r = 0;
  real_t rho_old  = 0;
};

DECLARE_STATE_TYPE( HyperbolicPolicy_PrimHydroMCState, 12 );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 0, rho );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 1, p );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 2, u );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 3, v );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 4, w );

DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 5, flux_x_l );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 6, flux_x_r );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 7, flux_y_l );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 8, flux_y_r );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 9, flux_z_l );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 10, flux_z_r );
DECLARE_STATE_GET( HyperbolicPolicy_PrimHydroMCState, 11, rho_old );


struct HyperbolicPolicy_Hydro_MC_Params
{
  static HyperbolicPolicy_Hydro_MC_Params from_configMap( ConfigMap& configMap )
  {
    return {
      .ndim   = configMap.getValue<int>("mesh", "ndim", 3),
      .gamma0 = configMap.getValue<real_t>("hydro", "gamma0", 1.4),
      .smallr = configMap.getValue<real_t>("hydro", "smallr", 1e-10),
      .smallp = configMap.getValue<real_t>("hydro", "smallp", 1e-10),
      .smallc = configMap.getValue<real_t>("hydro", "smallc", 1e-10),
    };
  }

  int ndim;
  real_t gamma0;
  real_t smallr;
  real_t smallp;
  real_t smallc;
};

class HyperbolicPolicy_State_Hydro_MC
{
private:
  int ndim;
  real_t gamma0;

  using CellIndex = ForeachCell::CellIndex;
  using FieldAccessor = UserData::FieldAccessor;
public:
  using PrimState = HyperbolicPolicy_PrimHydroMCState;
  using ConsState = HyperbolicPolicy_ConsHydroMCState;

  HyperbolicPolicy_State_Hydro_MC( const HyperbolicPolicy_Hydro_MC_Params& params )
  : ndim(params.ndim),
    gamma0(params.gamma0)
  {}

  using ConsVarIndex = HyperbolicPolicy_ConsHydroMCState::VarIndex;

  FieldAccessor getUin( UserData& U ) const
  {
    std::vector<FieldAccessor::FieldInfo> Uin_fieldinfo { 
      {"rho",      ConsVarIndex::Irho}, 
      {"e_tot",    ConsVarIndex::Ie_tot},
      {"rho_vx",   ConsVarIndex::Irho_vx},
      {"rho_vy",   ConsVarIndex::Irho_vy},
      {"rho_vz",   ConsVarIndex::Irho_vz},
      {"flux_x_l", ConsVarIndex::Iflux_x_l},
      {"flux_x_r", ConsVarIndex::Iflux_x_r},
      {"flux_y_l", ConsVarIndex::Iflux_y_l},
      {"flux_y_r", ConsVarIndex::Iflux_y_r},
      {"flux_z_l", ConsVarIndex::Iflux_z_l},
      {"flux_z_r", ConsVarIndex::Iflux_z_r},
      {"rho_old",  ConsVarIndex::Irho_old},
    };

    return U.getAccessor( Uin_fieldinfo );
  }

  FieldAccessor getUout( UserData& U ) const
  {
    std::vector<FieldAccessor::FieldInfo> Uout_fieldinfo { 
      {"rho_next",     ConsVarIndex::Irho}, 
      {"e_tot_next",   ConsVarIndex::Ie_tot},
      {"rho_vx_next",  ConsVarIndex::Irho_vx},
      {"rho_vy_next",  ConsVarIndex::Irho_vy},
      {"rho_vz_next",  ConsVarIndex::Irho_vz},
      {"flux_x_l",     ConsVarIndex::Iflux_x_l},
      {"flux_x_r",     ConsVarIndex::Iflux_x_r},
      {"flux_y_l",     ConsVarIndex::Iflux_y_l},
      {"flux_y_r",     ConsVarIndex::Iflux_y_r},
      {"flux_z_l",     ConsVarIndex::Iflux_z_l},
      {"flux_z_r",     ConsVarIndex::Iflux_z_r},
      {"rho_old",      ConsVarIndex::Irho_old},
    };
    return U.getAccessor( Uout_fieldinfo );
  }

  template < typename Array_t >
  KOKKOS_INLINE_FUNCTION
  ConsState getConsState( const Array_t& U, const CellIndex& iCell ) const
  {
    ConsState u;
    u.rho   = U.at(iCell, ConsState::VarIndex::Irho );
    u.rho_old  = U.at(iCell, ConsState::VarIndex::Irho_old );
    u.e_tot = U.at(iCell, ConsState::VarIndex::Ie_tot );
    u.rho_u = U.at(iCell, ConsState::VarIndex::Irho_vx );
    u.rho_v = U.at(iCell, ConsState::VarIndex::Irho_vy );
    u.rho_w = (ndim == 3 ? U.at(iCell, ConsState::VarIndex::Irho_vz ) : 0.0);

    u.flux_x_l = U.at(iCell, ConsState::VarIndex::Iflux_x_l );
    u.flux_x_r = U.at(iCell, ConsState::VarIndex::Iflux_x_r );
    u.flux_y_l = U.at(iCell, ConsState::VarIndex::Iflux_y_l );
    u.flux_y_r = U.at(iCell, ConsState::VarIndex::Iflux_y_r );
    u.flux_z_l = (ndim == 3 ? U.at(iCell, ConsState::VarIndex::Iflux_z_l ) : 0.0);
    u.flux_z_r = (ndim == 3 ? U.at(iCell, ConsState::VarIndex::Iflux_z_r ) : 0.0);
    return u;
  }

  template < typename Array_t >
  KOKKOS_INLINE_FUNCTION
  void setConsState( const Array_t& U, const CellIndex& iCell, const ConsState& u ) const
  {
    U.at(iCell, ConsState::VarIndex::Irho) = u.rho;
    U.at(iCell, ConsState::VarIndex::Irho_old) = u.rho_old;
    U.at(iCell, ConsState::VarIndex::Ie_tot) = u.e_tot;
    U.at(iCell, ConsState::VarIndex::Irho_vx) = u.rho_u;
    U.at(iCell, ConsState::VarIndex::Irho_vy) = u.rho_v;

    U.at(iCell, ConsState::VarIndex::Iflux_x_l) = u.flux_x_l;
    U.at(iCell, ConsState::VarIndex::Iflux_x_r) = u.flux_x_r;
    U.at(iCell, ConsState::VarIndex::Iflux_y_l) = u.flux_y_l;
    U.at(iCell, ConsState::VarIndex::Iflux_y_r) = u.flux_y_r;

    if (ndim == 3){
      U.at(iCell, ConsState::VarIndex::Irho_vz) = u.rho_w;
      U.at(iCell, ConsState::VarIndex::Iflux_z_l) = u.flux_z_l;
      U.at(iCell, ConsState::VarIndex::Iflux_z_r) = u.flux_z_r;
    }
  }

  template < typename Array_t >
  KOKKOS_INLINE_FUNCTION
  void atomic_addConsState( const Array_t& U, const CellIndex& iCell, const ConsState& u ) const
  {
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Irho), u.rho);
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Irho_old), u.rho_old);
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Ie_tot),  u.e_tot);
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Irho_vx), u.rho_u);
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Irho_vy), u.rho_v);

    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Iflux_x_l), u.flux_x_l);
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Iflux_x_r), u.flux_x_r);
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Iflux_y_l), u.flux_y_l);
    Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Iflux_y_r), u.flux_y_r);
    if (ndim == 3){
      Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Irho_vz), u.rho_w);
      Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Iflux_z_l), u.flux_z_l);
      Kokkos::atomic_add(&U.at(iCell, ConsState::VarIndex::Iflux_z_r), u.flux_z_r);
    }
  }

  template < typename Array_t >
  KOKKOS_INLINE_FUNCTION
  PrimState getPrimState( const Array_t& Q, const CellIndex& iCell ) const
  {
    PrimState q;
    q.rho = Q.at(iCell, PrimState::VarIndex::Irho );
    q.rho_old = Q.at(iCell, PrimState::VarIndex::Irho_old );
    q.p   = Q.at(iCell, PrimState::VarIndex::Ip );
    q.u   = Q.at(iCell, PrimState::VarIndex::Iu );
    q.v   = Q.at(iCell, PrimState::VarIndex::Iv );
    q.flux_x_l = Q.at(iCell, PrimState::VarIndex::Iflux_x_l);
    q.flux_x_r = Q.at(iCell, PrimState::VarIndex::Iflux_x_r);
    q.flux_y_l = Q.at(iCell, PrimState::VarIndex::Iflux_y_l);
    q.flux_y_r = Q.at(iCell, PrimState::VarIndex::Iflux_y_r);

    q.w   = (ndim == 3 ? Q.at(iCell, PrimState::VarIndex::Iw ) : 0.0);
    q.flux_z_l   = (ndim == 3 ? Q.at(iCell, PrimState::VarIndex::Iflux_z_l ) : 0.0);
    q.flux_z_r   = (ndim == 3 ? Q.at(iCell, PrimState::VarIndex::Iflux_z_r ) : 0.0);
    return q;
  }

  template < typename Array_t >
  KOKKOS_INLINE_FUNCTION
  void setPrimState( const Array_t& Q, const CellIndex& iCell, const PrimState& q ) const
  {
    Q.at(iCell, PrimState::VarIndex::Irho) = q.rho;
    Q.at(iCell, PrimState::VarIndex::Irho_old) = q.rho_old;
    Q.at(iCell, PrimState::VarIndex::Ip) = q.p;
    Q.at(iCell, PrimState::VarIndex::Iu) = q.u;
    Q.at(iCell, PrimState::VarIndex::Iv) = q.v;

    Q.at(iCell, PrimState::VarIndex::Iflux_x_l) = q.flux_x_l;
    Q.at(iCell, PrimState::VarIndex::Iflux_x_r) = q.flux_x_r;
    Q.at(iCell, PrimState::VarIndex::Iflux_y_l) = q.flux_y_l;
    Q.at(iCell, PrimState::VarIndex::Iflux_y_r) = q.flux_y_r;

    if (ndim == 3){
      Q.at(iCell, PrimState::VarIndex::Iw) = q.w;
      Q.at(iCell, PrimState::VarIndex::Iflux_z_l) = q.flux_z_l;
      Q.at(iCell, PrimState::VarIndex::Iflux_z_r) = q.flux_z_r;
    }
  }

  KOKKOS_INLINE_FUNCTION
  PrimState consToPrim( const ConsState& U ) const
  {
    real_t gamma0 = this->gamma0;

    const real_t Ek = 0.5 * (U.rho_u*U.rho_u+U.rho_v*U.rho_v+U.rho_w*U.rho_w)/U.rho;
    const real_t p = (U.e_tot - Ek) * (gamma0-1.0);
    return {U.rho, 
            p, 
            U.rho_u/U.rho, 
            U.rho_v/U.rho, 
            (ndim == 3 ? U.rho_w/U.rho : 0.0),
            U.flux_x_l,
            U.flux_x_r,
            U.flux_y_l,
            U.flux_y_r,
            (ndim == 3 ? U.flux_z_l : 0.0),
            (ndim == 3 ? U.flux_z_r : 0.0),
            U.rho_old};
  }

  KOKKOS_INLINE_FUNCTION
  ConsState primToCons( const PrimState& Q ) const
  {
    real_t gamma0 = this->gamma0;

    const real_t Ek = 0.5 * Q.rho * (Q.u*Q.u+Q.v*Q.v+Q.w*Q.w);
    const real_t E  = Ek + Q.p / (gamma0-1.0);
    return {Q.rho, 
            E, 
            Q.rho*Q.u, 
            Q.rho*Q.v, 
            (ndim ==3 ? Q.rho*Q.w : 0.0),
            Q.flux_x_l,
            Q.flux_x_r,
            Q.flux_y_l,
            Q.flux_y_r,
            (ndim == 3 ? Q.flux_z_l : 0.0),
            (ndim == 3 ? Q.flux_z_r : 0.0),
            Q.rho_old};
  }
};

class HyperbolicPolicy_RiemannSolver_Hydro_hllc_MC
{
private:
  struct Rparams {
    real_t gamma0;
    real_t smallr;
    real_t smallp;
    real_t smallc;
  } rparams;

public:
  using State = HyperbolicPolicy_State_Hydro_MC;
  using PrimState = State::PrimState;
  using ConsState = State::ConsState;

  HyperbolicPolicy_RiemannSolver_Hydro_hllc_MC( const HyperbolicPolicy_Hydro_MC_Params& params )
  : rparams( 
    {
      .gamma0 = params.gamma0,
      .smallr = params.smallr,
      .smallp = params.smallp,
      .smallc = params.smallc,
    })
  {}

  KOKKOS_INLINE_FUNCTION
  ConsState riemann_solver( PrimState qL, PrimState qR, ComponentIndex3D dir ) const
  {
    qL = swapComponents(qL, dir);
    qR = swapComponents(qR, dir);
    ConsState flux = riemann_hllc(qL, qR);
    flux = swapComponents(flux, dir);
    return flux;
  }

  KOKKOS_INLINE_FUNCTION
  ConsState riemann_solver( PrimState qL, PrimState qR, ComponentIndex3D dir, FlowDirection fdir ) const
  {
    qL = swapComponents(qL, dir);
    qR = swapComponents(qR, dir);
    ConsState flux = riemann_hllc(qL, qR, fdir);
    flux = swapComponents(flux, dir);
    return flux;
  }

private:

  KOKKOS_INLINE_FUNCTION
  PrimState swapComponents(const PrimState &q, ComponentIndex3D comp) const
  {
    switch( comp )
    {
      case IX:
        return q;
      case IY:
        return PrimState{q.rho, q.p, q.v, q.u, q.w, q.flux_y_l, q.flux_y_r, q.flux_x_l, q.flux_x_r, q.flux_z_l, q.flux_z_r, q.rho_old};
      case IZ:
        return PrimState{q.rho, q.p, q.w, q.v, q.u, q.flux_z_l, q.flux_z_r, q.flux_y_l, q.flux_y_r, q.flux_x_l, q.flux_x_r, q.rho_old};
      default:
        DYABLO_ASSERT_KOKKOS_DEBUG(false, "invalid component");
        return PrimState{};
    }
  }

  KOKKOS_INLINE_FUNCTION
  ConsState swapComponents(const ConsState &u, ComponentIndex3D comp) const
  {
    switch( comp )
    {
      case IX:
        return u;
      case IY:
        return ConsState{u.rho, u.e_tot, u.rho_v, u.rho_u, u.rho_w, u.flux_y_l, u.flux_y_r, u.flux_x_l, u.flux_x_r, u.flux_z_l, u.flux_z_r, u.rho_old};
      case IZ:
        return ConsState{u.rho, u.e_tot, u.rho_w, u.rho_v, u.rho_u, u.flux_z_l, u.flux_z_r, u.flux_y_l, u.flux_y_r, u.flux_x_l, u.flux_x_r, u.rho_old};
      default:
        DYABLO_ASSERT_KOKKOS_DEBUG(false, "invalid component");
        return ConsState{};
    }
  }

  KOKKOS_INLINE_FUNCTION
  ConsState riemann_hllc( PrimState qleft, PrimState qright) const
  {
    real_t gamma0 = rparams.gamma0;
    real_t smallr = rparams.smallr;
    real_t smallp = rparams.smallp;
    real_t smallc = rparams.smallc;

    const real_t entho = 1 / (gamma0 - 1);

    // Left variables
    real_t rl = fmax(qleft.rho, smallr);
    real_t pl = fmax(qleft.p, rl*smallp);
    real_t ul =      qleft.u;
    real_t vl =      qleft.v;
    real_t wl =      qleft.w;

    real_t ecinl = 0.5*rl*(ul*ul+vl*vl+wl*wl);
    real_t etotl = pl*entho+ecinl;
    real_t ptotl = pl;


    // Right variables
    real_t rr = fmax(qright.rho, smallr);
    real_t pr = fmax(qright.p, rr*smallp);
    real_t ur =      qright.u;
    real_t vr =      qright.v;
    real_t wr =      qright.w;

    real_t ecinr = 0.5*rr*(ur*ur+vr*vr+wr*wr);
    real_t etotr = pr*entho+ecinr;
    real_t ptotr = pr;
    
    // Find the largest eigenvalues in the normal direction to the interface
    real_t cfastl = SQRT(fmax(gamma0*pl/rl,smallc*smallc));
    real_t cfastr = SQRT(fmax(gamma0*pr/rr,smallc*smallc));

    // Compute HLL wave speed
    real_t SL = fmin(ul,ur) - fmax(cfastl,cfastr);
    real_t SR = fmax(ul,ur) + fmax(cfastl,cfastr);

    // Compute lagrangian sound speed
    real_t rcl = rl*(ul-SL);
    real_t rcr = rr*(SR-ur);
    
    // Compute acoustic star state
    real_t ustar    = (rcr*ur   +rcl*ul   +  (ptotl-ptotr))/(rcr+rcl);
    real_t ptotstar = (rcr*ptotl+rcl*ptotr+rcl*rcr*(ul-ur))/(rcr+rcl);

    // Left star region variables
    real_t rstarl    = rl*(SL-ul)/(SL-ustar);
    real_t etotstarl = ((SL-ul)*etotl-ptotl*ul+ptotstar*ustar)/(SL-ustar);
    
    // Right star region variables
    real_t rstarr    = rr*(SR-ur)/(SR-ustar);
    real_t etotstarr = ((SR-ur)*etotr-ptotr*ur+ptotstar*ustar)/(SR-ustar);
    
    // Sample the solution at x/t=0
    real_t ro, uo, ptoto, etoto;
    //real_t p_out;
    if (SL > 0) {
      ro=rl;
      uo=ul;
      ptoto=ptotl;
      etoto=etotl;
      //p_out=pl;
    } else if (ustar > 0) {
      ro=rstarl;
      uo=ustar;
      ptoto=ptotstar;
      etoto=etotstarl;
      //p_out=ptotstar;
    } else if (SR > 0) {
      ro=rstarr;
      uo=ustar;
      ptoto=ptotstar;
      etoto=etotstarr;
      //p_out=ptotstar;
    } else {
      ro=rr;
      uo=ur;
      ptoto=ptotr;
      etoto=etotr;
      //p_out=pr;
    }
      
    // Compute the Godunov flux
    ConsState flux;
    flux.rho   = ro*uo;
    flux.rho_u = ro*uo*uo+ptoto;
    flux.e_tot = (etoto+ptoto)*uo;
    if (flux.rho > 0) {
      flux.rho_v = flux.rho*qleft.v;
      flux.rho_w = flux.rho*qleft.w;
    } else {
      flux.rho_v = flux.rho*qright.v;
      flux.rho_w = flux.rho*qright.w;
    }
    return flux;
  }

  KOKKOS_INLINE_FUNCTION
  ConsState riemann_hllc( PrimState qleft, PrimState qright, FlowDirection fdir) const
  {
    ConsState flux = riemann_hllc(qleft, qright);

    

    switch (fdir) {
      case FlowDirection::LEFT:  flux.flux_x_l = flux.rho; break;
      case FlowDirection::RIGHT:  flux.flux_x_r = flux.rho; break;
      default:                   throw std::runtime_error(std::string("The flow doesn't have three directions."));
    }

    // switch (fdir) {
    //   case FlowDirection::LEFT:
    //     std::cout << "je suis ici (left)  avec rho " << flux.rho << " et flux_x_l=" << flux.flux_x_l << " et flux_x_r" << flux.flux_x_r << std::endl;
    //     break;
    //   case FlowDirection::RIGHT:
    //     std::cout << "je suis ici (right) avec rho " << flux.rho << " et flux_x_l=" << flux.flux_x_l << " et flux_x_r" << flux.flux_x_r << std::endl;
    //     break;
    //   default:
    //     std::cerr << "something bizarre happen" << std::endl;
    //     break;
    // }
    // if ((flux.flux_x_l) || (flux.flux_x_r))
    //   switch (fdir) {
    //     case FlowDirection::LEFT:
    //       std::cout << "je suis ici (left)  avec rho " << flux.rho << " et flux_x_l=" << flux.flux_x_l << " et flux_x_r" << flux.flux_x_r << std::endl;
    //       break;
    //     case FlowDirection::RIGHT:
    //       std::cout << "je suis ici (right) avec rho " << flux.rho << " et flux_x_l=" << flux.flux_x_l << " et flux_x_r" << flux.flux_x_r << std::endl;
    //       break;
    //     default:
    //       std::cerr << "something bizarre happen" << std::endl;
    //       break;

    //   }
      


    return flux;
  }

};

} // namespace dyablo

#include "HyperbolicPolicy_Hydro_BoundaryConditions.h"

namespace dyablo {

using HyperbolicPolicy_BoundaryConditions_Hydro_MC = HyperbolicPolicy_BoundaryConditions_Hydro<HyperbolicPolicy_State_Hydro_MC>;
using HyperbolicPolicy_BoundaryConditions_Hydro_DoubleMach_MC = HyperbolicPolicy_BoundaryConditions_Hydro_DoubleMach_Template<HyperbolicPolicy_State_Hydro_MC>;

using HyperbolicPolicy_BoundaryConditions_Hydro_dynamic = HyperbolicPolicy_BoundaryConditions_dynamic<
    HyperbolicPolicy_State_Hydro_MC,
    HyperbolicPolicy_BoundaryConditions_Hydro_MC,
    HyperbolicPolicy_BoundaryConditions_Hydro_DoubleMach_MC
  >;

class HyperbolicPolicy_Hydro_MC_impl
  : public HyperbolicPolicy_State_Hydro_MC,
    public HyperbolicPolicy_RiemannSolver_Hydro_hllc_MC,
    public HyperbolicPolicy_Slope_dynamic<HyperbolicPolicy_State_Hydro_MC>,
    public HyperbolicPolicy_BoundaryConditions_Hydro_dynamic
{
private:
  using CellIndex     = ForeachCell::CellIndex;
  using CellMetaData  = ForeachCell::CellMetaData;
  using State = HyperbolicPolicy_State_Hydro_MC;

  using RiemannSolver_t = HyperbolicPolicy_RiemannSolver_Hydro_hllc_MC;
  using Slope_t = HyperbolicPolicy_Slope_dynamic<HyperbolicPolicy_State_Hydro_MC>;
  using BoundaryConditions_t = HyperbolicPolicy_BoundaryConditions_Hydro_dynamic;

public:
  using PrimState = State::PrimState;
  using ConsState = State::ConsState;

  struct Params
  {
    HyperbolicPolicy_Hydro_MC_Params policy_params;
    BoundaryConditions_t::Params bc_params;
    Slope_t::Params slope_params;
  };

  static Params getParams( ConfigMap& configMap )
  {
    return Params{
      .policy_params = HyperbolicPolicy_Hydro_MC_Params::from_configMap(configMap),
      .bc_params = BoundaryConditions_t::getParams(configMap),
      .slope_params = Slope_t::getParams(configMap)
    };
  }

  HyperbolicPolicy_Hydro_MC_impl( const Params& params, const ScalarSimulationData& scalar_data )
  : HyperbolicPolicy_State_Hydro_MC(params.policy_params),
    RiemannSolver_t(params.policy_params),
    Slope_t(params.slope_params),
    BoundaryConditions_t(params.bc_params, scalar_data),
    smallr(params.policy_params.smallr),
    smallp(params.policy_params.smallp),
    negative_rho_count("negative_rho_count"),
    negative_p_count("negative_p_count")
  {}
private:
  real_t smallr, smallp;
  Kokkos::View<int> negative_rho_count, negative_p_count;

public:
  KOKKOS_INLINE_FUNCTION
  constexpr static bool has_postProcess()
  {return true;}

  KOKKOS_INLINE_FUNCTION
  ConsState postProcess( const ConsState &u ) const
  {
    real_t smallr = this->smallr;
    real_t smallp = this->smallp;
    PrimState q = this->consToPrim(u);
    if (q.rho < 0.0) {
      this->negative_rho_count()++; 
      // This inaccurate because of concurrency but it's always > 1 when there is an error 
      // Use atomic_inc if you need accurate results
      //Kokkos::atomic_inc( &this->negative_rho_count() );
      q.rho = smallr;
    }
    if (q.p < 0.0) {
      this->negative_p_count() ++;
      //Kokkos::atomic_inc( &this->negative_p_count() );
      q.p   = smallp;
    }
    ConsState u_pp = this->primToCons(q);
    // std::cout << u_pp.rho_old << std::endl;
    u_pp.flux_x_l /= u_pp.rho_old;
    u_pp.flux_x_r /= u_pp.rho_old;
    u_pp.flux_y_l /= u_pp.rho_old;
    u_pp.flux_y_r /= u_pp.rho_old;
    u_pp.flux_z_l /= u_pp.rho_old;
    u_pp.flux_z_r /= u_pp.rho_old;

    return u_pp;
  }

  void printWarnings() const
  {
    auto negative_rho_count = this->negative_rho_count;
    auto negative_p_count = this->negative_p_count;

    Kokkos::parallel_for( "Print Warnings", 1,
      KOKKOS_LAMBDA( int )
    {
      if( negative_rho_count() > 0 )
        Kokkos::printf( "Negative density detected\n");
      if( negative_p_count() > 0 )
        Kokkos::printf( "Negative pressure detected\n" );
    });
  }
};

using HyperbolicPolicy_Hydro_MC = HyperbolicPolicy_base< HyperbolicPolicy_Hydro_MC_impl >;

} //namespace dyablo
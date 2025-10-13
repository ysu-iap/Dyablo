#pragma once

#include "real_type.h"
#include "kokkos_shared.h"
#include "State_Ops.h"
#include "FieldManager.h"
#include "UserData.h"

namespace dyablo {


/**
 * @brief Structure holding conservative hydrodynamics variables
 **/ 
struct ConsHydroMCState {
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
  
  static std::vector<UserData::FieldAccessor::FieldInfo> getFieldsInfo()
  {
    return  { {"rho",      VarIndex::Irho}, 
              {"e_tot",    VarIndex::Ie_tot},
              {"rho_vx",   VarIndex::Irho_vx},
              {"rho_vy",   VarIndex::Irho_vy},
              {"rho_vz",   VarIndex::Irho_vz},
              {"flux_x_l", VarIndex::Iflux_x_l},
              {"flux_x_r", VarIndex::Iflux_x_r},
              {"flux_y_l", VarIndex::Iflux_y_l},
              {"flux_y_r", VarIndex::Iflux_y_r},
              {"flux_z_l", VarIndex::Iflux_z_l},
              {"flux_z_r", VarIndex::Iflux_z_r},
              {"rho_old",  VarIndex::Irho_old}, };
  }

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
  real_t rho_old = 0;
};

DECLARE_STATE_TYPE( ConsHydroMCState, 12 );
DECLARE_STATE_GET( ConsHydroMCState, 0, rho );
DECLARE_STATE_GET( ConsHydroMCState, 1, e_tot );
DECLARE_STATE_GET( ConsHydroMCState, 2, rho_u );
DECLARE_STATE_GET( ConsHydroMCState, 3, rho_v );
DECLARE_STATE_GET( ConsHydroMCState, 4, rho_w );
DECLARE_STATE_GET( ConsHydroMCState, 5, flux_x_l );
DECLARE_STATE_GET( ConsHydroMCState, 6, flux_x_r );
DECLARE_STATE_GET( ConsHydroMCState, 7, flux_y_l );
DECLARE_STATE_GET( ConsHydroMCState, 8, flux_y_r );
DECLARE_STATE_GET( ConsHydroMCState, 9, flux_z_l );
DECLARE_STATE_GET( ConsHydroMCState, 10, flux_z_r );
DECLARE_STATE_GET( ConsHydroMCState, 11, rho_old );


/**
 * @brief Structure holding primitive hydrodynamics variables
 */
struct PrimHydroMCState {
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
    Irho_old,
  };  
  
  static std::vector<UserData::FieldAccessor::FieldInfo> getFieldsInfo()
  {
    return  { {"rho",     VarIndex::Irho}, 
              {"e_tot",   VarIndex::Ip},
              {"rho_vx",  VarIndex::Iu},
              {"rho_vy",  VarIndex::Iv},
              {"rho_vz",  VarIndex::Iw},
              {"flux_x_l", VarIndex::Iflux_x_l},
              {"flux_x_r", VarIndex::Iflux_x_r},
              {"flux_y_l", VarIndex::Iflux_y_l},
              {"flux_y_r", VarIndex::Iflux_y_l},
              {"flux_z_l", VarIndex::Iflux_z_l},
              {"flux_z_r", VarIndex::Iflux_z_l},
              {"rho_old",  VarIndex::Irho_old} };
  }

  static FieldManager getFieldManager()
  {
    return FieldManager( {VarIndex::Irho, VarIndex::Ip, VarIndex::Iu, VarIndex::Iv, VarIndex::Iw, 
                          VarIndex::Iflux_x_l, VarIndex::Iflux_x_r, VarIndex::Iflux_y_l, VarIndex::Iflux_y_r, VarIndex::Iflux_z_l, VarIndex::Iflux_z_r, VarIndex::Irho_old } );
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
  real_t rho_old = 0;
};

DECLARE_STATE_TYPE( PrimHydroMCState, 12 );
DECLARE_STATE_GET( PrimHydroMCState, 0, rho );
DECLARE_STATE_GET( PrimHydroMCState, 1, p );
DECLARE_STATE_GET( PrimHydroMCState, 2, u );
DECLARE_STATE_GET( PrimHydroMCState, 3, v );
DECLARE_STATE_GET( PrimHydroMCState, 4, w );

DECLARE_STATE_GET( PrimHydroMCState, 5, flux_x_l );
DECLARE_STATE_GET( PrimHydroMCState, 6, flux_x_r );
DECLARE_STATE_GET( PrimHydroMCState, 7, flux_y_l );
DECLARE_STATE_GET( PrimHydroMCState, 8, flux_y_r );
DECLARE_STATE_GET( PrimHydroMCState, 9, flux_z_l );
DECLARE_STATE_GET( PrimHydroMCState, 10, flux_z_r );
DECLARE_STATE_GET( PrimHydroMCState, 11, rho_old );

/**
 * @brief Structure grouping the primitive and conservative hydro state as well
 *        as information on the number of fields to store per state
 */
struct HydroMCState {
  using PrimState = PrimHydroMCState;
  using ConsState = ConsHydroMCState;
  static constexpr size_t N = 12;
};

/**
 * @brief Returns a conservative state at a given cell index in an array
 * 
 * @tparam ndim the number of dimensions
 * @tparam Array_t the type of array where we are looking up
 * @tparam CellIndex the type of cell index used
 * 
 * @param U the array in which we are getting the state
 * @param iCell the index of the cell 
 * @return the hydro state at position iCell in U
 */
template< int ndim, 
          typename Array_t, 
          typename CellIndex >
KOKKOS_INLINE_FUNCTION
void getConservativeState(const Array_t& U, const CellIndex& iCell, ConsHydroMCState &res)
{
  res.rho   = U.at(iCell, ConsHydroMCState::VarIndex::Irho );
  res.e_tot = U.at(iCell, ConsHydroMCState::VarIndex::Ie_tot );
  res.rho_u = U.at(iCell, ConsHydroMCState::VarIndex::Irho_vx );
  res.rho_v = U.at(iCell, ConsHydroMCState::VarIndex::Irho_vy );
  res.rho_w = (ndim == 3 ? U.at(iCell, ConsHydroMCState::VarIndex::Irho_vz ) : 0.0);

  res.flux_x_l = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_l);
  res.flux_x_r = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_r);
  res.flux_y_l = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_l);
  res.flux_y_r = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_r);
  res.flux_z_l = (ndim == 3 ? U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_l ) : 0.0);
  res.flux_z_r = (ndim == 3 ? U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_r ) : 0.0);
  res.rho_old  = U.at(iCell, ConsHydroMCState::VarIndex::Irho_old );
}

/**
 * @brief Returns a primitive hydro state at a given cell index in an array
 * 
 * @tparam ndim the number of dimensions
 * @tparam Array_t the type of array where we are looking up
 * @tparam CellIndex the type of cell index used
 * 
 * @param U the array in which we are getting the state
 * @param iCell the index of the cell 
 * @return the hydro state at position iCell in U
 */
template< int ndim,
          typename Array_t, 
          typename CellIndex >
KOKKOS_INLINE_FUNCTION
void getPrimitiveState(const Array_t& U, const CellIndex& iCell, PrimHydroMCState &res)
{
  res.rho = U.at(iCell, PrimHydroMCState::VarIndex::Irho );
  res.p   = U.at(iCell, PrimHydroMCState::VarIndex::Ip );
  res.u   = U.at(iCell, PrimHydroMCState::VarIndex::Iu );
  res.v   = U.at(iCell, PrimHydroMCState::VarIndex::Iv );
  res.w   = (ndim == 3 ? U.at(iCell, PrimHydroMCState::VarIndex::Iw ) : 0.0);

  res.flux_x_l = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_l);
  res.flux_x_r = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_r);
  res.flux_y_l = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_l);
  res.flux_y_r = U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_r);
  res.flux_z_l = (ndim == 3 ? U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_l ) : 0.0);
  res.flux_z_r = (ndim == 3 ? U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_r ) : 0.0);
  res.rho_old  = U.at(iCell, PrimHydroMCState::VarIndex::Irho_old );
}

/**
 * @brief Stores a primitive hydro state in an array
 * 
 * @tparam ndim the number of dimensions
 * @tparam Array_t the type of array in which the primitive value is stored
 * @tparam CellIndex the type of cell index used
 * 
 * @param U the array where we are storing the state
 * @param iCell the index of cell
 * @param u the value to store in the array
 */
template <int ndim, typename Array_t, typename CellIndex >
KOKKOS_INLINE_FUNCTION
void setPrimitiveState( const Array_t& U, const CellIndex& iCell, PrimHydroMCState u) {
  U.at(iCell, PrimHydroMCState::VarIndex::Irho) = u.rho;
  U.at(iCell, PrimHydroMCState::VarIndex::Irho_old) = u.rho_old;
  U.at(iCell, PrimHydroMCState::VarIndex::Ip) = u.p;
  U.at(iCell, PrimHydroMCState::VarIndex::Iu) = u.u;
  U.at(iCell, PrimHydroMCState::VarIndex::Iv) = u.v;
  
  U.at(iCell, PrimHydroMCState::VarIndex::Iflux_x_l) = u.flux_x_l;
  U.at(iCell, PrimHydroMCState::VarIndex::Iflux_x_r) = u.flux_x_r;
  U.at(iCell, PrimHydroMCState::VarIndex::Iflux_y_l) = u.flux_y_l;
  U.at(iCell, PrimHydroMCState::VarIndex::Iflux_y_r) = u.flux_y_r;

  if (ndim == 3){
    U.at(iCell, PrimHydroMCState::VarIndex::Iw) = u.w;
    U.at(iCell, PrimHydroMCState::VarIndex::Iflux_z_l) = u.flux_z_l;
    U.at(iCell, PrimHydroMCState::VarIndex::Iflux_z_r) = u.flux_z_r;
  }
}

/**
 * @brief Stores a conservative hydro state in an array
 * 
 * @tparam ndim the number of dimensions
 * @tparam Array_t the type of array in which the primitive value is stored
 * @tparam CellIndex the type of cell index used
 * 
 * @param U the array where we are storing the state
 * @param iCell the index of cell
 * @param u the value to store in the array
 */
template <int ndim, typename Array_t, typename CellIndex >
KOKKOS_INLINE_FUNCTION
void setConservativeState( const Array_t& U, const CellIndex& iCell, ConsHydroMCState u) {
  U.at(iCell, ConsHydroMCState::VarIndex::Irho) = u.rho;
  U.at(iCell, ConsHydroMCState::VarIndex::Irho_old) = u.rho_old;
  U.at(iCell, ConsHydroMCState::VarIndex::Ie_tot) = u.e_tot;
  U.at(iCell, ConsHydroMCState::VarIndex::Irho_vx) = u.rho_u;
  U.at(iCell, ConsHydroMCState::VarIndex::Irho_vy) = u.rho_v;

  U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_l) = u.flux_x_l;
  U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_r) = u.flux_x_r;
  U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_l) = u.flux_y_l;
  U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_r) = u.flux_y_r;

  if (ndim == 3){
    U.at(iCell, ConsHydroMCState::VarIndex::Irho_vz) = u.rho_w;
    U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_l) = u.flux_z_l;
    U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_r) = u.flux_z_r;
  }
}

template <int ndim, typename Array_t, typename CellIndex >
KOKKOS_INLINE_FUNCTION
void atomic_add_ConservativeState( const Array_t& U, const CellIndex& iCell, ConsHydroMCState u) {
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Irho), u.rho);
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Irho_old), u.rho_old);
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Ie_tot), u.e_tot);
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Irho_vx), u.rho_u);
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Irho_vy), u.rho_v);

  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_l), u.flux_x_l);
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Iflux_x_r), u.flux_x_r);
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_l), u.flux_y_l);
  Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Iflux_y_r), u.flux_y_r);

  if (ndim == 3){
    Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Irho_vz), u.rho_w);
    Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_l), u.flux_z_l);
    Kokkos::atomic_add(&U.at(iCell, ConsHydroMCState::VarIndex::Iflux_z_r), u.flux_z_r);
  }
}

/**
 * @brief Converts from a hydro conservative state to a hydro primitive state
 * 
 * @tparam ndim the number of dimensions
 * 
 * @param U the initial conservative state
 * @param gamma0 adiabatic index
 * @return the primitive version of U
 */
template<int ndim>
KOKKOS_INLINE_FUNCTION
PrimHydroMCState consToPrim(const ConsHydroMCState &U, real_t gamma0) {
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

/**
 * @brief Converts from a hydro primitive state to a hydro conservative state
 * 
 * @tparam ndim the number of dimensions
 * 
 * @param Q the initial primitive state
 * @param gamma0 adiabatic index
 * @return the conservative version of Q
 */
template<int ndim>
KOKKOS_INLINE_FUNCTION
ConsHydroMCState primToCons(const PrimHydroMCState &Q, real_t gamma0) {
    const real_t Ek = 0.5 * Q.rho * (Q.u*Q.u+Q.v*Q.v+Q.w*Q.w);
    const real_t E  = Ek + Q.p / (gamma0-1.0);
    return {Q.rho, 
            E, 
            Q.rho*Q.u, 
            Q.rho*Q.v, 
            (ndim == 3 ? Q.rho*Q.w : 0.0),
            Q.flux_x_l,
            Q.flux_x_r,
            Q.flux_y_l,
            Q.flux_y_r,
            (ndim == 3 ? Q.flux_z_l : 0.0),
            (ndim == 3 ? Q.flux_z_r : 0.0),
            Q.rho_old};
}

/**
 * @brief Swaps a component in velocity with the X component. 
 *        The Riemann problem is always solved by considering an interface on the 
 *        X-axis. So when solving it for other components, those should be swapped 
 *        before and after solving the Riemann problem.
 *  
 * @param Q (IN/OUT) the primitive MHD state to modify
 * @param comp the component to swap with X
 */
KOKKOS_INLINE_FUNCTION
PrimHydroMCState swapComponents(const PrimHydroMCState &q, ComponentIndex3D comp) {
  switch( comp )
  {
    case IX:
      return q;
    case IY:
      return PrimHydroMCState{q.rho, q.p, q.v, q.u, q.w, q.flux_y_l, q.flux_y_r, q.flux_x_l, q.flux_x_r, q.flux_z_l, q.flux_z_r, q.rho_old};
    case IZ:
      return PrimHydroMCState{q.rho, q.p, q.w, q.v, q.u, q.flux_z_l, q.flux_z_r, q.flux_y_l, q.flux_y_r, q.flux_x_l, q.flux_x_r, q.rho_old};
    default:
      DYABLO_ASSERT_KOKKOS_DEBUG(false, "invalid component");
      return PrimHydroMCState{};
  }
}

/**
 * @brief Swaps a component in velocity with the X component. 
 *        The Riemann problem is always solved by considering an interface on the 
 *        X-axis. So when solving it for other components, those should be swapped 
 *        before and after solving the Riemann problem.
 *  
 * @param Q (IN/OUT) the primitive MHD state to modify
 * @param comp the component to swap with X
 */
KOKKOS_INLINE_FUNCTION
ConsHydroMCState swapComponents(const ConsHydroMCState &u, ComponentIndex3D comp) {
  switch( comp )
  {
    case IX:
      return u;
    case IY:
      return ConsHydroMCState{u.rho, u.e_tot, u.rho_v, u.rho_u, u.rho_w, u.flux_y_l, u.flux_y_r, u.flux_x_l, u.flux_x_r, u.flux_z_l, u.flux_z_r, u.rho_old};
    case IZ:
      return ConsHydroMCState{u.rho, u.e_tot, u.rho_w, u.rho_v, u.rho_u, u.flux_z_l, u.flux_z_r, u.flux_y_l, u.flux_y_r, u.flux_x_l, u.flux_x_r, u.rho_old};
    default:
      DYABLO_ASSERT_KOKKOS_DEBUG(false, "invalid component");
      return ConsHydroMCState{};
  }
}

} // namespace dyablo


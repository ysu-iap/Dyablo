#include "ParticleUpdate_base.h"
#include "ForeachParticle.h"
#include <Kokkos_Random.hpp>
namespace dyablo {

class ParticleUpdate_tracers_move_MC : public ParticleUpdate {
  using RNGPool = Kokkos::Random_XorShift64_Pool<>;
  using RNGType = RNGPool::generator_type;
  
  const int    seed;       // Seed number for the random perturbations
  RNGPool      rand_pool;  // Random pool for multi-mode perturbation

  public:
  ParticleUpdate_tracers_move_MC(
          ConfigMap& configMap,
          ForeachCell& foreach_cell,
          Timers& timers) 
  : foreach_cell(foreach_cell),
    foreach_particle(foreach_cell.get_amr_mesh(), configMap),
    timers(timers),
    data{
      .xmin = configMap.getValue<real_t>("mesh", "xmin", 0.0),
      .xmax = configMap.getValue<real_t>("mesh", "xmax", 1.0),      
      .ymin = configMap.getValue<real_t>("mesh", "ymin", 0.0),
      .ymax = configMap.getValue<real_t>("mesh", "ymax", 1.0),
      .zmin = configMap.getValue<real_t>("mesh", "zmin", 0.0),
      .zmax = configMap.getValue<real_t>("mesh", "zmax", 1.0),
      .ndim = configMap.getValue<int>("mesh", "ndim", 3),
    },
    seed(12345),
    rand_pool(seed*GlobalMpiSession::get_comm_world().MPI_Comm_rank()+1)
  {}

  ~ParticleUpdate_tracers_move_MC() {}

  void update( UserData& U, ScalarSimulationData& scalar_data) 
  {
    timers.get("ParticleUpdate_tracers_move_MC").start();

    const real_t dt = scalar_data.get<real_t>("dt");

    enum VarIndex_rhov{
      ID,IVX,IVY,IVZ,IFLUX_X_L,IFLUX_X_R,IFLUX_Y_L,IFLUX_Y_R,IFLUX_Z_L,IFLUX_Z_R
    };

    auto Uin = U.getAccessor( {{"rho", ID}, {"rho_vx", IVX},{"rho_vy", IVY},{"rho_vz", IVZ}, 
                                {"flux_x_l", IFLUX_X_L}, {"flux_x_r", IFLUX_X_R},
                                {"flux_y_l", IFLUX_Y_L}, {"flux_y_r", IFLUX_Y_R},
                                {"flux_z_l", IFLUX_Z_L}, {"flux_z_r", IFLUX_Z_R}});
    const ForeachParticle::ParticleArray& P = U.getParticleArray( "particles" );

    ForeachCell::CellMetaData cells = foreach_cell.getCellMetaData();

    foreach_particle.foreach_particle( "tracers_update_position_MC", P,
      KOKKOS_LAMBDA( const ForeachParticle::ParticleIndex& iPart )
    {
      const real_t x = P.pos(iPart, IX);
      const real_t y = P.pos(iPart, IY);
      const real_t z = P.pos(iPart, IZ);

      ForeachCell::CellIndex iCell = cells.getCellFromPos( {x,y,z} );

      real_t rho = Uin.at( iCell, ID );

      // std::cout << rho << std::endl;
      // P.pos(iPart, IX) += dt * Uin.at( iCell, IVX )/rho;
      // P.pos(iPart, IY) += dt * Uin.at( iCell, IVY )/rho;
      // P.pos(iPart, IZ) += dt * Uin.at( iCell, IVZ )/rho;

      // std::cout << "Je suis là 2" << std::endl;

      // Get mass of present cell
      auto sz = cells.getCellSize(iCell);
      const real_t Ax = sz[0];
      const real_t Ay = sz[1];
      const real_t Az = sz[2];
      
      // Get value of flux in each direction
      real_t fxm = Uin.at(iCell, IFLUX_X_L);
      real_t fxp = Uin.at(iCell, IFLUX_X_R);
      real_t fym = Uin.at(iCell, IFLUX_Y_L);
      real_t fyp = Uin.at(iCell, IFLUX_Y_R);
      real_t fzm = Uin.at(iCell, IFLUX_Z_L);
      real_t fzp = Uin.at(iCell, IFLUX_Z_R);

      fxm = -FMIN(fxm, 0);
      fxp = -FMIN(fxp, 0);
      fym = -FMIN(fym, 0);
      fyp = -FMIN(fyp, 0);
      fzm = -FMIN(fzm, 0);
      fzp = -FMIN(fzp, 0);

      real_t p_out = (fxm + fxp + fym + fyp);

      RNGType rand_gen = rand_pool.get_state();

      double r1 = rand_gen.drand(0., 1.);
      double r2 = rand_gen.drand(0., 1.);

      rand_pool.free_state(rand_gen);

      if ((0 > p_out) || (p_out > 1)) {std::cerr << "something bizarre happen for p_out, p_out vaut " << p_out << std::endl ; p_out = 1;}
      if (r1 >= p_out) return;

      struct FaceProb { int dir; int sign; real_t mass; };
      // std::cout << test << std::endl;
      
      FaceProb faces[6] = {
        {0,-1, fxm}, {0, 1, fxp},
        {1,-1, fym}, {1, 1, fyp},
        {2,-1, fzm}, {2, 1, fzp}
      };
      const int nfaces = 6;

      const Data&d = this->data;
      // std::cout << "nombre de face :" << nfaces << std::endl;
      // std::cout << "ndim :" << d.ndim << std::endl;
      
      int chosen = -1;
      real_t sum = p_out;
      for (int k=0; k<nfaces; ++k) {
        if (faces[k].mass < 0) {std::cerr << "something bizarre happen for face mass" << std::endl; continue;}

        real_t pk = faces[k].mass / sum;
        if (r2 < pk) { chosen = k; break; }
        r2 -= pk;
      }
      if (chosen < 0) std::cerr << "something bizarre happen for chosen variable" << std::endl;

      ForeachCell::CellIndex::offset_t off{0,0,0};
      off[faces[chosen].dir] = faces[chosen].sign;

      // std::cout << off[0] << " " << off[1] << " " << off[2] << std::endl;

      // auto shape = Uin.getShape();
      // ForeachCell::CellIndex iCellN = iCell.getNeighbor_ghost(off, shape);
      // std::cout << "Je suis là\n\n\n\n\n" << std::endl;

      // if (iCellN.is_valid()) {
        // std::cout << "output\n\n\n\n" << std::endl;
        // auto ctr = cells.getCellCenter(iCellN);
        P.pos(iPart, IX) += Ax * off[0];
        P.pos(iPart, IY) += Ay * off[1];
        P.pos(iPart, IZ) += Az * off[2];
        

        P.pos(iPart, IX) = fmod( (P.pos(iPart, IX) - d.xmin) + (d.xmax-d.xmin) , d.xmax-d.xmin) + d.xmin;
        P.pos(iPart, IY) = fmod( (P.pos(iPart, IY) - d.ymin) + (d.ymax-d.ymin) , d.ymax-d.ymin) + d.ymin;
        P.pos(iPart, IZ) = fmod( (P.pos(iPart, IZ) - d.zmin) + (d.zmax-d.zmin) , d.zmax-d.zmin) + d.zmin;
      // }
    });   

    timers.get("ParticleUpdate_tracers_move_MC").stop();
  }

private:
  ForeachCell& foreach_cell;
  ForeachParticle foreach_particle;
  Timers& timers;
public:
  struct Data {
    real_t xmin,xmax,ymin,ymax,zmin,zmax;
    int ndim;
  } data;
};

} // namespace dyablo

FACTORY_REGISTER( dyablo::ParticleUpdateFactory, 
                  dyablo::ParticleUpdate_tracers_move_MC, 
                  "ParticleUpdate_tracers_move_MC")

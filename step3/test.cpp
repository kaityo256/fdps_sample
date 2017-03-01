//------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <string>
#include <fstream>
#include <random>
#include <particle_simulator.hpp>
//------------------------------------------------------------------------
class FPLJ {
public:
  PS::F64vec p;
  PS::F64vec q;
  PS::F64vec getPos() const {
    return q;
  }
  void setPos(PS::F64vec nq) {
		q = nq;
  }
};
//------------------------------------------------------------------------
template<class System>
void
make_conf(System &sys, const PS::F64 L) {
  std::mt19937 mt(1);
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  const double V0 = 1.0;
  const int il = static_cast<int>(L);
  PS::S32 n_total = L * L * L;
  sys.setNumberOfParticleLocal(n_total);
  int n = 0;
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      for (int iz = 0; iz < L; iz++) {
        double z = ud(mt) * 2.0 - 1.0;
        double phi = ud(mt) * M_PI;
        sys[n].p.x = V0 * sqrt(1 - z * z) * cos(phi);
        sys[n].p.y = V0 * sqrt(1 - z * z) * sin(phi);
        sys[n].p.z = V0 * z;
        sys[n].q.x = ix + 0.5;
        sys[n].q.y = iy + 0.5;
        sys[n].q.z = iz + 0.5;
        n++;
      }
    }
  }
}
//------------------------------------------------------------------------
template<class System>
void
dump(const int rank, System &sys) {
	static int step = 0;
	char filename[256];
	sprintf(filename,"dump%02d_%03d.dat",rank,step);
	step++;
  std::ofstream ofs(filename);
  const PS::S32 n = sys.getNumberOfParticleLocal();
  for (PS::S32 i = 0; i < n ; i++) {
    ofs << sys[i].q.x << " ";
    ofs << sys[i].q.y << " ";
    ofs << sys[i].q.z << std::endl;
  }
}
//------------------------------------------------------------------------
template<class Tpsys>
void 
drift(Tpsys & system, const PS::F64 dt) {
	const PS::S32 n = system.getNumberOfParticleLocal();
	for(PS::S32 i = 0; i < n; i++) {
		system[i].q  += system[i].p * dt;
  }
}
//------------------------------------------------------------------------
int
main(int argc, char **argv) {
  PS::Initialize(argc, argv);
  const int rank = PS::Comm::getRank();
  const int procs = PS::Comm::getNumberOfProc();
  PS::ParticleSystem<FPLJ> lj_system;
  lj_system.initialize();

  const PS::F64 L = 10;
  PS::DomainInfo dinfo;
  dinfo.initialize();
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0, 0, 0), PS::F64vec(L, L, L));
  if (rank == 0) {
    make_conf(lj_system, L);
  } else {
    lj_system.setNumberOfParticleLocal(0);
  }
  const PS::F64 dt = 0.01;
  const int LOOP = 1000; 
	const int INTERVAL = 10;
  dinfo.decomposeDomainAll(lj_system);
  lj_system.exchangeParticle(dinfo);
	for(int i=0;i<LOOP;i++){
		if(i%INTERVAL == 0){
      dump(rank, lj_system);
      printf("(%d/%d) n= %d\n", rank, procs, lj_system.getNumberOfParticleLocal());
		}
		drift(lj_system, dt);
		lj_system.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(lj_system);
    lj_system.exchangeParticle(dinfo);
	}
  PS::Comm::barrier();
  PS::Finalize();
}
//------------------------------------------------------------------------

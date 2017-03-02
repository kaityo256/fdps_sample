//------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <string>
#include <fstream>
#include <random>
#include <particle_simulator.hpp>
//------------------------------------------------------------------------
const PS::F64 CUTOFF_LENGTH = 3.0;
//------------------------------------------------------------------------
class ForceLJ {
public:
  PS::F64vec force;
  PS::F64 potential;
  void clear() {
    force = 0.0;
    potential = 0.0;
  }
};
//------------------------------------------------------------------------
struct Atom {
  PS::F64vec p;
  PS::F64vec q;
};
//------------------------------------------------------------------------
class FPLJ {
public:
  PS::S32 id;
  PS::F64vec p;
  PS::F64vec q;
  PS::F64vec force;
  PS::F64 potential;
  PS::F64vec getPos() const {
    return q;
  }
  void setPos(PS::F64vec nq) {
    q = nq;
  }
  void copyFromForce(const ForceLJ & f) {
    force = f.force;
    potential = f.potential;
  }
};
//------------------------------------------------------------------------
class EPLJ {
public:
  PS::S32 id;
  PS::F64vec q;
  PS::F64vec getPos() const {
    return q;
  }
  void setPos(const PS::F64vec & nq) {
    q = nq;
  }
  PS::F64 getRSearch() const {
    return CUTOFF_LENGTH;
  }
  void copyFromFP(const FPLJ & fp) {
    id = fp.id;
    q = fp.q;
  }
};
//------------------------------------------------------------------------
class CalcForceLJ {
public:
  void operator() (const EPLJ * epi, const PS::S32 ni,
                   const EPLJ * epj, const PS::S32 nj,
                   ForceLJ *force) {
    const PS::F64 CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
    const PS::F64 CL6 = CL2 * CL2 * CL2;
    const PS::F64 CL12 = CL6 * CL6;
    const PS::F64 C0 = - 4.0 * (1.0 / CL12 - 1.0 / CL6);

    for (PS::S32 i = 0; i < ni; i++) {
      force[i].force = 0.0;
      force[i].potential = 0.0;
      for (PS::S32 j = 0; j < nj; j++) {
        if (epi[i].id == epj[j].id) continue;
        PS::F64vec dq = epj[j].q - epi[i].q;
        PS::F64 r2 = dq * dq;
        PS::F64 r6 = r2 * r2 * r2;
        PS::F64 r12 = r6 * r6;
        if (r2 > CL2) continue;
        PS::F64 df = (24.0 * r6 - 48.0) / (r6 * r6 * r2);
        force[i].force += dq * df;
        force[i].potential += (4.0 * (1.0 / r12 - 1.0 / r6) + C0) * 0.5;
      }
    }
  }
};
//------------------------------------------------------------------------
template<class Tpsys>
PS::F64
energy(const Tpsys &system) {
  PS::F64 e = 0.0;
  const PS::S32 n = system.getNumberOfParticleLocal();
  for (PS::S32 i = 0; i < n; i++) {
    FPLJ a = system[i];
    e += (a.p * a.p) * 0.5;
    e += a.potential;
  }
  return PS::Comm::getSum(e);
}
//------------------------------------------------------------------------
template<class Tpsys>
void
drift(Tpsys & system, const PS::F64 dt) {
  const PS::S32 n = system.getNumberOfParticleLocal();
  for (PS::S32 i = 0; i < n; i++) {
    system[i].q  += system[i].p * dt;
  }
}
//------------------------------------------------------------------------
template<class Tpsys>
void
kick(Tpsys & system, const PS::F64 dt) {
  const PS::S32 n = system.getNumberOfParticleLocal();
  for (PS::S32 i = 0; i < n; i++) {
    system[i].p  += system[i].force * dt;
  }
}
//------------------------------------------------------------------------
template<class C>
void
add_atom(PS::F64vec q, PS::F64vec c, PS::F64 r, PS::F64 xv, C &atoms) {
  Atom a;
  a.q = q;
  a.p = PS::F64vec(xv, 0.0, 0.0);
  PS::F64vec dq = q - c;
  PS::F64 r2 = dq * dq;
  if (r2 > r * r) return;
  atoms.push_back(a);
}
//------------------------------------------------------------------------
template<class C>
void
put_ball(PS::F64vec c, PS::F64 r, PS::F64 xv, double density, C &atoms) {
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const int is = static_cast<int>(2 * r / s);
  const double hs = s * 0.5;
  PS::F64vec dx(hs, 0, 0);
  PS::F64vec dy(0, hs, 0);
  PS::F64vec dz(0, 0, hs);
  for (int iz = 0; iz < is; iz++) {
    for (int iy = 0; iy < is; iy++) {
      for (int ix = 0; ix < is; ix++) {
        PS::F64vec b(ix * s - r, iy * s - r, iz * s - r);
        add_atom(b + c, c, r, xv, atoms);
        add_atom(b + c + dx, c, r, xv, atoms);
        add_atom(b + c + dy, c, r, xv, atoms);
        add_atom(b + c + dz, c, r, xv, atoms);
      }
    }
  }
}
//------------------------------------------------------------------------
template<class C>
void
export_cdview(int rank, int procs, C &system) {
  static int index = 0;
  char filename[256];
  sprintf(filename, "conf%03d.cdv", index);
  if (0 == rank) {
    std::ofstream ofs(filename);
    std::cout << filename << std::endl;
    ofs.close();
  }
  for (int r = 0; r < procs; r++) {
    PS::Comm::barrier();
    if (rank == r) {
      std::ofstream ofs(filename, std::ios::app);
      const PS::S32 n = system.getNumberOfParticleLocal();
      for (PS::S32 i = 0; i < n; i++) {
        ofs << i << " ";
        ofs << rank << " ";
        ofs << system[i].q.x << " ";
        ofs << system[i].q.y << " ";
        ofs << system[i].q.z << " ";
        ofs << std::endl;
      }
    }
  }
  index++;
}
//------------------------------------------------------------------------
template<class C>
int
make_conf(double density, PS::F64 L, int rank, C &lj_system) {
  std::vector<Atom> atoms;
  put_ball(PS::F64vec(-L * 0.5, 0, 0), L * 0.25, 2, density, atoms);
  put_ball(PS::F64vec(L * 0.5, 0, 0), L * 0.25, -2, density, atoms);
  const int N = atoms.size();
  if (rank == 0) {
    lj_system.setNumberOfParticleLocal(N);
    for (int i = 0; i < atoms.size() ; i++) {
      lj_system[i].id = i;
      lj_system[i].p = atoms[i].p;
      lj_system[i].q = atoms[i].q;
    }
  } else {
    lj_system.setNumberOfParticleLocal(0);
  }
  return N;
}
//------------------------------------------------------------------------
int
main(int argc, char **argv) {
  PS::Initialize(argc, argv);
  const int rank = PS::Comm::getRank();
  const int procs = PS::Comm::getNumberOfProc();
  PS::ParticleSystem<FPLJ> lj_system;
  lj_system.initialize();
  const PS::F64 L = 40;
  PS::DomainInfo dinfo;
  dinfo.initialize();
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(-L, -L*0.5, -L*0.5), PS::F64vec(L, L*0.5, L*0.5));
  PS::TreeForForceShort<ForceLJ, EPLJ, EPLJ>::Gather tree_lj;

  double density = 0.397;
  const int N = make_conf(density, L, rank, lj_system);
  tree_lj.initialize(N);
  const PS::F64 dt = 0.002;
  PS::F64 s_time = 0.0;
  const int LOOP = 10000;
  const int INTERVAL = 100;
  dinfo.decomposeDomainAll(lj_system);
  lj_system.exchangeParticle(dinfo);
  export_cdview(rank, procs, lj_system);

  for (int i = 0; i < LOOP; i++) {
    drift(lj_system, dt * 0.5);
    tree_lj.calcForceAllAndWriteBack(CalcForceLJ(), lj_system, dinfo);
    kick(lj_system, dt);
    drift(lj_system, dt * 0.5);
    lj_system.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(lj_system);
    lj_system.exchangeParticle(dinfo);
    if (i % INTERVAL == 0) {
      export_cdview(rank, procs, lj_system);
    }
    s_time += dt;
  }
  PS::Finalize();
}
//------------------------------------------------------------------------

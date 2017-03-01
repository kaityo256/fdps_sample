//------------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <particle_simulator.hpp>
int
main(int argc, char **argv) {
  PS::Initialize(argc, argv);
  std::cout << std::flush;
  PS::Comm::barrier();
  printf("%d/%d\n", PS::Comm::getRank(), PS::Comm::getNumberOfProc());
  PS::Finalize();
}
//------------------------------------------------------------------------

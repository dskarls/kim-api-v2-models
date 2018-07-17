#include "helper.hpp"

void ProcessVirialTerm(double const dEidr, double const rij,
    double const* const r_ij, int const i, int const j, VectorOfSizeSix virial)
{
  double const v = dEidr / rij;

  virial[0] += v * r_ij[0] * r_ij[0];
  virial[1] += v * r_ij[1] * r_ij[1];
  virial[2] += v * r_ij[2] * r_ij[2];
  virial[3] += v * r_ij[1] * r_ij[2];
  virial[4] += v * r_ij[0] * r_ij[2];
  virial[5] += v * r_ij[0] * r_ij[1];
}


void ProcessParticleVirialTerm(double const dEidr, double const rij,
    double const* const r_ij, int const i, int const j,
    VectorOfSizeSix* const particleVirial)
{
  double const v = dEidr / rij;
  VectorOfSizeSix vir;

  vir[0] = 0.5 * v * r_ij[0] * r_ij[0];
  vir[1] = 0.5 * v * r_ij[1] * r_ij[1];
  vir[2] = 0.5 * v * r_ij[2] * r_ij[2];
  vir[3] = 0.5 * v * r_ij[1] * r_ij[2];
  vir[4] = 0.5 * v * r_ij[0] * r_ij[2];
  vir[5] = 0.5 * v * r_ij[0] * r_ij[1];

  for (int k = 0; k < 6; ++k) {
    particleVirial[i][k] += vir[k];
    particleVirial[j][k] += vir[k];
  }
}

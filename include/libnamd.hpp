#ifndef LIBNAMD_hPP
#define LIBNAMD_hPP

#include <complex>
#include <Eigen/Core>
namespace constants {
  constexpr auto IM = std::complex<double>(0.0, 1.0);
}

namespace rhbi {
namespace libnamd {
/****
 * The following functions are for the general operator algebra.
 ****/
template <typename OP_MatrixType>
Eigen::MatrixXcd commutator(
    const OP_MatrixType &op1,
    const OP_MatrixType &op2) {
  return op1 * op2 - op2 * op1;
}

template <typename OP_MatrixType>
double expectation_value(
    const Eigen::MatrixXcd &rho,
    const OP_MatrixType &op) {
  return (rho * op).trace().real();
}

template <typename OP_VectorType>
double expectation_value(
    const Eigen::MatrixXcd &rho,
    const Eigen::VectorXcd &op) {
    return rho.diagonal().dot(op).real();
}

template <typename OP_MatrixType>
double expectation_value(
    const Eigen::VectorXcd &psi,
    const OP_MatrixType &op) {
    return (psi.adjoint() * op * psi).real();
}

template <typename OP_VectorType>
double expectation_value(
    const Eigen::VectorXcd &psi,
    const Eigen::VectorXcd &op) {
    return psi.dot(op).real();
}

/****
 * The following functions are for the time evoluation
 * of the density matrix and wavefunction.
 ****/
template <typename DC_MatrixType>
void rhs_density_matrix(
    const Eigen::MatrixXcd &rho,
    const Eigen::VectorXd &evals,
    const DC_MatrixType &v_dot_d,
    Eigen::MatrixXcd &drho_dt) {
  const size_t dim = rho.rows();
  for (size_t kk = 0; kk < dim; ++kk) {
    for (size_t jj = 0; jj < dim; ++jj) {
      drho_dt(kk, jj) += -constants::IM * (evals(kk) - evals(jj)) * rho(kk, jj);
      for (size_t ll = 0; ll < dim; ++ll) {
        drho_dt(kk, jj) += (v_dot_d(kk, ll) * rho(ll, jj) - rho(kk, ll) * v_dot_d(ll, jj));
      }
    }
  }
}

template <typename DC_MatrixType>
Eigen::MatrixXcd rhs_density_matrix(
    const Eigen::MatrixXcd &rho,
    const Eigen::VectorXcd &evals,
    const DC_MatrixType &v_dot_d) {
  Eigen::MatrixXcd drho_dt = Eigen::MatrixXcd::Zero(rho.rows(), rho.cols());
  rhs_density_matrix(rho, evals, v_dot_d, drho_dt);
  return drho_dt;
}

template <typename DC_MatrixType>
void rhs_wavefunction(
    const Eigen::VectorXcd &psi,
    const Eigen::VectorXcd &evals,
    const DC_MatrixType &v_dot_d,
    Eigen::VectorXcd &dpsi_dt) {
  const size_t dim = psi.size();
  for (size_t kk = 0; kk < dim; ++kk) {
    dpsi_dt(kk) += -constants::IM * evals(kk) * psi(kk);
    for (size_t ll = 0; ll < dim; ++ll) {
      dpsi_dt(kk) += -v_dot_d(kk, ll) * psi(ll);
    }
  }
}

template <typename DC_MatrixType>
Eigen::VectorXcd rhs_wavefunction(
    const Eigen::VectorXcd &psi,
    const Eigen::VectorXcd &evals,
    const DC_MatrixType &v_dot_d) {
  Eigen::VectorXcd dpsi_dt = Eigen::VectorXcd::Zero(psi.size());
  rhs_wavefunction(psi, evals, v_dot_d, dpsi_dt);
  return dpsi_dt;
}

}  // namespace libnamd
} // namespace rhbi

#endif  // LIBNAMD_HPP

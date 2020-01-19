/**
 *  Berelyium0 Driver
 */

/**
 *  Auxiliary Structure
 */ 


#include <gqcp/gqcp.hpp>
#include <queue>
struct Pair {
    double coeff;
    size_t index;
}

class myComparator 
{ 
public: 
    int operator() (const Pair& p1, const Pair& p2) 
    { 
        return p1.coeff > p2.coeff; 
    } 
}; 
  

int main() {

      // MIN HEAP
      std::priority_queue <Pair, std::vector<Pair>, myComparator > pq (1000); 
      std::cout<<pq.size();
      std::ofstream outfile ("berylium0.txt");

      GQCP::Nucleus Be (GQCP::elements::elementToAtomicNumber("Be"), 0, 0, 0);
      GQCP::Molecule Be_mol ({Be}, 0);
      RSpinorBasis<double, GQCP::GTOShell> spinor_basis (BeZero, "aug-cc-pvdz");

      auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, Be_mol);  // in an AO basis
      auto K = sq_hamiltonian.dimension();

      // Create a plain RHF SCF solver and solve the SCF equations
      GQCP::DIISRHFSCFSolver diis_scf_solver (sq_hamiltonian, spinor_basis, Be_mol);
      diis_scf_solver.solve();
      auto rhf = diis_scf_solver.get_solution();

      // Transform the Hamiltonian to the RHF orbital basis
      GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());

      GQCP::ProductFockSpace fock_space (K, Be_mol.numberOfElectrons()/2, Be_mol.numberOfElectrons()/2);  /

      // Create the FCI module
      GQCP::FCI fci (fock_space);
      GQCP::CISolver ci_solver (fci, sq_hamiltonian);

      // Solve Davidson
      GQCP::VectorX<double> initial_g = fock_space.HartreeFockExpansion();
      GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
      ci_solver.solve(davidson_solver_options);

      // Retrieve the eigenvalues
      auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();
      auto energy = fci_davidson_eigenvalue + GQCP::Operator::NuclearRepulsion(Be_mol).value();
      auto wave = ci_solver.makeWaveFunction();
      GQCP::RDMCalculator rdm_calculator(wave);
      GQCP::TwoRDMs<double> t = rdm_calculator.calculate2RDMs();
      GQCP::OneRDMs<double> o = rdm_calculator.calculate1RDMs();

      auto s2 = GQCP::calculateSpinSquared(o, t);

      outfile << "Eigenvalue:"<< energy << std::endl;
      outfile << "<S^2>:"<< s2 << std::endl;

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> saes(o);
      auto U = saes.eigenvectors();
      outfile << "Natural coefficients:"<< rhf.get_C() * U << std::endl;
      wave.basisTransform(U);
      const auto coefficients = wave.get_coefficients();
      for (size_t i = 0; i < fock_space.get_dimension(); i++) {}
      outfile.close();
      BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
      return 0;
   }

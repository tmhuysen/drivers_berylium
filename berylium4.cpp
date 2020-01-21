/**
 *  Berelyium0 Driver
 */

/**
 *  Auxiliary Structure
 */ 

#include <iomanip>      // std::setprecision
#include <gqcp/gqcp.hpp>
#include <queue>
struct Pair {
    double coeff;
    size_t index;
};

class myComparator 
{ 
public: 
    int operator() (const Pair& p1, const Pair& p2) 
    { 
        return std::abs(p1.coeff) >std::abs(p2.coeff); 
    } 
}; 
  

int main() {
   // MIN HEAP
   std::priority_queue <Pair, std::vector<Pair>, myComparator > pq; 
   std::ofstream outfile ("berylium4.txt");
   outfile<<std::setprecision(16);
   GQCP::Nucleus Be (GQCP::elements::elementToAtomicNumber("Be"), 0, 0, 0);
   GQCP::Molecule Be_mol ({Be}, -4);
   GQCP::Molecule Be_mol_HF ({Be}, -8);
   GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (Be_mol, "aug-cc-pvdz");

   auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, Be_mol);  // in an AO basis
   auto K = sq_hamiltonian.dimension();
   outfile << "K : "<< K << std::endl;
   // Create a plain RHF SCF solver and solve the SCF equations
   GQCP::DIISRHFSCFSolver diis_scf_solver (sq_hamiltonian, spinor_basis, Be_mol_HF);
   diis_scf_solver.solve();
   auto rhf = diis_scf_solver.get_solution();

   // Transform the Hamiltonian to the RHF orbital basis
   GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());

   GQCP::ProductFockSpace fock_space (K, Be_mol.numberOfElectrons()/2, Be_mol.numberOfElectrons()/2);  

   // Create the FCI module
   GQCP::FCI fci (fock_space);
   GQCP::CISolver ci_solver (fci, sq_hamiltonian);

   // Solve Davidson
   GQCP::VectorX<double> initial_g = fock_space.HartreeFockExpansion();
   GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
   davidson_solver_options.convergence_threshold = 1.0e-06;  // the tolerance on the norm of the residual vector
   davidson_solver_options.maximum_subspace_dimension = 5;
   ci_solver.solve(davidson_solver_options);

   // Retrieve the eigenvalues
   auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();
   auto energy = fci_davidson_eigenvalue + GQCP::Operator::NuclearRepulsion(Be_mol).value();
   auto wave = ci_solver.makeWavefunction();
   GQCP::RDMCalculator rdm_calculator(wave);
   GQCP::TwoRDMs<double> t = rdm_calculator.calculate2RDMs();
   GQCP::OneRDMs<double> o = rdm_calculator.calculate1RDMs();

   auto s2 = GQCP::calculateSpinSquared(o, t);

   outfile << "Eigenvalue:"<< energy << std::endl;
   outfile << "<S^2>:"<< s2 << std::endl;

   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(o.one_rdm);
   auto U = saes.eigenvectors();
   std::string type_of_coeff = "Natural coefficients";
   auto T = rhf.get_C();
   try {
      wave.basisTransform(U);
      T = T*U;
   } catch (...) {
      type_of_coeff =  "RHF coefficients:";
   }
   outfile << type_of_coeff <<  T << std::endl;
   const auto coefficients = wave.get_coefficients();
   for (size_t i = 0; i < fock_space.get_dimension(); i++) {
      pq.push(Pair {coefficients(i), i});
      if (pq.size() > 1000) {
         pq.pop();
      }
   }

   while (!pq.empty()) {
      const Pair& x = pq.top();
      outfile << x.index << " : "<< x.coeff << std::endl;
      pq.pop();
   }

   outfile.close();
   return 0;
}

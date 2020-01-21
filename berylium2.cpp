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
   std::ofstream outfile ("berylium2_2x.txt");
   outfile<<std::setprecision(16);
   GQCP::Nucleus Be (GQCP::elements::elementToAtomicNumber("Be"), 0, 0, 0);
   GQCP::Molecule Be_mol ({Be}, -2);
   GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (Be_mol, "aug-cc-pvdz");

   auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, Be_mol);  // in an AO basis
   auto K = sq_hamiltonian.dimension();
   outfile << "K:"<< K  << std::endl;
   // Create a plain RHF SCF solver and solve the SCF equations
   GQCP::DIISRHFSCFSolver diis_scf_solver (sq_hamiltonian, spinor_basis, Be_mol);
   diis_scf_solver.solve();
   auto rhf = diis_scf_solver.get_solution();
   outfile << "RHF Orbital energies:"<< std::endl << rhf.get_orbital_energies() << std::endl;

   // Transform the Hamiltonian to the RHF orbital basis
   GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());

   GQCP::ProductFockSpace fock_space (K, Be_mol.numberOfElectrons()/2, Be_mol.numberOfElectrons()/2);  
   const auto& fock_space_alpha = fock_space.get_fock_space_alpha();
   const auto& fock_space_beta = fock_space.get_fock_space_beta();

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
   auto wave = ci_solver.makeWavefunction();
   GQCP::RDMCalculator rdm_calculator(wave);
   GQCP::TwoRDMs<double> t = rdm_calculator.calculate2RDMs();
   GQCP::OneRDMs<double> o = rdm_calculator.calculate1RDMs();

   auto s2 = GQCP::calculateSpinSquared(o, t);

   outfile << "Eigenvalue:"<< energy << std::endl;
   outfile << "<S^2>:"<< s2 << std::endl;

   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(o.one_rdm);
   auto U = saes.eigenvectors(); 
   outfile << "D_a:"<< std::endl << o.one_rdm_aa << std::endl;
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(o.one_rdm_aa);
   outfile << "Y_a:"<< std::endl << saes2.eigenvalues() << std::endl;
   outfile << "D_b:"<< std::endl << o.one_rdm_bb << std::endl;
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes3(o.one_rdm_bb);
   outfile << "Y_b:"<< std::endl << saes3.eigenvalues() << std::endl;

   outfile << "D_d:"<< std::endl << o.spinDensityRDM() << std::endl;
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes4( o.spinDensityRDM() );
   outfile << "Y_d:"<< std::endl << saes4.eigenvalues() << std::endl;

   outfile << "D:"<< std::endl << o.one_rdm << std::endl;
   outfile << "Y:"<< std::endl << saes.eigenvalues() << std::endl;

   outfile << "RHF coefficients:"<< std::endl << rhf.get_C()<< std::endl;
   outfile << "T coefficients:"<< std::endl << rhf.get_C() * U << std::endl;

   wave.basisTransform(GQCP::TransformationMatrix<double>(U));
   const auto coefficients2 = wave.get_coefficients();
   for (size_t i = 0; i < fock_space.get_dimension(); i++) {
      pq.push(Pair {coefficients2(i), i});
      if (pq.size() > 1000) {
         pq.pop();
      }
   }

   while (!pq.empty()) {
      const Pair& x = pq.top();
      size_t alpha_index = x.index/fock_space_alpha.get_dimension();
      size_t beta_index = x.index - alpha_index * fock_space_alpha.get_dimension();
      const GQCP::ONV alpha = fock_space_alpha.makeONV(alpha_index);
      const GQCP::ONV beta = fock_space_alpha.makeONV(beta_index);
      outfile << x.index << " : " << alpha.asString() << " | " << beta.asString() << " : "<< x.coeff << std::endl;
      pq.pop();
   }

   outfile.close();
   return 0;
}

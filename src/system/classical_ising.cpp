#include "classical_ising.h"
#include "../algorithm/sa.h"
#include "../utility/union_find.h"

#include <cassert>
#include <cmath>

namespace openjij {
    namespace system {

        ClassicalIsing::ClassicalIsing(const graph::Dense<double>& interaction, const graph::Spins& spins)
            : spins(spins), interaction(interaction), urd{0.0, 1.0}{
                //random number generator
                std::random_device rd;
                mt = std::mt19937(rd());
                uid = std::uniform_int_distribution<>{0,(int)spins.size()-1};
                //disable error check to improve performance
                this->interaction.set_err_check(false);
        }

        ClassicalIsing::ClassicalIsing(const graph::Dense<double>& interaction, const graph::Spins& spins, const uint_fast32_t seed)
            : spins(spins), interaction(interaction), urd{0.0, 1.0}{
                //random number generator
                mt = std::mt19937(seed);
                uid = std::uniform_int_distribution<>{0,(int)spins.size()-1};
                //disable error check to improve performance
                this->interaction.set_err_check(false);
        }

        ClassicalIsing::ClassicalIsing(const graph::Dense<double>& interaction)
            :ClassicalIsing(interaction, interaction.gen_spin()) {
        }

        ClassicalIsing::ClassicalIsing(const graph::Dense<double>& interaction, const uint_fast32_t seed)
            :ClassicalIsing(interaction, interaction.gen_spin(seed)) {
        }

        void ClassicalIsing::initialize_spins(){
            spins = interaction.gen_spin();
        }

        void ClassicalIsing::initialize_spins(const uint_fast32_t seed){
            spins = interaction.gen_spin(seed);
        }

        void ClassicalIsing::set_spins(const graph::Spins& initial_spins){
            spins = initial_spins;
        }

        double ClassicalIsing::update(const double beta, const std::string& algo){
            double totaldE = 0;
            const size_t num_spins = spins.size();

            // CASE: single flip
            if(algo == "single_spin_flip" or algo == ""){
                //default updater (single_spin_flip)
                for(auto i = 0u; i < num_spins; ++i){
                    size_t index = uid(mt);
                    //do metropolis
                    double dE = 0;
                    for(auto&& adj_index : interaction.adj_nodes(index)){
                        dE += -2 * spins[index] * (index != adj_index ? (interaction.J(index, adj_index) * spins[adj_index]) : interaction.h(index));
                    }

                    //metropolis
                    if(exp(-beta*dE) > urd(mt)){
                        spins[index] *= -1;
                        totaldE += dE;
                    }
                }
            }
            // CASE: Swendsen Wang
            else if (algo == "swendsen_wang") {
                utility::UnionFindTree clusters(num_spins);
                
                graph::Spins spin_candidates(num_spins);
                for (auto i = 0u; i < num_spins; i++) {
                    spin_candidates[i] = urd(mt) < 0.5 ? -1: 1;
                }
                // make clusters
                const double unite_rate = 1.0 - std::exp(-2.0 * beta);
                for (auto i = 0u; i < num_spins; ++i) {
                    // find adjacent spins
                    for (auto&& adj_index : interaction.adj_nodes(i)) {
                        // are signs of both spins the same?
                        // or
                        // is proberbility equal to or lager than p if the sign of each spin is diffelent?
                        if (spins[i]*spins[adj_index] > 0 || urd(mt) >= unite_rate) continue;
                        // unite clusters of both spins
                        clusters.unite(i, adj_index);
                    }
                }

                // update states in each claster
                for (auto i = 0u; i < num_spins; ++i) {
                    const int cluster = clusters.get_root(i);

                    spins[i] = spin_candidates[cluster];
                }
                totaldE = interaction.calc_energy(spins);
            }

            return totaldE;
        }

        void ClassicalIsing::simulated_annealing(const double beta_min, const double beta_max, const size_t step_length, const size_t step_num, const std::string& algo){
            algorithm::SA sa(beta_min, beta_max, step_length, step_num);
            //do simulated annealing
            sa.run(*this, algo);
        }

        void ClassicalIsing::simulated_annealing(const Schedule& schedule, const std::string& algo) {
            algorithm::SA sa(schedule);
            //do simulated annealing
            sa.run(*this, algo);
        }

        const graph::Spins ClassicalIsing::get_spins() const{
            return spins;
        }
    } // namespace system
} // namespace openjij

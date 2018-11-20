//
// Created by shahnoor on 9/20/18.
//

#include "test_percolation.h"
#include "../types.h"
#include "../percolation/percolation.h"
#include "../util/time_tracking.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace std;

/**
 * This function simulates for all fundamental data for site percolation on square lattice.
 * 1. cluster size distribution by sites and bonds
 * 2. critical point data
 * 3. entropy jump
 * 4. entropy and order parameter
 *
 * length of the lattice and the ensemble size must be provided via command line arguments.
 * @param argc
 * @param argv
 */
void simulate_site_percolation(int argc, char **argv) {
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);

    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v9 lattice_percolation(length, true);

    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ",\"ensemble_size\":" << ensemble_size
                << ",\"signature\":\"" << lattice_percolation.getSignature() << "\""
                << "}" ;

    string tm = getCurrentTime();

    string filename_s = lattice_percolation.getSignature() + "_cluster_by_site_" + to_string(length) + '_' + tm;
    string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond_" + to_string(length) + '_' + tm;
    string filename_critical = lattice_percolation.getSignature() + "_critical_" + to_string(length) + '_' + tm;

    string filename_entropy_order_parameter = lattice_percolation.getSignature()  + to_string(length) + '_' + tm;

    filename_s += ".csv";
    filename_b += ".csv";
    filename_critical += ".csv";
    filename_entropy_order_parameter += ".csv";

    ofstream fout_s(filename_s);
    // JSON formated header
    fout_s << '#' << header_info.str() << endl;
    fout_s << "#each line is an independent realization" << endl;
    fout_s << "#each line contains information about all clusters at critical point" << endl;
    fout_s << "#cluster size is measured by number of sites in it" << endl;

    ofstream fout_b(filename_b);
    // JSON formated header
    fout_b << '#' << header_info.str() << endl;
    fout_b << "#each line is an independent realization" << endl;
    fout_b << "#each line contains information about all clusters at critical point" << endl;
    fout_b << "#cluster size is measured by number of bonds in it" << endl;

    ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << endl;
    fout_critical << "#data at critical occupation probability or pc" << endl;
    fout_critical << "#<pc>,<sites in wrapping cluster>,<bonds in wrapping cluster>" << endl;

    // simulation starts here
    value_type counter{};
    vector<double> entropy(lattice_percolation.maxIterationLimit());
    vector<double> nob_wraping(lattice_percolation.maxIterationLimit()),
            nob_largest(lattice_percolation.maxIterationLimit());

    for(value_type i{} ; i != ensemble_size ; ++i){

        lattice_percolation.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        bool wrapping_written{false};
        while (true){
            successful = lattice_percolation.occupy();
            if(successful) {
                entropy[counter] += lattice_percolation.entropy();
                nob_wraping[counter] += lattice_percolation.numberOfBondsInTheWrappingClusters();
                nob_largest[counter] += lattice_percolation.numberOfBondsInTheLargestCluster_v2();
                lattice_percolation.jump();
                if(!wrapping_written && lattice_percolation.detectWrapping()){
                    fout_critical << lattice_percolation.occupationProbability() << ","
                            << lattice_percolation.numberOfSitesInTheWrappingClusters() << ","
                            << lattice_percolation.numberOfBondsInTheWrappingClusters()  << endl;

                    vector<value_type> site, bond;

                    lattice_percolation.get_cluster_info(site, bond);

                    for(value_type j{}; j != site.size(); ++j){
                        fout_s << site[j] << ',';
                    }
                    for(value_type j{}; j != bond.size(); ++j){
                        fout_b << bond[j] <<',';
                    }


                    fout_s << endl;
                    fout_b << endl;
                    wrapping_written = true;
                }


                ++counter;
            }
            if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }

        auto t_end = std::chrono::system_clock::now();
        cout << "Iteration " << i
             //                 << " . Thread " << std::this_thread::get_id()
             << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;

//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }

    fout_b.close();
    fout_s.close();
    fout_critical.close();



    ofstream fout(filename_entropy_order_parameter);
    fout << '#' << header_info.str() << endl;
    fout << "#<p>,<H(p,L)>,<P1(p,L)>,<P2(p,L)>" << endl;
    fout << "#p = occupation probability" << endl;
    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << endl;
    fout << "#P1(p,L) = Order parameter = (number of bonds in largest cluster) / (total number of bonds)" << endl;
    fout << "#P2(p,L) = Order parameter = (number of bonds in spanning or wrapping cluster) / (total number of bonds)" << endl;
    fout << "#C(p,L) = Specific heat = -T dH/dT" << endl;
    fout << "#X(p,L) = Susceptibility = dP/dp" << endl;
    fout << "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)" << endl;
    for(size_t i{}; i < lattice_percolation.maxIterationLimit(); ++i){
        fout << (i+1) / double(lattice_percolation.maxIterationLimit()) << ",";
        fout << entropy[i] / double(ensemble_size) << ",";
        fout << nob_largest[i] / double(ensemble_size /* lattice_percolation.maxBonds()*/) << ",";
        fout << nob_wraping[i] / double(ensemble_size /* lattice_percolation.maxBonds()*/) ;
        fout << endl;
    }
    fout.close();
}



#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <thread>
#include <mutex>

#include "lattice/lattice.h"
#include "percolation/percolation.h"
#include "util/time_tracking.h"
#include "util/printer.h"

#include "tests/test_percolation.h"


using namespace std;


/****
 *  All the function that is run in main
 * @param argc
 * @param argv
 */
void run_in_main(int argc, char** argv){
    // Simulating site percolation for l=0
    simulate_site_percolation_T<SitePercolation_ps_v9>(argc, argv);

    // Simulating site percolation for l=1
    simulate_site_percolation_T<SitePercolationBallisticDeposition_L1_v2>(argc, argv); // 2018.11.03

    // Simulating site percolation for l=2
    simulate_site_percolation_T<SitePercolationBallisticDeposition_L2_v2>(argc, argv); // 2018.11.03
}




/**************************************
 *  The main function
 *
 ***************************************/
int main(int argc, char** argv) {
    cout << currentTime() << endl;

    cout << "Compiled on " << __DATE__ << "\t at " << __TIME__ << endl;
    std::cout << "Percolation in a Square Lattice" << std::endl;
    auto t_start = std::chrono::system_clock::now();

    time_t seed = time(NULL);
    srand(seed);    // seeding

    run_in_main(argc, argv);

    auto t_end= std::chrono::system_clock::now();
    std::chrono::duration<double> drtion = t_end - t_start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(t_end);
    cout << "Program finished at " << std::ctime(&end_time) << endl;
    std::cout << "Time elapsed "   << getFormattedTime(drtion.count()) << std::endl;
    return 0;
}



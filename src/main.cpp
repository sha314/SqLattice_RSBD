
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

    int l = atoi(argv[1]);
    value_type length = atoi(argv[2]);
    value_type ensemble_size = atoi(argv[3]);

    if(l==1) {
        cout << "Simulating site percolation for l=1" << endl;
        simulate_site_percolation_T<SitePercolationBallisticDeposition_L1_v2>(length, ensemble_size); // 2018.11.03
    }
    else if(l==2) {
        cout << "Simulating site percolation for l=2" << endl;
        simulate_site_percolation_T<SitePercolationBallisticDeposition_L2_v2>(length, ensemble_size); // 2018.11.03
    }else{
        cout << "Simulating site percolation for l=0" << endl;
        simulate_site_percolation_T<SitePercolation_ps_v9>(length, ensemble_size);
    }
}




/**************************************
 *  The main function
 *
 ***************************************/
int main(int argc, char** argv) {

    cout << "Running started at : " << currentTime() << endl;
    cout << "Compiled on        : " << __DATE__ << "\t at " << __TIME__ << endl;
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



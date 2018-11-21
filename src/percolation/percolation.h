//
// Created by shahnoor on 10/5/2017.
//

#ifndef SITEPERCOLATION_PERCOLATION_H
#define SITEPERCOLATION_PERCOLATION_H


#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <climits>
#include <fstream>


#include "../types.h"
#include "../lattice/lattice.h"
#include "../index/index.h"

#include <random>


/**
 * The Square Lattice Percolation class
 */
class SqLatticePercolation{
    // constants
    value_type  _length;
    value_type _max_number_of_bonds;
    value_type _max_number_of_sites;
    char type{'0'}; // percolation type. 's' -> site percolation. 'b' -> bond percolation
protected:

    // structural variables of lattice
    SqLattice _lattice;

    value_type _index_sequence_position{};

    std::vector<Cluster> _clusters;   // check and remove reapeted index manually
    // every birthTime we create a cluster we assign an set_ID for them
    int _cluster_id{};

    value_type min_index{};
    value_type max_index{};

    double _occuption_probability {};
    // entropy
    double _entropy{};
    double _entropy_current{};
    double _entropy_previous{};
    double _largest_jump_entropy{}; // lrgest jump in entropy
    double _entropy_jump_pc{}; // at what pc there have been the largest jump
    size_t _cluster_count{};
    value_type _bonds_in_cluster_with_size_two_or_more{0};   // total number of bonds in the clusters. all cluster has bonds > 1
    bool _reached_critical = false; // true if the system has reached critical value

    value_type _total_relabeling{};
    double time_relabel{};
    value_type _number_of_occupied_sites{};
    value_type _number_of_occupied_bonds{};
    value_type _max_iteration_limit{};
    std::vector<value_type> randomized_index;
    std::random_device _random_device;
    std::mt19937 _random_generator;

    void set_type(char t){type = t;} // setting percolation type
public:
    static constexpr const char* signature = "SqLatticePercolation";

    virtual ~SqLatticePercolation() = default;
    SqLatticePercolation(value_type length);
    void reset();
    void set_cluster_measuring_unit(int i){
        std::cout << "Cluster measuring unit = " << ((i==0) ? "bond" : "site")
                  << " : line " << __LINE__ << std::endl;
    }

    bool occupy();
    value_type length() const { return _length;}
    value_type maxSites() const {return _max_number_of_sites;}
    value_type maxBonds() const { return _max_number_of_bonds;}

    /*********
     *  I/O functions
     */
    virtual void viewCluster();
    virtual void viewClusterExtended();
    virtual void view_bonds(){
        _lattice.view_bonds();
    }
    virtual void viewLattice(){
        _lattice.view_sites();

    }

    /**
 * Also shows the cluster index of the sites
 */
    virtual void viewLatticeExtended(){
        _lattice.view_sites_extended();

    }

    /**
    * Displays group ids of sites in a matrix form
    */
    virtual void viewLatticeByID(){
        _lattice.view_sites_by_id();
        _lattice.view_bonds_by_id();
    }

    virtual void viewSiteByID(){
        _lattice.view_sites_by_id();
    }

    virtual void viewBondByID(){
        _lattice.view_bonds_by_id();
    }

    virtual void viewSiteByRelativeIndex(){
        _lattice.view_sites_by_relative_index();
    }
    virtual void viewBondByRelativeIndex(){
        _lattice.view_bonds_by_relative_index_v4();
    }

    virtual void viewByRelativeIndex(){
        _lattice.view_by_relative_index();
    }

    virtual void view(){
        _lattice.view();
    }

    virtual double occupationProbability() const { return _occuption_probability;}
    virtual double entropy() { return _entropy_current;}
    double entropy_by_site(); // for future convenience. // the shannon entropy. the full calculations. time consuming
    double entropy_by_bond(); // for future convenience. // the shannon entropy. the full calculations. time consuming
    double orderParameter();
    size_t numberOfcluster() const {return _cluster_count;}


    void get_cluster_info(
            std::vector<value_type> &site,
            std::vector<value_type> &bond
    );

    char get_type() const {return type;} // get percolation type
    virtual value_type maxIterationLimit() {return _max_iteration_limit;};

    double get_relabeling_time() const {return time_relabel;}
    value_type relabeling_count() const {return _total_relabeling;}
};


/**
 * Site Percolation by Placing Sites
 *
 * version 9
 *
 * First it randomizes the site index list then use it.
 * Paradigm Shift:
 * Does not delete cluster only makes it empty so that index and id remains the same.
 * This way Searching for index of the cluster using id can be omitted.
 *
 * Feature :
 * 1. Can turn on and off both horizontal and boundary condition
 *
 * 2. Uses class Cluster_v2 for storing clusters
 *
 * 3. Uses Group_ID for Bonds and Sites to identify that they are in the same cluster
 *
 * 4. Occupation probability is calculated by sites,
 *      i.e., number of active sites divided by total number of sites
 *
 * 5. Spanning is calculated by number of bonds in a spanning clusters with periodicity turned off,
 *      i.e., number of bonds in the spanning clusters divided by total number of bonds
 *
 * 6. Unweighted relabeling is ommited in this version ??
 *
 * 7. Runtime is significantly improved. For example, if L=200 program will take ~1 min to place all sites.
 *
 * 8. Unnecessary methods of previous version is eliminated
 *
 * 9. Checking spanning by keeping track of boundary sites is implemented
 *
 * 10. last modified cluster id can be obtained from @var _last_placed_site
 *
 *
 */
class SitePercolation_ps_v9 : public SqLatticePercolation{
protected:
    // flags to manipulate method
    bool _periodicity{false};

    value_type min_index; // minimum index = 0
    value_type max_index; // maximum index = length - 1

    // index sequence
    std::vector<Index> index_sequence;  // initialized once
    std::vector<value_type> randomized_index;

    // every birthTime we create a cluster we assign an set_ID for them
    int _cluster_id{};
    value_type _index_last_modified_cluster{};  // id of the last modified cluster

    // order parameter calculation ingradiants
    // id of the cluster which has maximum number of bonds. used to calculate order parameter
    value_type _number_of_bonds_in_the_largest_cluster{};
    value_type _number_of_sites_in_the_largest_cluster{};   // might be useful later

    Index _last_placed_site;    // keeps track of last placed site

    /**************
     * Spanning variables
     ************/
    /*Holds indices on the edges*/
    std::vector<Index> _top_edge, _bottom_edge, _left_edge, _right_edge;

    std::vector<Index> _spanning_sites;
    std::vector<Index> _wrapping_sites;
    std::vector<value_type> number_of_sites_to_span;
    std::vector<value_type> number_of_bonds_to_span;

    value_type _total_relabeling{};

    /*****************************************
     * Private Methods
     ******************************************/
    void relabel_sites(const std::vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab) ;

    double time_relabel{};
public:
    static constexpr const char* signature = "SitePercolation_ps_v8";

    ~SitePercolation_ps_v9() = default;
    SitePercolation_ps_v9() = default;
    SitePercolation_ps_v9(SitePercolation_ps_v9 & ) = default;
    SitePercolation_ps_v9(SitePercolation_ps_v9 && ) = default;
    SitePercolation_ps_v9(value_type length, bool periodicity=true);

    SitePercolation_ps_v9& operator=(SitePercolation_ps_v9 & ) = default;
//    SitePercolation_ps_v8&& operator=(SitePercolation_ps_v8 && ) = default;
    double get_relabeling_time() {return time_relabel;}
    value_type relabeling_count() const {return _total_relabeling;}

    virtual void reset();



    /*************************************************
     * Flags
     ************************************************/
    bool periodicity() const {return _periodicity;}


    /***********************************************
     * Properties of Percolation class
     ***********************************************/
    std::string getSignature();


    void numberOfActiveSites() const {std::cout << "Number of active sites " << _number_of_occupied_sites << std::endl;}
    double activeSites() const { return _number_of_occupied_sites;}

    value_type count_number_of_active_site();

    int birthTimeOfSpanningCluster() const;
    int birthTimeOfACluster(int id) const;


    /****************************************************************
     * Calculations
     ***************************************************************/

    void add_entropy_for_bond(value_type index);
    void subtract_entropy_for_bond(const std::set<value_type> &found_index_set, int base=-1);

    /*************************************************
     * Site placing methods
     *
     ************************************************/
    virtual bool occupy();
    value_type placeSite(Index site,
                         std::vector<Index>& neighbor_sites,
                         std::vector<BondIndex>& neighbor_bonds);
    value_type placeSite(Index site);
    value_type placeSite_weighted(Index site); // uses weighted relabeling by first identifying the largest cluster
    value_type placeSite_weighted(Index site,
                                  std::vector<Index>& neighbor_sites,
                                  std::vector<BondIndex>& neighbor_bonds);

    Index selectSite(); // selecting site

    void connection_v2(Index site, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);

    /*************************************************
     * Relabeling methods
     *************************************************/
    // applicable to weighted relabeling
    value_type relabel(value_type index_1, value_type index_2);
    void relabel_sites(const Cluster&  clstr, int id);
    void relabel_sites_v4(Index root_a, const Cluster& clstr_b); // relative index is set accordingly
    void relabel_sites_v5(Index root_a, const Cluster& clstr_b); // relative index is set accordingly
    void relabel_sites_v6(Index root_a, const Cluster& clstr_b, int id); // relative index is set accordingly
    void relabel_bonds(const Cluster&  clstr, int id);


    /**********************************************
     * Information about current state of Class
     **********************************************/
    double numberOfOccupiedSite() const { return _number_of_occupied_sites;}
    double occupationProbability() const { return double(_number_of_occupied_sites)/maxSites();}
    double entropy(); // the shannon entropy

    double orderParameter() const;  // number of bonds in the largest cluster / total number of bonds
    double orderParameter_v2() const;  // number of bonds in the largest cluster / total number of bonds

    value_type numberOfBondsInTheLargestCluster();
    value_type numberOfBondsInTheLargestCluster_v2();
    value_type numberOfBondsInTheSpanningCluster();
    value_type numberOfSitesInTheLargestCluster();

    double numberOfSitesInTheSpanningClusters()  ;
    double numberOfBondsInTheSpanningClusters()  ;
    value_type numberOfSitesInTheSpanningClusters_v2()  ;
    value_type numberOfBondsInTheSpanningClusters_v2()  ;

    value_type numberOfSitesInTheWrappingClusters()  ;
    value_type numberOfBondsInTheWrappingClusters()  ;

    value_type wrappingClusterSize() {return numberOfBondsInTheWrappingClusters();}


    /***********************************
     * Spanning Detection
     **********************************/
    bool detectSpanning_v6(const Index& site);

    void save_index_if_in_boundary_v2(const Index& site);
    bool check_if_id_matches(Index site, const std::vector<Index> &edge);
    bool check_if_id_matches_and_erase(Index site, std::vector<Index> &edge);
    bool isSpanned() const { return _reached_critical;}

    bool detectWrapping();

    /************************************
     *  Tracker
     *  Must be called each time a site is placed
     ***********************************/
    void track_numberOfBondsInLargestCluster();
    void track_numberOfSitesInLargestCluster();

    /*********************************
     * I/O functions
     * Printing Status
     ********************************/
    Index lastPlacedSite() const { return _last_placed_site;}

    void spanningIndices() const;
    void wrappingIndices() const;

    /***********************************************
     * Visual data for plotting
     *********************************************/
    // lattice visual data for python
    void writeVisualLatticeData(const std::string& filename, bool only_spanning=true);

protected:
    void initialize();
    void initialize_index_sequence();
    void randomize_v2(); // better random number generator

    std::set<value_type> find_index_for_placing_new_bonds(const std::vector<Index> &neighbors);
    int find_cluster_index_for_placing_new_bonds(const std::vector<Index> &neighbors, std::set<value_type> &found_indices);

    value_type manage_clusters(
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters(
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site,
            int base_id // since id and index is same
    );

    bool anyActiveSite(value_type r, value_type c, value_type delta);
    bool anyActiveSpanningSite(value_type row, value_type col, value_type delta);

public:
    // on test
    IndexRelative getRelativeIndex(Index root, Index site_new);

};


/******************************************************************************
 * Explosive site percolation in square lattice with sum rule and product rule
 */
class SitePercolationExplosive: public SitePercolation_ps_v9{

public:
    ~SitePercolationExplosive() = default;
    SitePercolationExplosive(value_type length);
    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_explosive_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        s += std::to_string(length());
        return s;
    }

};

/*******************************************************************************
 * Site Percolation Ballistic Deposition
 * Extended from SitePercolation_ps_v9
 * *************************************************************/
class SitePercolationBallisticDeposition_v2: public SitePercolation_ps_v9{
protected:
    // elements of @indices_tmp will be erased if needed but not of @indices
    std::vector<value_type> indices;
    std::vector<value_type> indices_tmp;
public:
    static constexpr const char* signature = "SitePercolation_BallisticDeposition_v2";
    virtual ~SitePercolationBallisticDeposition_v2(){
        indices.clear();
        indices_tmp.clear();
    };
    SitePercolationBallisticDeposition_v2(value_type length, bool periodicity);

    virtual bool occupy();

    /************************************
     * Site selection methods
     */
    Index select_site(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_1nn(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_2nn(std::vector<Index> &sites, std::vector<BondIndex> &bonds);


    void reset();
    void initialize_indices();
//    void randomize_index();

    virtual std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }


    /***********************************
     * occupy upto 1st nearset neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied skip the rest setps and start next iteration Else occupy it
     *
     *
     */

    value_type placeSite_1nn_v2();


    /*********************************
     * occupy upto 2nd nearest neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied, select the next neighbor in the direction of motion Else occupy it.
     * If the 2nd nearest neighbor in the direction of motion is also occupied then skip the rest of the steps
     *      and start the next iteration
     *
     */

    value_type placeSite_2nn_v1();

};

/***********
 * Only L1
 */
class SitePercolationBallisticDeposition_L1_v2: public SitePercolationBallisticDeposition_v2{
public:
    ~SitePercolationBallisticDeposition_L1_v2() = default;
    SitePercolationBallisticDeposition_L1_v2(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v2(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {
//        value_type v = placeSite_1nn_v0(); // debugging version
            value_type v = placeSite_1nn_v2();
            _occuption_probability = occupationProbability(); // for super class
            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L1";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};

/*********************
 *
 */
class SitePercolationBallisticDeposition_L2_v2: public SitePercolationBallisticDeposition_v2{
public:
    ~SitePercolationBallisticDeposition_L2_v2() = default;
    SitePercolationBallisticDeposition_L2_v2(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v2(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {

//            value_type v = placeSite_2nn_v0();
            value_type v = placeSite_2nn_v1();
            _occuption_probability = occupationProbability(); // for super class

            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L2";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};

#endif //SITEPERCOLATION_PERCOLATION_H


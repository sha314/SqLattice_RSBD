//
// Created by shahnoor on 9/4/18.
//


#include <cstdlib>
#include <climits>
#include <unordered_set>
#include <mutex>

#include "percolation.h"

#include "../util/printer.h"
#include <omp.h>
#include <thread>
#include <algorithm>

#include "../util/time_tracking.h"

using namespace std;



/**
 *
 * @param length       : length of the lattice
 * @param impure_sites : number of impure sites. cannot be greater than length*length
 */
SitePercolation_ps_v9::SitePercolation_ps_v9(value_type length, bool periodicity)
        :SqLatticePercolation(length)
{
    std::cout << "Constructing SitePercolation_ps_v9 object : line " << __LINE__ << endl;
    SqLatticePercolation::set_type('s');

    _periodicity = periodicity;
    _index_sequence_position = 0;
    _lattice = SqLattice(length, true, false, false, true);   // since it is a site percolation all bonds will be activated by default

    min_index = 0;
    max_index = length - 1;

    index_sequence.resize(maxSites());
    randomized_index.resize(maxSites());
    _max_iteration_limit = maxSites();

    initialize_index_sequence();
    initialize();
    randomize_v2();  // randomize the untouched_site_indices
}


/**
 *
 */
void SitePercolation_ps_v9::initialize() {

    // to improve performence
    number_of_sites_to_span.reserve(maxSites());
    number_of_bonds_to_span.reserve(maxSites());

    _top_edge.reserve(length());
    _bottom_edge.reserve(length());
    _left_edge.reserve(length());
    _right_edge.reserve(length());

//    randomized_index_sequence = index_sequence;
}


/**
 * Called only once when the object is constructed for the first time
 */
void SitePercolation_ps_v9::initialize_index_sequence() {
    value_type m{}, n{};
    for (value_type i{}; i != index_sequence.size(); ++i) {
        randomized_index[i] = i;
        index_sequence[i] = Index(m, n);
        ++n;
        if (n == length()) {
            n = 0;
            ++m;
        }
    }
    //for (value_type i{}; i != index_sequence.size(); ++i) {cout << index_sequence[i] << endl;}
}


/**
 * Reset all calculated values and then call initiate()
 * to initiallize for reuse
 *
 * caution -> it does not erase _calculation_flags, for it will be used for calculation purposes
 */
void SitePercolation_ps_v9::reset() {
    SqLatticePercolation::reset();
    // variables
    _number_of_occupied_sites = 0;
    _index_sequence_position = 0;
//    _first_spanning_cluster_id = -1;
    _cluster_id = 0;
//    bonds_in_largest_cluster = 0;
//    sites_in_largest_cluster = 0;

    // containers
//    randomized_index_sequence.clear();    // reseted in the initialize function

//    _number_of_occupied_sites.clear();
//    _entropy_by_bond.clear();
    number_of_sites_to_span.clear();
    number_of_bonds_to_span.clear();
//    spanning_cluster_ids.clear();
    _spanning_sites.clear();
    _wrapping_sites.clear();
//    wrapping_cluster_ids.clear();

    _bonds_in_cluster_with_size_two_or_more = 0;
//    _id_largest_cluster = 0;

//    _id_last_modified_cluster = -1;

    _index_last_modified_cluster = 0;  // id of the last modified cluster
//    _index_largest_cluster = 0;
    _number_of_bonds_in_the_largest_cluster = 0;
    _number_of_sites_in_the_largest_cluster = 0;
//    _cluster_id_set.clear();

    // clearing edges
    _top_edge.clear();
    _bottom_edge.clear();
    _left_edge.clear();
    _right_edge.clear();

    initialize();
    randomize_v2();
    time_relabel = 0;
    _total_relabeling = 0;
}


/**
 * Randomize the indices
 */
void SitePercolation_ps_v9::randomize_v2(){

    std::shuffle(randomized_index.begin(), randomized_index.end(), _random_generator);
//    cout << "Index sequence : " << randomized_index_sequence << endl;
}


/*************************************************
 * Calculation methods
 *
 ***********************************/

/*
 * Instead of calculating entropy for 1000s of cluster in every iteration
 * just keep track of entropy change, i.e.,
 * how much to subtract and how much to add.
 */
/**
 * Must be called before merging the clusters
 * @param found_index_set
 */
void SitePercolation_ps_v9::subtract_entropy_for_bond(const set<value_type> &found_index, int base){
    double nob, mu_bond, H{};
    if(base >= 0){
        nob = _clusters[base].numberOfBonds();
        mu_bond = nob / maxBonds();
        H += log(mu_bond) * mu_bond;
    }
    for(auto x : found_index){
        nob = _clusters[x].numberOfBonds();
        mu_bond = nob / maxBonds();
        H += log(mu_bond) * mu_bond;
    }
    _entropy -= -H;
}



/**
 * Must be called after merging the clusters
 * Cluster length is measured by bonds
 * @param index
 */
void SitePercolation_ps_v9::add_entropy_for_bond(value_type index){
    double nob = _clusters[index].numberOfBonds();
    double mu_bond = nob / maxBonds();
    double H = log(mu_bond) * mu_bond;
    _entropy += -H;
}



/**
 * Condition: must be called each time a site is placed
 */
void SitePercolation_ps_v9::track_numberOfBondsInLargestCluster() {

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfBonds() > _number_of_bonds_in_the_largest_cluster){
        _number_of_bonds_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfBonds();
    }

}

/**
 *
 */
void SitePercolation_ps_v9::track_numberOfSitesInLargestCluster(){

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfSites() > _number_of_sites_in_the_largest_cluster){
        _number_of_sites_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfSites();
    }
}


/**
 *
 * Find one row from _cluster to place 4 or less new bonds
 * Also remove the matched index values, because they will be inserted later.
 * This gives an advantage, i.e., you don't need to perform a checking.
 * @param hv_bonds
 * @return a set
 */
set<value_type>
SitePercolation_ps_v9::find_index_for_placing_new_bonds(const vector<Index> &neighbors) {
    set<value_type> found_index_set;    // use set to prevent repeated index
    for (auto n: neighbors) {
        int id = _lattice.getSite(n).get_groupID();
        if(id >=0) {
            found_index_set.insert(value_type(id));
        }
    }

    return found_index_set;
}


/**
 *
 * @param neighbors         :
 * @param found_index_set   : index of the clusters that will be merged together.
 *                            Does not contain the base cluster index or id.
 * @return                  : id of the base cluster
 */
int
SitePercolation_ps_v9::find_cluster_index_for_placing_new_bonds(
        const vector<Index> &neighbors, std::set<value_type> &found_index_set
){
    found_index_set.clear();
    value_type size{}, tmp{}, index, base{ULONG_MAX};
    int base_id{-1};
    int id;
    for (auto n: neighbors) {
        id = _lattice.getSite(n).get_groupID();
        if(id >=0) {
            index = value_type(id);
            tmp = _clusters[index].numberOfSites();
            if(tmp > size){
                size = tmp;
                base_id = id;
                base = index;
            }

            found_index_set.insert(index);

        }
    }
    found_index_set.erase(base);
    return base_id;
}


/**
 * Last placed site is added to a cluster. If this connects other clusters then merge all
 * cluster together to get one big cluster. All sites that are part of the other clusters
 * are relabled according to the id of the base cluster.
 * @param found_index_set : index of the clusters that are neighbors of the last placed site
 * @param hv_bonds        : bonds that connects the last placed site and its neighbors
 *                          and which are not part of any cluster of size larger than one
 * @param site            : last placed site
 * @param base_id         : id of the base cluster
 * @return
 */
value_type SitePercolation_ps_v9::manage_clusters(
        const set<value_type> &found_index_set,
        vector<BondIndex> &hv_bonds,
        Index &site,
        int base_id
)
{


    if (base_id != -1) {
        value_type base = value_type(base_id); // converting here
        _clusters[base].addSiteIndex(site);
        int id_base = _clusters[base].get_ID();
        vector<Index> neibhgors = _lattice.get_neighbor_site_indices(site);
        // find which of the neighbors are of id_base as the base cluster
        IndexRelative r;
        for(auto n: neibhgors){
            if(_lattice.getSite(n).get_groupID() == id_base){
                // find relative index with respect to this site
                r = getRelativeIndex(n, site);
                break; // since first time r is set running loop is doing no good
            }
        }

        // put_values_to_the_cluster new values in the 0-th found index
        _clusters[base].insert(hv_bonds);
        _lattice.getSite(site).relativeIndex(r);
        _lattice.getSite(site).set_groupID(id_base); // relabeling for 1 site

        // merge clusters with common values from all other cluster        // merge clusters with common values from all other cluster


        for(value_type ers: found_index_set){

            _total_relabeling += _clusters[ers].numberOfSites(); // only for debugging purposes
            // perform relabeling on the sites
            relabel_sites_v5(site, _clusters[ers]);

            // store values of other found indices to the cluster
            _clusters[base].insert_v2(_clusters[ers]);
            _cluster_count--; // reducing number of clusters
            _clusters[ers].clear(); // emptying the cluster

        }
        _index_last_modified_cluster = base;


    } else {
        // create new element for the cluster
        _clusters.push_back(Cluster(_cluster_id));
        value_type _this_cluster_index = _clusters.size() -1;
        _lattice.getSite(site).set_groupID(_cluster_id); // relabeling for 1 site
        _cluster_count++; // increasing number of clusters
        _cluster_id++;
        _clusters.back().insert(hv_bonds);
        _clusters[_this_cluster_index].addSiteIndex(site);
        _index_last_modified_cluster = _this_cluster_index;   // last cluster is the place where new bonds are placed

    }
    return _index_last_modified_cluster;
}





/**
 * Relative index of site_new with respect to root
 * @param root
 * @param site_new
 * @return
 */
IndexRelative SitePercolation_ps_v9::getRelativeIndex(Index root, Index site_new){
//    cout << "Entry \"SitePercolation_ps_v9::getRelativeIndex\" : line " << __LINE__ << endl;
    int delta_x = -int(root.column_) + int(site_new.column_); // if +1 then root is on the right ??
    int delta_y = int(root.row_) - int(site_new.row_); // if +1 then root is on the top ??


    // normalizing delta_x
    if(delta_x > 1){
        delta_x /= -delta_x;
    }
    else if(delta_x < -1){
        delta_x /= delta_x;
    }

    // normalizing delta_y
    if(delta_y > 1){
        delta_y /= -delta_y;
    }else if(delta_y < -1){
        delta_y /= delta_y;
    }

    IndexRelative indexRelative_root = _lattice.getSite(root).relativeIndex();
//    cout << "Relative index of root " << indexRelative_root << endl;
//    cout << "Delta x,y " << delta_x << ", " << delta_y << endl;
    IndexRelative r =  {indexRelative_root.x_ + delta_x, indexRelative_root.y_ + delta_y};
//    cout << "Relative index of site_new " << r << endl;
    return r;
}



/**
 * Take a bond index only if the corresponding site is active
 * takes longer? time than version 1?, i.e.,  connection()
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void SitePercolation_ps_v9::connection_v2(Index site, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor)
{

    value_type prev_column  = (site.column_ + length() - 1) % length();
    value_type prev_row     = (site.row_ + length() - 1) % length();
    value_type next_row     = (site.row_ + 1) % length();
    value_type next_column  = (site.column_ + 1) % length();

    if(!_periodicity){
        // without periodicity
        if (site.row_ == min_index) { // top edge including corners
            if(site.column_ == min_index){
                // upper left corner

                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {next_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, site.row_, site.column_});
                }

                return;

            }
            else if(site.column_ == max_index){
                // upper right corner

                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, prev_column};
                site_neighbor[1] = {next_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, site.row_, site.column_});
                }

                return;
            }
            else{
                // top edge excluding corners
                site_neighbor.resize(3);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {site.row_, prev_column};
                site_neighbor[2] = {next_row, site.column_};

                bond_neighbor.reserve(4);
                if(!_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[2]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
                }

                return;

            }
        }
        else if (site.row_ == max_index) { // bottom edge including corners
            if (site.column_ == min_index) {
                // lower left corner
                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {prev_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
                }


                return;

            } else if (site.column_ == max_index) {
                // lower right corner
                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, prev_column};
                site_neighbor[1] = {prev_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
                }

                return;

            } else {
                // bottom edge excluding corners
                //  bottom edge
                site_neighbor.resize(3);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {site.row_, prev_column};
                site_neighbor[2] = {prev_row, site.column_};

                bond_neighbor.reserve(3);
                if(!_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[2]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
                }

                return;
            }
        }
            /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
        else if (site.column_ == min_index) { // left edge not in the corners
            site_neighbor.resize(3);
            site_neighbor[0] = {site.row_, next_column};
            site_neighbor[1] = {next_row, site.column_};
            site_neighbor[2] = {prev_row, site.column_};

            bond_neighbor.reserve(3);
            if(!_lattice.getSite(site_neighbor[0]).isActive()) {
                bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
            }
            if(!_lattice.getSite(site_neighbor[1]).isActive()){
                bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
            }
            if(!_lattice.getSite(site_neighbor[2]).isActive()){
                bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
            }

            return;
        }
        else if (site.column_ == max_index) {
            // right edge no corners

            site_neighbor.resize(3);
            site_neighbor[0] = {site.row_, prev_column};
            site_neighbor[1] = {next_row, site.column_};
            site_neighbor[2] = {prev_row, site.column_};

            bond_neighbor.reserve(3);
            if(!_lattice.getSite(site_neighbor[0]).isActive()){
                bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
            }
            if(!_lattice.getSite(site_neighbor[1]).isActive()){
                bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
            }
            if(!_lattice.getSite(site_neighbor[2]).isActive()){
                bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
            }

            return;
        }

    }
    // 1 level inside the lattice
    // not in any the boundary
    site_neighbor.resize(4);
    site_neighbor[0] = {site.row_, next_column};
    site_neighbor[1] = {site.row_, prev_column};
    site_neighbor[2] = {next_row, site.column_};
    site_neighbor[3] = {prev_row, site.column_};

    bond_neighbor.reserve(4);
    if(!_lattice.getSite(site_neighbor[0]).isActive()) {
        bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
    }
    if(!_lattice.getSite(site_neighbor[1]).isActive()){
        bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
    }
    if(!_lattice.getSite(site_neighbor[2]).isActive()){
        bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
    }
    if(!_lattice.getSite(site_neighbor[3]).isActive()) {
        bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
    }

}



/**
 *
 * @param site
 * @param edge
 * @return
 */
bool SitePercolation_ps_v9::check_if_id_matches(Index site, const vector<Index> &edge){
    for(auto s :edge){
        if(_lattice.getSite(site).get_groupID() == _lattice.getSite(s).get_groupID()){
            // no need to put the site here
//            cout << "Site " << site << " and Id " << _lattice.getSite(site).set_groupID()
//                 << " is already in the edge : line " << __LINE__ << endl;
            return true;
        }
    }
    return false;
}



/***********************************************
 *  Placing sites
 *
 *****************************************/

/**
 * All site placing method in one place
 *
 * @return true if operation is successfull
 */
bool SitePercolation_ps_v9::occupy() {
    if(_index_sequence_position >= maxSites()){
        return false;
    }
    Index site = selectSite();
    placeSite_weighted(site);
    _occuption_probability = occupationProbability(); // for super class
    return true;
}

/***
 * Index of the selected site must be provided with the argument
 *
 * Wrapping and spanning index arrangement is enabled.
 * Entropy is calculated smoothly.
 * Entropy is measured by site and bond both.
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v9::placeSite_weighted(Index current_site) {
    // randomly choose a site
    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }

    _last_placed_site = current_site;
    _lattice.activate_site(current_site);
    ++_number_of_occupied_sites;
    // find the bonds for this site
    vector<BondIndex> bonds;
    vector<Index>     sites;
    connection_v2(current_site, sites, bonds);
    _bonds_in_cluster_with_size_two_or_more += bonds.size();

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set;
    int  base_id = find_cluster_index_for_placing_new_bonds(sites, found_index_set);

    subtract_entropy_for_bond(found_index_set, base_id);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters(
            found_index_set, bonds, current_site, base_id
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change
    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    track_numberOfSitesInLargestCluster();
    return merged_cluster_index;
}

/***
 * Index of the selected site must be provided with the argument
 *
 * Wrapping and spanning index arrangement is enabled.
 * Entropy is calculated smoothly.
 * Entropy is measured by site and bond both.
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v9::placeSite_weighted(
        Index current_site,
        vector<Index>& neighbor_sites,
        vector<BondIndex>& neighbor_bonds
) {
    // randomly choose a site
    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }
    _bonds_in_cluster_with_size_two_or_more += neighbor_bonds.size();
    _last_placed_site = current_site;
    _lattice.activate_site(current_site);
    ++_number_of_occupied_sites;
    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set;
    int  base_id = find_cluster_index_for_placing_new_bonds(neighbor_sites, found_index_set);
    subtract_entropy_for_bond(found_index_set, base_id);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters(
            found_index_set, neighbor_bonds, current_site, base_id
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change
    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    track_numberOfSitesInLargestCluster();
    return merged_cluster_index;
}



/**
 *
 * @return
 */
Index SitePercolation_ps_v9::selectSite(){
//    Index current_site = randomized_index_sequence[_index_sequence_position]; // old
    value_type index = randomized_index[_index_sequence_position];
    Index current_site = index_sequence[index]; // new process
    ++_index_sequence_position;
    return current_site;
}


/***************************************************
 * View methods
 ****************************************/


/**
 *
 */
void SitePercolation_ps_v9::spanningIndices() const {
    cout << "Spanning Index : id" << endl;
    for(auto i: _spanning_sites){
        cout << "Index " << i << " : id " << _lattice.getSite(i).get_groupID() << endl;
    }
}

void SitePercolation_ps_v9::wrappingIndices() const {
    cout << "Wrapping Index : id : relative index" << endl;
    for(auto i: _wrapping_sites){
        cout << "Index " << i << " : id "
             << _lattice.getSite(i).get_groupID()
             << " relative index : " << _lattice.getSite(i).relativeIndex() << endl;
    }
}


/****************************************
 * Spanning Detection
 ****************************************/


/**
 * success : gives correct result
 * length       time
 * 200          7.859000 sec
 * 500          2 min 18.874000 sec
 * only check for the cluster id of the recently placed site
 * @param site : Check spanning for this argument
 * @return
 */
bool SitePercolation_ps_v9::detectSpanning_v6(const Index& site) {
//    cout << "Entry -> detectSpanning_v4() : line " << __LINE__ << endl;
    if(_periodicity) {
        cout << "Cannot detect spanning if _periodicity if ON: line " << __LINE__ << endl;
        return false;
    }
    if(_reached_critical ){
        return true;  // we have already reached critical point
    }

    // first check if the site with a cluster id is already a spanning site
    for(const Index& ss: _spanning_sites){
        if(_lattice.getSite(ss).get_groupID() == _lattice.getSite(site).get_groupID()){
//            cout << "Already a spanning site : line " << __LINE__ << endl;
            return true;
        }
    }

    // only check for the newest site placed
    if(site.row_ == min_index){ // top index
        if(!check_if_id_matches(site, _top_edge)) {
            _top_edge.push_back(site);
        }
    }
    else if(site.row_ == max_index){
        if(!check_if_id_matches(site, _bottom_edge)){
            _bottom_edge.push_back(site);
        }
    }

    // checking column indices for Left-Right boundary
    if(site.column_ == min_index){ // left edge
        if(!check_if_id_matches(site, _left_edge)) {
            _left_edge.push_back(site);
        }
    }
    else if(site.column_ == max_index){
        if(!check_if_id_matches(site, _right_edge)) {
            _right_edge.push_back(site);
        }
    }

    if(_number_of_occupied_sites < length()){
//        cout << "Not enough site to span : line " << __LINE__ << endl;
        return false;
    }


    vector<Index>::iterator it_top = _top_edge.begin();
    vector<Index>::iterator it_bot = _bottom_edge.begin();
    bool found_spanning_site = false;
    int id = _lattice.getSite(site).get_groupID();

    if(_top_edge.size() < _bottom_edge.size()){
        // if matched found on the smaller edge look for match in the larger edge
        for(; it_top < _top_edge.end(); ++it_top){
            if(id == _lattice.getSite(*it_top).get_groupID()){
                for(; it_bot < _bottom_edge.end(); ++it_bot){
                    if(id == _lattice.getSite(*it_bot).get_groupID()){
                        // match found !
                        if(!check_if_id_matches(*it_top ,_spanning_sites)) {
                            _reached_critical = true;
                            _spanning_sites.push_back(*it_top);
                        }
                        found_spanning_site = true;
                        _bottom_edge.erase(it_bot);
                    }
                }

                if(found_spanning_site){
                    found_spanning_site = false;
                    _top_edge.erase(it_top);
                }

            }
        }
    }else{
        for (; it_bot < _bottom_edge.end(); ++it_bot) {
            if (id == _lattice.getSite(*it_bot).get_groupID()) {
                for (; it_top < _top_edge.end(); ++it_top) {
                    if (id == _lattice.getSite(*it_top).get_groupID()) {
                        // match found !
                        if (!check_if_id_matches(*it_top, _spanning_sites)) {
                            _reached_critical = true;
                            _spanning_sites.push_back(*it_top);
                        }
                        found_spanning_site = true;
                        _top_edge.erase(it_top);
                    }
                }
                if(found_spanning_site){
                    found_spanning_site = false;
                    _bottom_edge.erase(it_top);
                }
            }
        }

    }

    found_spanning_site = false;
    vector<Index>::iterator it_lft = _left_edge.begin();
    vector<Index>::iterator it_rht = _right_edge.begin();

    if(_left_edge.size() < _right_edge.size()){
        for(; it_lft < _left_edge.end(); ++it_lft) {
            if (id == _lattice.getSite(*it_lft).get_groupID()) {
                for (; it_rht < _right_edge.end(); ++it_rht) {
                    if (id == _lattice.getSite(*it_rht).get_groupID()) {
                        if (!check_if_id_matches(*it_lft, _spanning_sites)) {
                            _spanning_sites.push_back(*it_lft);
                            _reached_critical = true;
                        }
                        found_spanning_site = true;
                        _right_edge.erase(it_rht);
                    }
                }
                if (found_spanning_site) {
                    found_spanning_site = false;
                    _left_edge.erase(it_lft);
                }
            }
        }
    }else{
        for (; it_rht < _right_edge.end(); ++it_rht) {
            if (id == _lattice.getSite(*it_rht).get_groupID()) {
                for(; it_lft < _left_edge.end(); ++it_lft) {
                    if (id == _lattice.getSite(*it_lft).get_groupID()) {
                        if (!check_if_id_matches(*it_lft, _spanning_sites)) {
                            _spanning_sites.push_back(*it_lft);
                            _reached_critical = true;
                        }
                        found_spanning_site = true;
                        _left_edge.erase(it_lft);
                    }
                }
                if (found_spanning_site) {
                    found_spanning_site = false;
                    _right_edge.erase(it_rht);
                }
            }
        }
    }


    // now do the matching with left and right for horizontal spanning
    // meaning new site is added to _spanning_sites so remove them from top and bottom edges



    // filter spanning ids


    return !_spanning_sites.empty();

}



/***********************************
 * Wrapping Detection
 **********************************/
/**
 * Wrapping is detected here using the last placed site
 * @return bool. True if wrapping occured.
 */
bool SitePercolation_ps_v9::detectWrapping() {
    Index site = lastPlacedSite();
    // only possible if the cluster containing 'site' has sites >= length of the lattice
    if(_number_of_occupied_sites < length()){
        return false;
    }

    if(_reached_critical){
        return true; // reached critical in previous step
    }
    // check if it is already a wrapping site
    int id = _lattice.getSite(site).get_groupID();
    int tmp_id{};
    for (auto i: _wrapping_sites){
        tmp_id = _lattice.getSite(i).get_groupID();
        if(id == tmp_id ){
            return true;
        }
    }

    // get four neighbors of site always. since wrapping is valid if periodicity is implied
    vector<Index> sites = _lattice.get_neighbor_site_indices(site);

    if(sites.size() < 2){ // at least two neighbor of  site is required
        return false;
    }else{
        IndexRelative irel = _lattice.getSite(site).relativeIndex();
//        cout << "pivot's " << site << " relative " << irel << endl;
        IndexRelative b;
        for (auto a:sites){
            if(_lattice.getSite(a).get_groupID() != _lattice.getSite(site).get_groupID()){
                // different cluster
                continue;
            }
            // belongs to the same cluster
            b = _lattice.getSite(a).relativeIndex();
//            cout << "neibhbor " << a << " relative " << b << endl;
            if(abs(irel.x_ - b.x_) > 1 || abs(irel.y_ - b.y_) > 1){
//                cout << "Wrapping : line " << __LINE__ << endl;
                _wrapping_sites.push_back(site);
                _reached_critical = true;
                return true;
            }
        }
    }
    // if %_wrapping_indices is not empty but wrapping is not detected for the current site (%site)
    // that means there is wrapping but not for the %site
    return !_wrapping_sites.empty();
}

/********************************************************************
 * Relabeling
 *
 *********************************************************************/
/**
 * Relabels site and also reassign relative index to the relabeled sites
  *
  * @param site_a  : last added site index of the base cluster
  * @param clstr_b : 2nd cluster, which to be merged withe the root
  */
void SitePercolation_ps_v9::relabel_sites_v5(Index site_a, const Cluster& clstr_b) {
    const vector<Index> sites = clstr_b.getSiteIndices();
    int id_a = _lattice.getSite(site_a).get_groupID();
    int id_b = clstr_b.get_ID();
    Index b = clstr_b.getRootSite();

    // get four site_b of site_a
    vector<Index> sites_neighbor_a = _lattice.get_neighbor_site_indices(site_a);
    Index site_b;
    IndexRelative relative_index_b_after;
    bool flag{false};
    // find which site_b has id_a of clstr_b
    for(auto n: sites_neighbor_a){
        if(id_b == _lattice.getSite(n).get_groupID()){
            // checking id_a equality is enough. since id_a is the id_a of the active site already.
            relative_index_b_after = getRelativeIndex(site_a, n);
            site_b = n;
            flag = true;
            break;
        }
    }
    if(!flag){
        cout << "No neibhgor found! : line " << __LINE__ << endl;
    }

    IndexRelative relative_site_a = _lattice.getSite(site_a).relativeIndex();
    // with this delta_a and delta_y find the relative index of site_b while relative index of site_a is known
    IndexRelative relative_site_b_before = _lattice.getSite(site_b).relativeIndex();
    int delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    int delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;
    relabel_sites(sites, id_a, delta_x_ab, delta_y_ab);
}



void SitePercolation_ps_v9::relabel_sites(const vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab)  {
    int x, y;
    Index a;
    IndexRelative relative_site__a;
    for (value_type i = 0; i < sites.size(); ++i) {
        a = sites[i];
        _lattice.getSite(a).set_groupID(id_a);
        relative_site__a = _lattice.getSite(a).relativeIndex();
        x = relative_site__a.x_ + delta_x_ab;
        y = relative_site__a.y_ + delta_y_ab;
        _lattice.getSite(a).relativeIndex(x, y);
    }
}



/**********************************************
 * Information about current state of Class
 **********************************************/

/**
 * Entropy calculation is performed here. The fastest method possible.
 * Cluster size is measured by bond.
 * @return current entropy of the lattice
 */
double SitePercolation_ps_v9::entropy() {
    double H{};
    double number_of_cluster_with_size_one = maxBonds() - _bonds_in_cluster_with_size_two_or_more;
//    cout << " _bonds_in_cluster_with_size_two_or_more " << _bonds_in_cluster_with_size_two_or_more << " : line " << __LINE__ << endl;
    double mu = 1.0/double(maxBonds());
    H += number_of_cluster_with_size_one * log(mu) * mu;
    H *= -1;
    _entropy_current =  _entropy + H;
    return _entropy_current;
}



/**
 * Only applicable if the number of bonds in the largest cluster is calculated when occupying the lattice.
 * Significantly efficient than the previous version numberOfBondsInTheLargestCluster()
 * @return
 */
value_type SitePercolation_ps_v9::numberOfBondsInTheLargestCluster_v2() {
//    return _clusters[_index_largest_cluster].numberOfBonds();
    return _number_of_bonds_in_the_largest_cluster;
}



/**
 *
 * @return
 */
value_type SitePercolation_ps_v9::numberOfSitesInTheLargestCluster() {
    value_type  len{}, nob{};
    for(auto c: _clusters){
        nob = c.numberOfSites();
        if (len < nob){
            len = nob;
        }
    }
    _number_of_sites_in_the_largest_cluster = len;
    return len;
}


/**********************************
 * Spanning methods
 **********************************/

/**
 *
 * @return
 */
value_type SitePercolation_ps_v9::numberOfSitesInTheSpanningClusters_v2() {

    if(! _spanning_sites.empty()){
        int id = _lattice.getSite(_spanning_sites.front()).get_groupID();
        if(id >= 0) {
            return _clusters[id].numberOfSites();
        }
    }
    return 0;
}


/**
 *
 * @return
 */
value_type SitePercolation_ps_v9::numberOfBondsInTheSpanningClusters_v2() {
    if(!_spanning_sites.empty()){
//        cout << "number of spanning sites " << _spanning_sites.size() << " : line " << __LINE__ << endl;
        int id = _lattice.getSite(_spanning_sites.front()).get_groupID();
        if(id >= 0) {
            return _clusters[id].numberOfBonds();
        }
    }
    return 0;
}

/**
 *
 * @return
 */
value_type SitePercolation_ps_v9::numberOfSitesInTheWrappingClusters(){
    value_type nos{};
    int id{};
    for(auto i: _wrapping_sites){
        id = _lattice.getSite(i).get_groupID();
        if(id >= 0) {
            nos += _clusters[id].numberOfSites();
        }
    }
    return nos;
}

/**
 *
 * @return
 */
value_type SitePercolation_ps_v9::numberOfBondsInTheWrappingClusters(){
    value_type nob{};
    int id{};
    for(auto i: _wrapping_sites){
        id = _lattice.getSite(i).get_groupID();
        if(id >= 0) {
            nob += _clusters[id].numberOfBonds();
        }
    }
    return nob;
}


std::string SitePercolation_ps_v9::getSignature() {
    string s = "sq_lattice_site_percolation";
    if(_periodicity)
        s += "_periodic_";
    else
        s += "_non_periodic_";
    return s;
}

/**
 *
 * @param filename
 * @param only_spanning
 */
void SitePercolation_ps_v9::writeVisualLatticeData(const string &filename, bool only_spanning) {
    std::ofstream fout(filename);
    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length()
                << ",\"signature\":\"" << getSignature() << "\""
                << ",\"x\":\"" << lastPlacedSite().column_ << "\""
                << ",\"y\":\"" << lastPlacedSite().row_ << "\""
                << "}" ;

    fout << "#" << header_info.str() << endl;
    fout << "#<x>,<y>,<color>" << endl;
    fout << "# color=0 -means-> unoccupied site" << endl;
    int id{-1};
    if(!_spanning_sites.empty()){
        id = _lattice.getSite(_spanning_sites.front()).get_groupID();
    }
    else if(!_wrapping_sites.empty()){
        id = _lattice.getSite(_wrapping_sites.front()).get_groupID();
    }

    if(only_spanning){
        if(id < 0){
            cerr << "id < 0 : line " << __LINE__ << endl;
        }
        vector<Index> sites = _clusters[id].getSiteIndices();
        for(auto s: sites){
            fout << s.column_ << ',' << s.row_ << ',' << id << endl;
        }
    }
    else {
        for (value_type y{}; y != length(); ++y) {
            for (value_type x{}; x != length(); ++x) {
                id = _lattice.getSite({y, x}).get_groupID();
                if(id != -1) {
                    fout << x << ',' << y << ',' << id << endl;
                }
            }
        }
    }
    fout.close();
}


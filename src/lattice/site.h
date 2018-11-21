//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_SITE_H
#define SITEPERCOLATION_SITE_H

#include <array>
#include <set>
#include <vector>
#include <iostream>
#include <memory>

#include "../index/index.h"
#include "../types.h"



/**
 * single Site of a lattice
 */
struct Site{
    /**
     * if true -> site is placed.
     * if false -> the (empty) position is there but the site is not (required for site percolation)
     */
    bool _status{false};
    int _group_id{-1};
    Index _id{};

    //relative distance from the root site. {0,0} if it is the root site
    //for detecting wrapping
    IndexRelative _relative_index{0,0};


public:

    ~Site()                 = default;
    Site()                  = default;
    Site(const Site&)             = default;
    Site(Site&&)            = default;
    Site& operator=(const Site&)  = default;
    Site& operator=(Site&&) = default;

    Site(Index id, value_type length){
        // I have handle _neighbor or corner points and edge points carefully
        if(id.row_ >= length || id.column_ >= length){
            std::cout << "out of range : line " << __LINE__ << std::endl;
        }
        _id.row_ = id.row_;
        _id.column_ = id.column_;
    }


    bool isActive() const { return _status;}
    void activate(){ _status = true;}
    void deactivate() {
        _relative_index = {0,0};
        _group_id = -1;
        _status = false;
    }
    Index ID() const { return  _id;}
    /*
     * Group get_ID is the set_ID of the cluster they are in
     */
    int     get_groupID() const {return _group_id;}
    void    set_groupID(int g_id) {_group_id = g_id;}

    std::stringstream getSite() const {
        std::stringstream ss;
        if(isActive())
            ss << _id;
        else
            ss << "(*)";
        return ss;
    }


    void relativeIndex(IndexRelative r){
        _relative_index = r;
    }

    void relativeIndex(int x, int y){
        _relative_index = {x,y};
    }

    IndexRelative relativeIndex() const {return _relative_index;}
};

std::ostream& operator<<(std::ostream& os, const Site& site);
bool operator==(Site& site1, Site& site2);


#endif //SITEPERCOLATION_SITE_H

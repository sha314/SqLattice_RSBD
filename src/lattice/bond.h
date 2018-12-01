//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_BOND_H
#define SITEPERCOLATION_BOND_H


#include <ostream>
#include <iostream>
#include <sstream>

#include "../index/index.h"
#include "../types.h"



/**
 * A bond has two end
 * say a 5x5 lattice bond between end1 (0,0) and end2 (0,1)
 * if _status is false -> bond is not there
 *
 */
struct Bond{
    // check if active or not
    bool _status{false};
    value_type _length;
    int _group_id{-1};

    BondType bondType;
    //relative distance from the root site. {0,0} if it is the root site
    IndexRelative _relative_index{0,0};
    Index _end1;
    Index _end2;
    BondIndex _id;

    ~Bond() = default;
    Bond() = default;
    Bond(Index end1, Index end2, value_type length){
        _end1.row_ = end1.row_;
        _end1.column_ = end1.column_;
        _end2.row_ = end2.row_;
        _end2.column_ = end2.column_;

        _length = length;
        // check if the bond is valid
        if(_end1.row_ == _end2.row_){
            bondType = BondType::Horizontal;
            // means x_ values are equal
            if(_end1.column_ > _end2.column_){
                if(_end1.column_ == _length-1 && _end2.column_ ==0){
                    // do nothing
                }
                else{
                    // sort them out
                    _end1.column_ = end2.column_;
                    _end2.column_ = end1.column_;
                }
            }
            else if(_end1.column_ < _end2.column_){
                if(_end1.column_ == 0 && _end2.column_ == _length-1){
                    _end1.column_ = end2.column_;
                    _end2.column_ = end1.column_;
                }
            }
        }
        else if(_end1.column_ == _end2.column_){
            bondType = BondType::Vertical;
            // means y_ values are equal
            if(_end1.row_ > _end2.row_){
                if(_end1.row_ == _length-1 && _end2.row_ ==0){
                }
                else{
                    // sort them out
                    _end1.row_ = end2.row_;
                    _end2.row_ = end1.row_;
                }
            }
            else if(_end1.row_ < _end2.row_){
                if(_end1.row_ == 0 && _end2.row_ == _length-1){
                    _end1.row_ = end2.row_;
                    _end2.row_ = end1.row_;
                }
            }
        }
        else{
            std::cout << '(' << _end1.row_ << ',' << _end1.column_ << ')' << "<->"
                    << '(' << _end2.row_ << ',' << _end2.column_ << ')'
                    << " is not a valid bond : line " << __LINE__ << std::endl;
        }

        _id = BondIndex(bondType, _end1.row_, _end1.column_);  // unsigned long
    }

    std::vector<Index> getSites() const { return {_end1, _end2};}

    Index id() const {
        return _end1;
    }

    BondIndex ID() const {
        return _id;
    }

    void activate() {_status = true;}
    void deactivate() {
        _relative_index = {0,0};
        _group_id = -1;
        _status = false;
    }
    bool isActive() const { return _status;}
/*
* Group get_ID is the set_ID of the cluster they are in
*/
    int get_groupID() const {return _group_id;}
    void set_groupID(int g_id) {_group_id = g_id;}

    std::stringstream getBondString() const {
        std::stringstream ss;
        if(isActive()) {
            // place '-' for horizontal bond and '|' for vertical bong
            if(bondType == BondType::Horizontal) {
                ss << '(' << _end1 << "<->" << _end2 << ')';
            }
            else {
                ss << '(' << _end1 << "<|>" << _end2 << ')';
            }
        }
        else
            ss << "(**)";
        return ss;
    }

    bool isHorizontal() const { return bondType == BondType ::Horizontal;}
    bool isVertical()   const { return bondType == BondType ::Vertical;}

    void relativeIndex(IndexRelative r){
        _relative_index = r;
    }

    void relativeIndex(int x, int y){
        _relative_index = {x,y};
    }

    IndexRelative relativeIndex() const {return _relative_index;}
};



std::ostream&   operator<<(std::ostream& os, const Bond& bond);
bool            operator==(Bond a, Bond b);
bool            operator<(const Bond& bond1, const Bond& bond2);
bool            operator>(const Bond& bond1, const Bond& bond2);


#endif //SITEPERCOLATION_BOND_H

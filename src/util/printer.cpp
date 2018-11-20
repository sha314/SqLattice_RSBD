//
// Created by shahnoor on 8/16/18.
//

#include "printer.h"
#include <iostream>

using namespace std;


void print_h_barrier(size_t n, const string& initial, const string& middles, const string& end){
    cout << initial;
    for(size_t i{}; i < n ; ++i){
        cout << middles;
    }
    cout << end; // end of barrier
}

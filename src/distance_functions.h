#pragma once
#define __STDC_FORMAT_MACROS
template <class T>
bool is_in_neighborhood8(Position * p1, Position * p2){
    if(labs(p1->x - p2->x) > 1) return false;
    if(labs(p1->y - p2->y) > 1) return false;
    if(labs(p1->z - p2->z) > 1) return false;
    return true;
}

bool all_together(Position * p1, Position * p2){
    return true;
}
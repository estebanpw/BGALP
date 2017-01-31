#pragma once

#define READBUF 1000
#define POOL_SIZE 1024*1024*1024 //1 MB 
#define LIST_POOL_SIZE 128

template <typename T>
struct Unit{
    T elem;
    struct Unit * next;
};

struct Problem{
    int x;
};
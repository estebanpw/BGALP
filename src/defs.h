#pragma once

#define READBUF 1000
#define POOL_SIZE 1024*1024*1024 //1 MB 
#define LIST_POOL_SIZE 128

template <typename T>
struct Unit{
    T elem;
    struct Unit * next;
};

struct Position{
    uint64_t x;
    uint64_t y;
    uint64_t z;
    Position(uint64_t x = 0, uint64_t y = 0, uint64_t z = 0): x(x), y(y), z(z) {}
};
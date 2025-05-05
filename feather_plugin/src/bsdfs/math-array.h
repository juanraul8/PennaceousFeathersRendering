#pragma once
#ifndef _MATH_ARRAY_H_
#define _MATH_ARRAY_H_
#endif

#include <array>
#include <list>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <iostream>

template<typename NumA, typename NumB, std::size_t N, typename F>
auto zipWith(const F& f, const std::array<NumA,N>& a, const std::array<NumB,N>& b) {
    std::array<std::decay_t<decltype(f(a[0],b[0]))>,N> r;
    std::transform(a.begin(),a.end(),b.begin(),r.begin(),f);
    return r;
}

template<typename NumA, std::size_t N, typename F>
auto map(const F& f, const std::array<NumA,N>& a) {
    std::array<std::decay_t<decltype(f(a[0]))>,N> r;
    std::transform(a.begin(),a.end(),r.begin(),f);
    return r;
}

template<typename NumA, typename NumIni, std::size_t N, typename F>
auto foldl(const F& f, const NumIni& ini, const std::array<NumA,N>& a) {
    return std::accumulate(a.begin(),a.end(), ini, f);
}

template<typename NumA, typename NumB, std::size_t N>
auto operator+(const std::array<NumA,N>& a, const std::array<NumB,N>& b) {
    return zipWith([] (const NumA& x, const NumB& y) { return x+y; },a,b);
}

template<typename NumA, typename NumB, std::size_t N>
auto operator-(const std::array<NumA,N>& a, const std::array<NumB,N>& b) {
    return zipWith([] (const NumA& x, const NumB& y) { return x-y; },a,b);
}

template<typename NumA, typename NumB, std::size_t N>
auto operator*(const std::array<NumA,N>& a, const std::array<NumB,N>& b) {
    return zipWith([] (const NumA& x, const NumB& y) { return x*y; },a,b);
}

template<typename NumA, typename NumB, std::size_t N>
auto operator/(const std::array<NumA,N>& a, const std::array<NumB,N>& b) {
    return zipWith([] (const NumA& x, const NumB& y) { return x/y; },a,b);
}

template<typename NumA, typename NumP,std::size_t N>
auto operator*(const std::array<NumA,N>& a, const NumP& f) {
    return map([f] (const NumA& x) { return f*x; }, a);
}

template<typename NumP, typename NumA,std::size_t N>
auto operator*(const NumP& f, const std::array<NumA,N>& a) {
    return a*f;
}

template<typename NumA, typename NumP,std::size_t N>
auto operator/(const std::array<NumA,N>& a, const NumP& f) {
    return map([f] (const NumA& x) { return x/f; }, a);
}

template<typename NumA, typename NumP,std::size_t N>
auto operator/(const NumP& f, const std::array<NumA,N>& a) {
    return map([f] (const NumA& x) { return f/x; }, a);
}


template<typename Num, std::size_t N>
std::array<Num,N> operator-(const std::array<Num,N>& a) {
    return a*Num(-1);
}

template<typename T>
auto sqr(const T& v) {
    return v*v;
}

template<typename Num, std::size_t N>
Num sum(const std::array<Num,N>& a) {
    return foldl([] (const Num& x, const Num& y) { return x+y; }, Num(0), a);
}

template<typename NumA, typename NumB, std::size_t N>
auto dot(const std::array<NumA,N>& a, const std::array<NumB,N>& b) {
    return sum(a*b);
}

template<typename NumA, typename NumB>
auto cross(const std::array<NumA,3>& a, const std::array<NumB,3>& b) {
    return std::array<NumA,3>{
        a[1]*b[2]-a[2]*b[1],
        a[2]*b[0]-a[0]*b[2],
        a[0]*b[1]-a[1]*b[0]
    };
}

template<typename Num, std::size_t N>
auto length(const std::array<Num,N>& a) {
    return std::sqrt(dot(a,a));
}

template<typename Num, std::size_t N>
auto normalized(const std::array<Num,N>& a) {
    return a/length(a);
}

template<typename Num, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<Num,N>& a) {
    os<<"[ ";
    for (std::size_t i = 0; i<N; ++i) os<<a[i]<<" ";
    os<<"]";
    return os;
}


template<typename Num>
Num determinant(const std::array<std::array<Num,2>,2>& m) {
    return m[0][0]*m[1][1] - m[0][1]*m[1][0];    
}

template<typename Num>
Num determinant(const std::array<std::array<Num,3>,3>& m) {
    return m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1]
        - m[0][2]*m[1][1]*m[2][0] - m[0][1]*m[1][0]*m[2][2] - m[0][0]*m[1][2]*m[2][1];   
}

template<typename Num>
Num adjoint(const std::array<std::array<Num,3>,3>& m, std::size_t i, std::size_t j) {
    return m[(j+1)%3][(i+1)%3]*m[(j+2)%3][(i+2)%3] - m[(j+1)%3][(i+2)%3]*m[(j+2)%3][(i+1)%3];
}

template<typename Num>
std::array<std::array<Num,3>,3> inverse(const std::array<std::array<Num,3>,3>& m) {
    Num f = 1.0f/determinant(m);
    return std::array<std::array<Num,3>,3>{
        std::array<Num,3>{f*adjoint(m,0,0),f*adjoint(m,0,1),f*adjoint(m,0,2)},
        std::array<Num,3>{f*adjoint(m,1,0),f*adjoint(m,1,1),f*adjoint(m,1,2)},
        std::array<Num,3>{f*adjoint(m,2,0),f*adjoint(m,2,1),f*adjoint(m,2,2)}
    };
}

template<typename Num, std::size_t N> 
std::array<std::array<Num,3>,N> matrix(const std::array<Num,N>& c0, const std::array<Num,N>& c1, const std::array<Num,N>& c2) {
    std::array<std::array<Num,3>,N> m;
    for (std::size_t i = 0; i<N; ++i) {
        m[0][i] = c0[i];
        m[1][i] = c1[i];
        m[2][i] = c2[i];
    } 
    return m;   
}

template<typename NumA, typename NumB, std::size_t AH, std::size_t N>
auto product(const std::array<std::array<NumA,N>,AH>& m, const std::array<NumB,N>& v) {
    return map([&v] (const std::array<NumA,N>& row) { return dot(row,v); }, m);
}

template<typename Num, std::size_t AH, std::size_t N, std::size_t BW>
auto matrix_product(const std::array<std::array<Num,N>,AH>& m1, const std::array<std::array<Num,BW>,N>& m2) {
    std::array<std::array<Num,BW>,AH> s;
    for (std::size_t r = 0; r<AH; ++r)
        for (std::size_t c = 0; c<BW; ++c) {
            s[r][c] = Num(0);
            for (std::size_t i = 0; i<N;++i)
                s[r][c] += m1[r][i]*m2[i][c];
        }
    return s;
}





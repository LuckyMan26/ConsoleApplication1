// Halide tutorial lesson 1: Getting started with Funcs, Vars, and Exprs

// This lesson demonstrates basic usage of Halide as a JIT compiler for imaging.

// On linux, you can compile and run it like so:
// g++ lesson_01*.cpp -g -I <path/to/Halide.h> -L <path/to/libHalide.so> -lHalide -lpthread -ldl -o lesson_01 -std=c++17
// LD_LIBRARY_PATH=<path/to/libHalide.so> ./lesson_01

// On os x:
// g++ lesson_01*.cpp -g -I <path/to/Halide.h> -L <path/to/libHalide.so> -lHalide -o lesson_01 -std=c++17
// DYLD_LIBRARY_PATH=<path/to/libHalide.dylib> ./lesson_01

// If you have the entire Halide source tree, you can also build it by
// running:
//    make tutorial_lesson_01_basics
// in a shell with the current directory at the top of the halide
// source tree.

// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"

// We'll also include stdio for printf.
#include <stdio.h>

Halide::Func matrix_power_p(Halide::Func& input, int exponent, int cols) {
    Halide::Var x, y,d;
    Halide::Func result;
    if (exponent == 0) {
        result(x, y,d) = select(x == y, 1, 0);
    }
    else if (exponent > 0) {
        int e = exponent - 1;
        result(x, y,d) = input((x - e) % cols, y,d);

    }

    else {
        int e = (exponent % cols + cols) - 1;
        result(x, y, d) = input((x + e) % cols, y,d);
    }


    return result;
}
Halide::Func matrix_multiplication(Halide::Func& L, Halide::Func& P) {
    Halide::Func res;
    Halide::RDom r(0, 4, "r");
    Halide::Var x, y,d;
    res(x, y,d) = Halide::sum(L(r, y,d) * P(x, r,d));
    std::cout << std::endl;
    Halide::Buffer<int> b = res.realize({ 4,4,4 });
 
    return res;
}
Halide::Func matrix_power_q(Halide::Func& input, int exponent, int rows) {
    Halide::Var x, y,d;
    Halide::Func result;
    if (exponent == 0) {
        result(x, y,d) = select(
        (x==y),1,
            0
        );
    }
    else if (exponent > 0) {
        int e = exponent - 1;
        result(x, y,d) = input(x, (y + e) % rows,d);

    }

    else {
        int e = (exponent % rows + rows) - 1;
        result(x, y,d) = input(x, (y - e) % rows,d);

    }

  
  
    return result;
}
int main(int argc, char** argv) {
  
    Halide::Var s, t, x, y, d;
    Halide::Func  X, B, M;

  
    int m, n;
    m = 4;
    n = 4;
    int Mi, N;
    Mi = 4;
    N = 4;
   
 
    Halide::Func P;
    Halide::Func Q;

    P(x, y,d) = select(
        (x == 0 && y < m - 1), 0,
        (x == 0 && y == m - 1), 1,
        (y == m - 1 && x > 1), 0,
        (x > 0 && y < m - 1 && y == x - 1), 1,

        0
    );
    Q(x, y,d) = select(
        (y == 0 && x < n - 1), 0,
        (y == 0 && x == n - 1), 1,
        (y > 0 && x < n - 1 && y - 1 == x), 1,

        0
    );

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Halide::Func power_p = matrix_power_p(P, -i, Mi);
            Halide::Func power_q = matrix_power_q(Q, -j, N);
            B(x, y, s, t, d) = matrix_multiplication(power_p, power_q)(x, y, d);
        }
    }
    M(x, y, s, t, d) = 0;
   
    M(x, y, s, t, d) = select(
        (s == 0 && t == 0), B(x, y, 0, 0, d),
        (s > 0 && t == 0), B(x, y, s, 0, d) + M(x, y, s - 1, 0, d),
        (s == 0 && t > 0),B(x, y, 0, t, d) + M(x, y, 0, t-1, d),
        (s > 0 && t > 0),  M(x, y, s-1, t , d) + M(x, y, s, t-1, d)-M(x,y,s-1,t-1,d)+B(x,y,s,t,d),
        0);
    
   

 
    Halide::RDom r(0, m, 0, n);
    Halide::Func L, G, K, J;
    L(x, y, d) = M(x,y,m-1,n-1,d)-M(x,y,m-1,n-t-1,d)-M(x,y,m-s-1,n-1,d) + M(x,y,m-s-1,n-t-1,d); 
    G(x, y, d) = M(x,y,m-1,n-y-1,d) - M(x,y,m-s-1,n-1,d); 
    K(x, y, d) = M(x, y, m - s- 1, n - 1, d) - M(x, y, m - s - 1, n -t- 1, d);
    J(x, y, d) = M(x, y, m - s - 1, n -t- 1, d);

  
    Halide::Func Z;
    Z(s, t,d) = L(x,y,d) + G(x, y, d) + K(x, y, d) + J(x, y, d);
    Halide::Func res = Z.vectorize(d);
    return 0;
}

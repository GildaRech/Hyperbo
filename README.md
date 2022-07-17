

# Hyperbo.    

**Hyperbo** is an open-source python package, a mathematical package that implements arithmetical computations on hyperbola. This work in particular, focuses on research findings on lattice points on hyperbolas in relation with Fermat factorization equation in [Lattice Points on Fermat Factorization Method](https://www.hindawi.com/journals/jmath/2022/6360264/).  .

It implements arithmetical results in relation with the Fermat factorization equation. This includes lattice points and solution computations, cardinals, group operations like addition and scalar multiplications on different algebraic structures, plots and general algebraic information and homomorphism relations.   

> This package considers mainly the following hyperbola forms to which any other parametrization can result through a morphism:    
            
![](https://github.com/GildaRech/Hyperbo/blob/main/img2.png?raw=true) 
---
The Library has 3 classes: **Common**, **H** and **B**.

- [x] **Common**:   

This class contains common methods used in hyperbola parametrizations mainly $H_n$ and $B_n$.  

    It provides its own verified autonomous functions for primality testing, gcd, inverse modulo, factor, and prime factorization.
    This is for the purpose of reducing dependencies to external libraries.   
    
    FUNCTIONS:   
    
            ```rabinMiller(x)```: function that returns True or False whether the integer x is prime or not using RabinMiller primality test algorithm.  
                            e.g: rabinMiller(19) returns True.            
            ```is_prime(x)```: function that returns True or False if the given integer x is prime or not. This function is optimized and therefore has a 
                            better complexity. 
                         e.g: is_prime(20) returns False.
            ```pgcd(x1, x2)```: function that returns the greatest common divisor of x1 and x2.
                          eg: pgcd(10, 7) returns 1
            ```extended_euclidean(x1, x2)```: function that returns the gcd and bezout coefficients x, y such that x1*x+x2*y = gcd(x1, x2)
                                        e.g: extended_euclidean(11, 13) returns (1, 6, -5): (gcd, coeff_x, coeff_y)
            ```inverse_modulo(x n)```: function that returns the inverse of x modulo n.
                                 eg: inverse_modulo(5, 11) returns -2
            ```facto(x)```: function that returns the list of factors of the integer x with their degrees of multiplicity.                     
                      eg: facto(60) returns [(2, 2), (3, 1), (5, 1)]
            ```pfactors(x)```: function that returns a list of prime factors of x.
                         e.g: pfactors(15) returns [3, 5]
            ```is_square(x)```: function that returns True or False whether the integer x is prime or not.
                         e.g: is_square(16) returns True
            ```is_same(list)```: function that returns True or False whether the elements of the list are the same.
                          e.g: is_same([2, 2, 2]) returns True and is_same([2, 3, 4, 3]) retunrs False
            ```pair_sort(list)```: function that returns the sorted list of tuples by first element.
                             eg: pair_sort([(4, 6), (2, 5), (5, 7)]) returns [(2, 5), (4, 6), (5, 7)]
---
- [x] **H**:   
 
 This class implements methods used in hyperbola parametrizations $H_n$.    
 
            It provides methods related to the object $H_n$.   
            
            FUNCTIONS:   
            
                     ```is_fermat_solvable```: property that returns True or False whether the Fermat equation has a solution or not.
                                         e.g: for n=15 returns True.  
                     ```info```: property that returns the general info on H_n
                     ```is_in_H(x)```: function that returns True or False whether if the point x is in H or not.
                     ```negativPoints(P)```: function that returns negative points through Symmetry to P on Hn(x, y).
                     ```card```: property that returns the cardinal of H_n(x, y).
                     ```points```: property that returns the list of points of H_n over self.S structure
                     ```add(P,Q)```: function that adds two points P and Q on H_n(x, y).
                     ```double(P)```: function that doubles a point P on H_n(x, y).
                     ```mul(k, P)```: function that multiplies a point P by a scalar k on H_n(x, y).
                     ```plot```: property that plots points on H_n(x, y). 
  ---
- [x] **B**:   

This class implements methods used in hyperbola parametrizations $B_n$.    

            It provides methods related to the object $B_n$.    
            
            FUNCTIONS:    
            
                     ```info```: property that prints the general info about the object B_n(x, y) over self.S structure.
                     ```is_in_B(P)```: function that checks whether a point P is in B_n(x, y).
                     ```nbr_pointsS4```: property that returns the number of points on B_n over Z4.
                     ```_points```: property that returns points on B_n over Z4.
                     ```U(i)```: function that returns the i term of the sequence U(i). i represents the number of primes.
                           E.g: U(2) returns 5
                    ```card```: property that returns the cardinal of B_n(x, y).
                    ```add(P, Q)```: function that adds two points P and Q on B_n(x, y).
                    ```double(P)```: function that doubles a point on B_n(x, y).
                    ```mul(k, P)```: function that multiplies a point P by a scalar k on B_n(x, y).
                    ```card_sum```: property that returns the sum S_n of cardinals on B_n(x, y).
                    ```_productp```: property that returns the product of prime divisors of n that make up n on B_n(x, y).
                    ```pointsZ4```: property that returns points on B_n over Z4 using algebraic results on B_n.
                    ```negativPoints(l)```: function that returns the negative points on B_n(x, y) by symmetry from points in the list l. 
                    ```points```: function that returns points on B_n over different algebraic structures.
                    ```plot```: property that plots points on B_n(x, y).

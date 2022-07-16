

# Hyperbo
**Hyperbo** is an open-source python package, a mathematical package that implements arithmetical computations on hyperbola. This work in particular, focuses on research findings on lattice points on hyperbolas in relation with Fermat factorization equation.

It implements arithmetical results in relation with the Fermat factorization equation. This includes lattice points and solution computations, cardinals, group operations like addition and scalar multiplications on different algebraic structures, plots and general algebraic information and homomorphism relations.

            
![](https://github.com/GildaRech/Hyperbo/blob/main/img2.png?raw=true) 
---
The Library has 3 classes: **Common**, **H** and **B**.

- [x] **Common**:
.
This class contains common methods used in hyperbola parametrizations mainly H_n and B_n.  

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
- [x] **H**
           

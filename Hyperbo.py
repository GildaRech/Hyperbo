#! -*- coding: utf-8 -*-
import sys, random
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import pylab



class Common:
    """ 
    This class contains common methods used in hyperbola parametrizations mainly H_n and B_n.
    It provides its own verified autonomous functions for primality testing, gcd, inverse modulo, factor, and prime factorization.
    This is for the purpose of reducing dependencies to external libraries.
    FUNCTIONS:
            rabinMiller(x): function that returns True or False whether the integer x is prime or not using RabinMiller primality test algorithm.
                            e.g: rabinMiller(19) returns True.    
            is_prime(x): function that returns True or False if the given integer x is prime or not. This function is optimized and therefore has a better complexity. 
                         e.g: is_prime(20) returns False.
            pgcd(x1, x2): function that returns the greatest common divisor of x1 and x2.
                          eg: pgcd(10, 7) returns 1
            extended_euclidean(x1, x2): function that returns the gcd and bezout coefficients x, y such that x1*x+x2*y = gcd(x1, x2)
                                        e.g: extended_euclidean(11, 13) returns (1, 6, -5): (gcd, coeff_x, coeff_y)
            inverse_modulo(x n): function that returns the inverse of x modulo n.
                                 eg: inverse_modulo(5, 11) returns -2
            facto(x): function that returns the list of factors of the integer x with their degrees of multiplicity.                     
                      eg: facto(60) returns [(2, 2), (3, 1), (5, 1)]
            pfactors(x): function that returns a list of prime factors of x.
                         e.g: pfactors(15) returns [3, 5]
            is_square(x): function that returns True or False whether the integer x is prime or not.
                         e.g: is_square(16) returns True
            is_same(list): function that returns True or False whether the elements of the list are the same.
                          e.g: is_same([2, 2, 2]) returns True and is_same([2, 3, 4, 3]) retunrs False
            pair_sort(list): function that returns the sorted list of tuples by first element.
                             eg: pair_sort([(4, 6), (2, 5), (5, 7)]) returns [(2, 5), (4, 6), (5, 7)]      
    """
    def __init__(self, n) -> None:
        self.n=n
        pass
    
    def rabinMiller(self, nbre:int) -> bool:
        """ function that implements the rabinMiller primality test.
         Args:
            n (int): integers that one wants to check the primality
         Returns:
            bool: True or False whether the integer n is prime or not  
         Eg: rabinMiller(100) returns False.
        """
        self.nbre=nbre
        "returns True or False whether prime or not"
        s = self.nbre - 1
        t = 0
        while s % 2 == 0:
                s = s // 2
                t += 1
        for essai in range(5): 
            a = random.randrange(2, self.nbre - 1)
            v = pow(a, s, self.nbre)
            if v != 1:
                i = 0
                while v != (self.nbre - 1):
                    if i == t - 1:
                        return False
                    else:
                        i = i + 1
                        v = (v ** 2) % self.nbre
        return True   
    
    def is_prime(self, nbre:int) -> bool:
        """ function that checks the primality of an integer.
            It particularly uses the rabinMiller and the below list of first primes less than 1000.
         Args:
            n (int): integers that one wants to check the primality
         Returns:
            bool: True or False whether the integer n is prime or not  
         Eg: is_prime(31) returns True.
        """
        self.nbre=nbre
        "faster for performance purpose"
        if (self.nbre < 2):
            return False
        self.primesN = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
                    109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
                    233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
                    367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
                    499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641,
                    643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787,
                    797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
                    947, 953, 967, 971, 977, 983, 991, 997]
        if self.nbre in self.primesN:
            return True
        for prime in self.primesN:
            if (self.nbre % prime == 0):
                    return False
        return self.rabinMiller(self.nbre)

    def pgcd(self, a, b)->int:
        """ function that returns the greatest common divisor (gcd) of two given integer inputs
         Args:
            a, b (int): integers that one wants to compute the gcd
         Returns:
            int: the gcd of a and b 
         Eg: pgcd(5, 10) returns 5.
        """
        if not a>b: a, b =b, a
        while b != 0:
            r= a%b; a=b; b=r
        return a

    def extended_euclidean(self, a, b) -> tuple:
        """ function that returns the gcd and bezout coefficients x and y such that ax+by=pgcd(a, b)=d.
         Args:
            a, b (int): integers that one wants to compute the bezout coefficients
         Returns:
            tuple: (d, x, y), d the gcd of a and b; x and y the coefficients respectively of a and b. 
         Eg: extended_euclidean(11, 13) returns (1, 6, -5).
        """
        if b==0:
            d=a; x=1; y=0
            return (d, x, y)
        x2=1; x1=0; y2=0; y1=1
        while b>0:
            q=int(a/b); r=a-q*b; x=x2-q*x1; y=y2-q*y1
            a=b; b=r; x2=x1; x1=x; y2=y1; y1=y
        d=a; x=x2; y=y2
        return (d, x, y)

    def inverse_modulo(self, a, n) -> int:
        """ function that returns the inverse of a given integer in a finite field of characteristic n ie the inverse of the integer a modulo n.
         Args:
            a, n (int): a is the integer that one wants to compute the inverse and n is the integer modulo
         Returns:
            int: the inverse of a modulo n
         Eg: inverse_modulo(5, 11) returns -2.
        """
        if self.pgcd(a, n)>1:
            return "Inverse of "+str(a)+" mod "+str(n)+" does NOT Exist"
        return self.extended_euclidean(a, n)[1]
    
    def facto(self, n) -> list:
        """ function that returns a list of factors of n with their degrees of multiplicity.
         Args:
            n (int): an integer that one wants the prime factorization
         Returns:
            list: the list of prime factors with their order of multiplicity
         Eg: factor(60) returns [(2, 2), (3, 1), (5, 1)], or factor(33) returns [(3, 1), (11, 1)].
        """
        self.t=self.pfactors(n)
        self.factor = []
        for p in self.t:
            if not (p, self.t.count(p)) in self.factor: self.factor.append((p, self.t.count(p)))
            #for i in range(self.t.count(p)): self.t.remove(p)
        return self.factor
    
    def pfactors(self, n) -> list: 
        """ function that returns prime factors of n
         Args:
            n (int): an integer that one wants to check if square or not
         Returns:
            list: the list of prime factors
         Eg: pfactors(15) returns [3, 5], or pfactors(20) returns [2, 2, 5].
        """
        self.l = []; self.n1=n
        while self.n1%2==0:
            self.l.append(2)
            self.n1=self.n1/2
        for i in range(3, int(sqrt(self.n))+1, 2):
            while self.n1%i==0:
                self.l.append(i)
                self.n1=self.n1/i
        if self.n1>2:
            self.l.append(int(self.n1))           
        self.l.sort()
        return self.l

    def is_square(self, x:int) -> bool:
        """ function that returns True or False whether the integer x is a square or not
         Args:
            x (int): an integer that one wants to check if square or not
         Returns:
            bool: True or False whether the given integer is a square or not
         Eg: is_square(25) returns True, or is_square(10) returns False since 10 is not a square.
        """
        self.x=x
        if "/" in str(self.x):
            xp1=int(str(self.x)[:str(self.x).index("/")]); xp2=int(str(self.x)[str(self.x).index("/")+1:])
            if int(sqrt(xp1))**2==xp1 and int(sqrt(xp2))**2==xp2: return True
        if int(sqrt(int(self.x)))**2==int(self.x):
            return True
        return False

    def is_diff(self, l:list) -> bool:
        """ function that returns True if the elements of the list are different with no repetition, or False if not.
         Args:
            l (list): list of elements
         Returns:
            bool: True or False whether there are repeated elements or not
         Eg: is_diff([2, 5, 9, 10]) returns True, or is_diff([2, 5, 2, 10]) returns False since 2 is repeated two times
        """
        self.l=l
        for i in self.l:
            if self.l.count(i)>1:
                return False
        return True

    def is_same(self, l:list) -> bool:
        """ function that returns True if the elements of the list are the same, or False if not
         Args:
            l (list): list of elements
         Returns:
            bool: True or False whether the elements of the list are identical or not
         Eg: is_same([5, 5, 5]) returns True, or is_same([2, 5, 7, 2]) returns False since 2 is repeated two times
        """
        self.l=l
        if self.l.count(self.l[0])==len(self.l):
            return True
        return False


    def pair_sort(self, l)->list:
        """ function that sorts the list of tuples by first element..
         Args:
            l (list): list of tuples to sort
         Returns:
            list: the ordered tuples
         Eg: pair_sort([(4, 6), (2, 5), (5, 7)]) returns [(2, 5), (4, 6), (5, 7)]
        """
        self.l=l; self.leng=len(self.l)
        for i in range(self.leng):
            for j in range(self.leng-i-1):
                if self.l[j+1][0]<self.l[j][0]:
                    self.a = self.l[j]; self.l[j] = self.l[j+1]; self.l[j+1]=self.a
        return self.l

class H(Common):
    """
            This class implements methods used in hyperbola parametrizations H_n. 
            It provides methods related to the object H_n.
            FUNCTIONS:
                     is_fermat_solvable: property that returns True or False whether the Fermat equation has a solution or not.
                                         e.g: for n=15 returns True.  
                     info: property that returns the general info on H_n
                     is_in_H(x): function that returns True or False whether if the point x is in H or not.
                    ..... 
    """
    def __init__(self, n:int, S:str) -> None:
        super().__init__(n)
        self.n=n
        self.S=S
        if not self.S in ["Z", "Q", "Z+"] and not self.S.startswith("F"):
            print(str(self.S)+" is Not a valid algebraic structure for object H_"+str(self.n)+"(x, y), allowed are Z, Q, Z+ or Fp \n")
            sys.exit()
        if self.S=="Q":
            self.morphism="There exists a morphism over Q, f: B_"+str(self.n)+"   ------>      H_"+str(self.n)+"(Q)\n                                  (x, y)  |-----> ((x-2*"+str(self.n)+")/2*"+str(self.n)+", y/2*"+str(self.n)+") \n"
        else:
            self.morphism=""
        pass

    @property
    def is_fermat_solvable(self)->bool:
        """ property that checks the solvability of the Fermat's equation.
         Args:
            None (it considers self.n from the constructor)
         Returns:
            bool: True if x**2-y**2=n is solvable or False if not
        """
        if (self.n-6)%4==0:
            return False
        return True

    @property
    def info(self):
        """ property that prints the general info about the object H_n(x, y).
         Args:
            None
         Returns:
            int: the cardinal of H_n over self.S structure
        """
        self.start="\n________________________General Info on H_"+str(self.n)+": x^2-y^2="+str(self.n)+" over "+str(self.S)+ "________________________\n\n"
        self.form="furthermore H_"+str(self.n)+" is isomorphic to the hyperbola x^2/a^2-y^2/b^2 = 1 with a=b=sqrt("+str(self.n)+")"
        self.group="It forms a group with the additive law defined as for P+Q=(Xp*Xq+Yp*Yq, Xp*Yq+Xq*Yp),\nwith neutral element O=(1, 0).\n"
        if self.is_fermat_solvable==False:
            self.status=str(self.n)+" is not Factorizable by Fermat Method ie cannot be represented as difference of two squares"
        else:
            self.status=str(self.n)+" is Factorizable by Fermat Method ie can be written as difference of two squares"
        self.inf=str(self.start)+str(self.status)+"\n"+str(self.form)+"\n"+str(self.group)+str(self.morphism)
        print(self.inf)

    def is_in_H(self, x):
        """ function that checks whether a point is in H_n(x, y).

        Args:
            x (tuple): a point that one wants to check if it is in H_n

        Returns:
            bool: True if x is in H_n, False else
        """
        self.x=x
        if self.is_square(self.x**2-int(self.n))==True:
            return True
        return False

    def negativPoints(self, P) ->list: 
        """ function that returns negative points through Symmetry to P on Hn(x, y).
        Args:
            P (point): the point that one wants to compute negative points on Hn
        Returns:
            list: list of negative points
        """
        self.P=P;self.rep=[]
        self.rep.append((-self.P[0], -self.P[1])); self.rep.append((self.P[0], -self.P[1])); self.rep.append((-self.P[0], self.P[1])); 
        return self.rep
    
    @property
    def card(self):
        """ property that returns the cardinal of H_n(x, y).
         Args:
            None
         Returns:
            int: the cardinal of H_n over self.S structure
        """
        if self.is_fermat_solvable==False: return 0
        else:
            if self.is_prime(self.n)==True and self.n > 2:
                if self.S=="Z": return 4
                elif self.S=="Z+": return 1
                elif self.S=="Q": return "Undifined. Infinite points"
                else: return Exception("Undifined. "+str(self.S)+" not defined")
            if self.n=="2":
                if self.S=="Z" or self.S=="Z+": return 0
                elif self.S=="Q": return "Undifined."
                else: return Exception("Undifined. "+str(self.S)+" not defined")
            if self.is_diff(self.pfactors(self.n))==True and len(self.pfactors(self.n))==2 and 2 not in self.pfactors(self.n):
                if self.S=="Z": return 8
                elif self.S=="Z+": return 2
                elif self.S=="Q": return "Undifined. Infinite points"
                else: return Exception("Undifined. "+str(self.S)+" not defined")
            
    @property
    def points(self):
        """ property that returns points of H_n(x, y).
         Args:
            None
         Returns:
            list: list of points of H_n over self.S structure
        """
        if self.is_fermat_solvable==False: return "Solution is empty set"
        else:
            if self.is_prime(self.n)==True and self.n > 2:
                if self.S=="Z": return [((self.n+1)/2, (self.n-1)/2)]+self.negativPoints(((self.n+1)/2, (self.n-1)/2))
                elif self.S=="Z+": return ((self.n+1)/2, (self.n-1)/2)
                elif self.S=="Q": return "Undifined. Infinite points"
                else: return Exception("Undifined. "+str(self.S)+" not defined")
            if self.n=="2":
                if self.S=="Z" or self.S=="Z+": return "empty set"
                elif self.S=="Q": return "Undifined."
                else: return Exception("Undifined. "+str(self.S)+" not defined")
            if self.is_diff(self.pfactors(self.n))==True and len(self.pfactors(self.n))==2 and 2 not in self.pfactors(self.n):
                if self.S=="Z": return [((self.pfactors(self.n)[0]+self.pfactors(self.n)[1])/2, (self.pfactors(self.n)[1]-self.pfactors(self.n)[0])/2)]+[((self.n+1)/2, (self.n-1)/2)]+self.negativPoints(((self.n+1)/2, (self.n-1)/2))+self.negativPoints(((self.pfactors(self.n)[0]+self.pfactors(self.n)[1])/2, (self.pfactors(self.n)[1]-self.pfactors(self.n)[0])/2))
                elif self.S=="Z+": return [((self.pfactors(self.n)[0]+self.pfactors(self.n)[1])/2, (self.pfactors(self.n)[1]-self.pfactors(self.n)[0])/2)]+[((self.n+1)/2, (self.n-1)/2)]
                elif self.S=="Q": return "Undifined. Infinite points"
                else: return Exception("Undifined. "+str(self.S)+" not defined")
            else:
                pts=[]
                for x in range(int(sqrt(self.n))+1, int(self.n)**2):
                    if self.is_in_H(x)==True:
                        pts.append((x, int(sqrt(x**2-int(self.n)))))
                return pts


    def add(self, P, Q):
        """ function that adds two points on H_n(x, y).
        Args:
            P, Q (tuples): P and Q points in H_n
        Returns:
            point: result of the addition of P and Q 
        """
        xp, yp, xq, yq = P[0], P[1], Q[0], Q[1]
        if self.S in ["Z", "Q"]: return (xp*xq+yp*yq, xp*yq+xq*yp)
        elif self.S.startswith("F"): 
            self.p=int(self.S[1:])
            if "/" in str(xp):xp=int(str(xp)[:str(xp).index("/")])*self.inverse_modulo(int(str(xp)[str(xp).index("/")+1:]), self.p)
            if "/" in str(yp):yp=int(str(yp)[:str(yp).index("/")])*self.inverse_modulo(int(str(yp)[str(yp).index("/")+1:]), self.p)
            if "/" in str(xq):xq=int(str(xq)[:str(xq).index("/")])*self.inverse_modulo(int(str(xq)[str(xq).index("/")+1:]), self.p)
            if "/" in str(yq):yq=int(str(yq)[:str(yq).index("/")])*self.inverse_modulo(int(str(yq)[str(yq).index("/")+1:]), self.p)
            return ((xp*xq+yp*yq)%self.p, (xp*yq+xq*yp)%self.p)
        else:
            return Exception("Undifined. "+self.S+" Invalid struture.")

    def double(self, P):
        """ function that doubles a point on H_n(x, y).
        Args:
            P (tuple): P point in H_n
        Returns:
            tuple: result of doubling of P ie 2P 
        """
        xp, yp = P[0], P[1]
        if self.S in ["Z", "Q"]: return (xp**2+yp**2, 2*xp*yp)
        elif self.S.startswith("F"): 
            self.p=int(self.S[1:])
            if "/" in str(xp):xp=int(str(xp)[:str(xp).index("/")])*self.inverse_modulo(int(str(xp)[str(xp).index("/")+1:]), self.p)
            if "/" in str(yp):yp=int(str(yp)[:str(yp).index("/")])*self.inverse_modulo(int(str(yp)[str(yp).index("/")+1:]), self.p)
            return ((xp**2+yp**2)%self.n, (2*xp*yp)%self.n)
        else:
            return Exception("Undifined. "+self.S+" Invalid struture.")
    
    def mul(self, k, P):
        """ function that multiplies a point by a scalar on H_n(x, y).
        Args:
            k (int): a scalar
            P (tuple): P point in H_n
        Returns:
            tuple: result of multiplication of P by k ie kP 
        """
        self.k=k
        if self.k==0 : raise Exception("Invalid multiplicator k")
        self.k=bin(self.k)[2:]; self.k=str(self.k); Q=P
        for i in range(1, len(self.k)):
            Q=self.double(Q); 
            if self.k[i]=="1": Q=self.add(Q, P)
        return Q

    @property
    def plot(self, points=False):
        """ property that plots points on H_n(x, y).
        Args:
            None
        Returns:
            plot: the plot
        """
        print(self.points)
        fig = pylab.gcf()
        fig.canvas.manager.set_window_title('Hyperbo v1.0')
        x = np.linspace(-(self.n+1)-self.n, (self.n+1)+self.n)
        y = np.linspace(-(self.n+1)-self.n, self.n+self.n)
        x, y = np.meshgrid(x, y)
        plt.contour(x, y, (x**2-y**2), [self.n])
        if self.points != "Solution is empty set" and self.points != None: plt.scatter([x[0] for x in self.points ], [x[1] for x in self.points])
        plt.title("Curve of H_{} over {}".format(self.n, self.S))
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

class B(Common):
    """
            This class implements methods used in hyperbola parametrizations B_n. 
            It provides methods related to the object B_n.
            FUNCTIONS:
            ....
    """
    def __init__(self, n, S) -> None:
        super().__init__(n)
        self.n=n; self.S=S
        if not self.S in ["Z", "Q", "Z+", "Z4"] and not self.S.startswith("F"):
            print("Not valid algebraic structure, allowed Z, Q, Z+, Z4 or Fp \n")
            sys.exit()
        if self.S=="Q":
            self.p1="\nwhich inverse is            f^-1: H_"+str(self.n)+"("+str(self.S)+") ------> B_"+str(self.n)+"(x, y)\n                                   (X, Y) |-----> (2*"+str(self.n)+"*(X+1), 2*"+str(self.n)+"*Y)\n"
            self.morphism="There exists a morphism over Q, f: B_"+str(self.n)+"   ------>      H_"+str(self.n)+"(Q)\n                                  (x, y)  |-----> ((x-2*"+str(self.n)+")/2*"+str(self.n)+", y/2*"+str(self.n)+") "+self.p1
        else:
            self.morphism=f"There isn't a morphism Between H_{self.n} and B_{self.n} over {self.S}"

    
    @property
    def info(self):
        """ property that prints the general info about the object B_n(x, y).
         Args: 
            None
         Returns:
            int: the cardinal of H_n over self.S structure
        """
        self.start="\n________________________General Info on B_"+str(self.n)+": y^2=x^2-4*"+str(self.n)+"*x over "+str(self.S)+ " ________________________\n\n"
        if self.S=="Z" or self.S=="Q":self.group="It forms a group with the additive law defined as for P+Q=(1/(2*"+str(self.n)+")*((Xp-2*"+str(self.n)+")*(Xq-2*"+str(self.n)+")+Yp*Yq)+2*"+str(self.n)+", (1/(2*"+str(self.n)+")*(Yp*(Xq-2*"+str(self.n)+")+Yp*(Yq-2*"+str(self.n)+")) \n with neutral element O=(4*"+str(self.n)+", 0).\n"
        else: self.group="B_{} does not form a group over {}, But nevertheless ".format(self.n, self.S)
        if not self.S in ["Z", "Z+", "Z4"]: self.form="This structure is isomorphic to the hyperbola x^2/a^2-y^2/b^2 = 1 with a=b=sqrt("+str(self.n)+"). \n"
        else: self.form=""
        self.inf=str(self.start)+str(self.group)+"\n"+str(self.form)+str(self.morphism)
        print(self.inf)
        return ""

    def is_in_B(self, x):
        """ function that checks whether a point is in B_n(x, y).
        
        Args:
            x (tuple): a point that one wants to check if it is in B_n

        Returns:
            bool: True if x is in B_n, False else
        """
        self.x=x
        if self.is_square(self.x**2-4*int(self.n)*self.x)==True:
            return True
        return False
    
    @property
    def nbr_pointsS4(self):
        """ property that returns the number of points on B_n over Z4.

        Args:
            None: it considers self.n from the constructor

        Returns:
            int: The number of points
        """
        self.nbr_ptsS4=0
        for i in range(4*int(self.n), (int(self.n)+1)**2+1):
            if self.is_in_B(i)==True:
                self.nbr_ptsS4+=1
        return self.nbr_ptsS4
    
    @property
    def _points(self):
        """ property that returns points on B_n over Z4.

        Args:
            None, it considers self.n from the constructor

        Returns:
            list: points on B_n over Z4
        """
        self.pt=[]
        for i in range(4*int(self.n), (int(self.n)+1)**2+1):
            if self.is_in_B(i)==True:
                self.pt.append((i, int(sqrt(i**2-4*int(self.n)*i))))
        return self.pt
  
    def U(self, i:int):
        """ function that returns the i term of the sequence U(i). i represents the number of primes 

        Args:
            i (int): the integer to compute U(i)

        Returns:
            int: The element corresponding to U(i)
        E.g: U(2) returns 5
        """
        self.i=i
        if self.i==0:
            return 1
        return 3*self.U(self.i-1)-1

    def _U(self, i:int):
        self.i=i
        return 6*(2*self.U(self.i-1)-1)

    @property
    def card(self):
        """ property that returns the cardinal of B_n(x, y).
         Args:
            None
         Returns:
            int: the cardinal of B_n over self.S structure
        """
        self.fact, self.is_diff, self.is_same = self.pfactors(self.n), self.is_diff, self.is_same
        if self.is_diff(self.fact)==True and self.S=="Z4":
            return self.U(len(self.fact))
        elif self.is_diff(self.fact)==True and self.S=="Z+":
            return self.U(len(self.fact))+1
        elif self.is_diff(self.fact)==True and self.S=="Z":
            return self._U(len(self.fact))
        if self.is_same(self.fact)==True and  self.S=="Z4":
            return len(self.fact)+1
        elif self.is_same(self.fact)==True and  self.S=="Z+":
            return len(self.fact)+1+1
        elif self.is_same(self.fact)==True and  self.S=="Z":
            return 4*len(self.fact)+2
        else:
            return self.nbr_pointsS4 if self.S=="Z4" else self.nbr_pointsS4+1 if self.S=="Z+" else 4*self.nbr_pointsS4-2 if self.S=="Z" else "cardinal of B_"+str(int(self.n))+" is not defined over "+str(self.S)+" or is infinite"

    def add(self, P, Q):
        """ function that adds two points on B_n(x, y).
        Args:
            P, Q (tuples): P and Q points in B_n
        Returns:
            point: result of the addition of P and Q 
        """
        xp, yp, xq, yq = P[0], P[1], Q[0], Q[1]
        x=((xp-2*self.n)*(xq-2*self.n)+yp*yq)/(2*self.n)+2*self.n
        y=(yp*(xq-2*self.n)+yq*(xp-2*self.n))/(2*self.n)
        return (x, y)

    def double(self, P):
        """ function that doubles a point on B_n(x, y).
        Args:
            P (tuple): P point in B_n
        Returns:
            tuple: result of doubling of P ie 2P 
        """
        xp, yp = P[0], P[1]
        x=((xp-2*self.n)**2+yp**2)/(2*self.n)+2*self.n
        y=(yp*(xp-2*self.n))/self.n
        return (x, y)

    def mul(self, k, P):
        """ function that multiplies a point by a scalar on B_n(x, y).
        Args:
            k (int): a scalar
            P (tuple): P point in B_n
        Returns:
            tuple: result of multiplication of P by k ie kP 
        """
        self.k, self.P = k, P
        if self.k==0 : raise Exception("Invalid multiplicator k")
        self.k_bin=bin(self.k)[2:]
        self.k_bin=str(self.k_bin); Q=self.P 
        for i in range(1, len(self.k_bin)):
            Q = self.double(Q)
            if self.k_bin[i]=="1":
                Q=self.add(Q, self.P)
        return (Q)

    @property
    def card_sum(self):
        """ property that returns the sum S_n of cardinals on B_n(x, y).
         Args:
            None: It considers self.n from the constructor
         Returns:
            int: the sum of cardinals of B_n over self.S structure
        """
        if not self.is_diff(self.pfactors(self.n))==False:
            self.leng=len(self.pfactors(self.n))
            return self.leng/2-3*((1-3**self.leng)/4)
        return Exception(str(self.n)+" does not have all prime divisors distincts. Sum of cardinals not defined")
    
    @property
    def _productp(self) ->list:
        """ property that returns the product of prime divisors of n that make up n on B_n(x, y).
         Args:
            None: It considers self.n from the constructor
         Returns:
            list: list of tuples representing the primes that product make up self.n
        """
        _rep=[]
        for self.ki in self.pfactors(self.n):
            for self.kj in self.pfactors(self.n):
                if self.ki>self.kj and self.ki*self.kj==self.n:
                    if not (self.ki, self.kj) in _rep: _rep.append(((self.ki+self.kj)**2, self.ki**2-self.kj**2))
        return _rep

    @property
    def pointsZ4(self) ->list:
        """ property that returns points on B_n over Z4 using algebraic results on B_n.
        Args:
            None, it considers self.n from the constructor

        Returns:
            list: points on B_n over Z4
        """
        if not len(self.pfactors(self.n))>2:
            _points=[((self.k+2)*self.n+self.n/self.k, (self.k**2-1)*self.n/self.k) for self.k in self.pfactors(self.n)+[1, self.n]]+self._productp
            __points=[p for p in _points if not _points.count(p)>1]
            return self.pair_sort(__points)
        return self._points

    def negativPoints(self, l) ->list :
        """ function that returns the negative points on B_n(x, y) by symmetry from points in l.
         Args:
            l (list): list of points in B_n 
         Returns:
            list: the list of negative points of points in l 
        """
        self.l=l;self.rep=[]
        for P in self.l:
             if not P[0]==4*self.n: self.rep.append((P[0], -P[1])), self.rep.append((-P[0]+4*self.n, P[1])), self.rep.append((-P[0]+4*self.n, -P[1]))
        return self.rep
    
    @property
    def points(self):
        """ function that returns points on B_n over different algebraic structures.
         Args:
            None: It considers self.n from the constructor 
         Returns:
            list: the list of negative points of points in l 
        """
        if self.S=="Z4":
            return self.pointsZ4
        elif self.S=="Z+":
            return [(0, 0)]+self.pointsZ4
        elif self.S=="Z":
            return self.pair_sort(self.negativPoints(self.pointsZ4))+[(0, 0)]+self.pointsZ4
        elif self.S=="Q":
            return "B_"+str(self.n)+" has infinite solutions over "+str(self.S)
        else:
            return "Undifined structure. "+str(self.S)+" not defined."
    
    @property
    def plot(self, points=False):
        """ property that plots points on B_n(x, y).
        Args:
            None
        Returns:
            plot: the plot
        """
        fig = pylab.gcf()
        fig.canvas.manager.set_window_title('Hyperbo v1.0')
        x = np.linspace(-(self.n+1)**2-self.n, (self.n+1)**2+self.n)
        y = np.linspace(-(self.n+1)**2-self.n, self.n**2+self.n)
        x, y = np.meshgrid(x, y)
        plt.contour(x, y, (x**2-4*self.n*x-y**2), [0])
        if self.points != "Solution is empty set" and self.points != None: plt.scatter([x[0] for x in self.points ], [x[1] for x in self.points])
        plt.title("Curve of B_{} over {}".format(self.n, self.S))
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

    

#H(2020, "Z").info()
#p=self.pfactors(60)
#print(p)
#print(B(8, "Z").pfactors())
#print(B(60, "Z").card)
#P=(80, 40); Q=(64, 16); R=(108, 72)
#print(B(15, "Z").add(R, P))
#print(B(15, "Z").mul(100, P))
#l=[(80.0, 40.0), (108.0, 72.0), (60.0, 0.0), (256.0, 224.0), (64, 16)]
#print(B(210, "Z4").points)
#print(self.inverse_modulo(3, 7))
#P=(17/15, 8/15); Q=(5/3, 4/3); PP=(1, 12); QQ=(10, 9)
#print(H(18, "F19").add(PP, QQ))$
#print(self.is_square(25/16))
#print(H(100000000000000000000000000, "Z").card)
#print(B(15, "Z").info)
#print(self.facto(30))
#print(B(15, "Z").card_sum)
#print(Common(15).extended_euclidean(11, 13)) 
#print(Common(15).inverse_modulo(5, 11))
#print(Common(60).facto(60))
#print(Common(60).pfactors(60))
#print(Common(15).pfactors(15))
#print(Common(15).is_square(16))
#print(Common(15).pair_sort([(4, 6), (2, 5), (5, 7)]))
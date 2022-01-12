#! -*- coding: utf-8 -*-
import sys, random
from math import sqrt


class Common:
    def __init__(self) -> None:
        pass
    def rabinMiller(self, nbre:int) -> bool:
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

    def facto(self, func):
        """function that returns a list of factors of n with their degree of multiplicity"""
        self.func=func
        self.t=self.func(self.n)
        self.factor = []
        for p in self.t:
            self.factor.append((p, self.t.count(p)))
            self.t.remove(p)
        return self.factor

    #@facto
    def pfactors(self, n) ->list: 
        """function that returns prime factors of n"""
        self.l = []; self.n=n
        while self.n%2==0:
            self.l.append(2)
            self.n=self.n/2
        for i in range(3, int(sqrt(self.n))+1, 2):
            while self.n%i==0:
                self.l.append(i)
                self.n=self.n/i
        if self.n>2:
            self.l.append(int(self.n))           
        self.l.sort()
        return self.l

    def is_square(self, x:int) -> bool:
        self.x=x
        if int(sqrt(self.x))**2==x:
            return True
        return False

    def is_diff(self, l:list) -> bool:
        self.l=l
        for i in self.l:
            if self.l.count(i)>1:
                return False
        return True

    def is_same(self, l:list) -> bool:
        self.l=l
        if self.l.count(self.l[0])==len(self.l):
            return True
        return False

    def pair_sort(self, l)->list:
        self.l=l; self.leng=len(self.l)
        for i in range(self.leng):
            for j in range(self.leng-i-1):
                if self.l[j+1][0]<self.l[j][0]:
                    self.a = self.l[j]; self.l[j] = self.l[j+1]; self.l[j+1]=self.a
        return self.l

class H(Common):
    def __init__(self, n:int, S:str) -> None:
        self.n=n
        self.S=S
        if not self.S in ["Z", "Q", "Z+", "Z4"] and not self.S.startswith("F"):
            print("Not valid algebraic structure, allowed Z, Q, Z+, Z4 or Fp \n")
            sys.exit()
        if self.S=="Q":
            self.morphism="There exists a morphism over Q, f: B_"+str(self.n)+"   ------>      H_"+str(self.n)+"(Q)\n                                  (x, y)  |-----> ((x-2*"+str(self.n)+")/2*"+str(self.n)+", y/2*"+str(self.n)+") \n"
        else:
            self.morphism=""
        pass
    
    @property
    def is_fermat_solvable(self)->bool:
        if (self.n-6)%4==0:
            return False
        return True

    @property
    def info(self):
        self.start="\n________________________General Info on H_"+str(self.n)+": x^2-y^2="+str(self.n)+" over "+str(self.S)+ "________________________\n\n"
        self.form="furthermore H_"+str(self.n)+" is isomorphic to the hyperbola x^2/a^2-y^2/b^2 = 1 with a=b=sqrt("+str(self.n)+")"
        self.group="It forms a group with the additive law defined as for P+Q=(Xp*Xq+Yp*Yq, Xp*Yq+Xq*Yp),\nwith neutral element O=(1, 0).\n"
        if self.is_fermat_solvable()==False:
            self.status=str(self.n)+" is not Factorizable by Fermat Method ie cannot be represented as difference of two squares"
        else:
            self.status=str(self.n)+" is Factorizable by Fermat Method ie can be written as difference of two squares"
        self.inf=str(self.start)+str(self.status)+"\n"+str(self.form)+"\n"+str(self.group)+str(self.morphism)
        print(self.inf)
    
    @property
    def card(self):
        pass

    def add(self, P, Q):
        xp, yp, xq, yq = P[0], P[1], Q[0], Q[1]
        pass

    def mul(self, k, P):
        self.k=k; xp=P[0]; yp=P[1]
        pass

    
#H(10, "F11").info()

class B:
    def __init__(self, n, S) -> None:
        self.n=n; self.S=S
        if not self.S in ["Z", "Q", "Z+", "Z4"] and not self.S.startswith("F"):
            print("Not valid algebraic structure, allowed Z, Q, Z+, Z4 or Fp \n")
            sys.exit()
        if self.S=="Q":
            self.p1="\nwhich inverse is            f^-1: H_"+str(self.n)+"("+str(self.S)+") ------> B_"+str(self.n)+"(x, y)\n                                   (X, Y) |-----> (2*"+str(self.n)+"(X+1), 2*"+str(self.n)+"Y)\n"
            self.morphism="There exists a morphism over Q, f: B_"+str(self.n)+"   ------>      H_"+str(self.n)+"(Q)\n                                  (x, y)  |-----> ((x-2*"+str(self.n)+")/2*"+str(self.n)+", y/2*"+str(self.n)+") "+self.p1
        else:
            self.morphism=""
        pass

    @property
    def info(self):
        self.start="\n________________________General Info on B_"+str(self.n)+": y^2=x^2-4*"+str(self.n)+"*x over "+str(self.S)+ " ________________________\n\n"
        self.form="B_"+str(self.n)+" is isomorphic to the hyperbola x^2/a^2-y^2/b^2 = 1 with a=b=sqrt("+str(self.n)+")"
        self.group="It forms a group with the additive law defined as for P+Q=(1/(2*"+str(self.n)+")*((Xp-2*"+str(self.n)+")*(Xq-2*"+str(self.n)+")+Yp*Yq)+2*"+str(self.n)+", (1/(2*"+str(self.n)+")*(Yp*(Xq-2*"+str(self.n)+")+Yp*(Yq-2*"+str(self.n)+")) \nwith neutral element O=(4*"+str(self.n)+", 0).\n"
        self.inf=str(self.start)+str(self.form)+"\n"+str(self.group)+str(self.morphism)
        print(self.inf)

    def is_in_B(self, x):
        self.x=x
        if Common().is_square(self.x**2-4*int(self.n)*self.x)==True:
            return True
        return False
    
    @property
    def nbr_pointsS4(self):
        self.nbr_ptsS4=0
        for i in range(4*int(self.n), (int(self.n)+1)**2+1):
            if self.is_in_B(i)==True:
                self.nbr_ptsS4+=1
        return self.nbr_ptsS4

    @property
    def _points(self):
        self.pt=[]
        for i in range(4*int(self.n), (int(self.n)+1)**2+1):
            if self.is_in_B(i)==True:
                self.pt.append((i, int(sqrt(i**2-4*int(self.n)*i))))
        return self.pt
  

    def U(self, i:int):
        self.i=i
        if self.i==0:
            return 1
        return 3*self.U(self.i-1)-1

    def _U(self, i:int):
        self.i=i
        return 6*(2*self.U(self.i-1)-1)

    @property
    def card(self):
        """Cardinals"""
        self.fact, self.is_diff, self.is_same = Common().pfactors(self.n), Common().is_diff, Common().is_same
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
        """Adding on B_n"""
        xp, yp, xq, yq = P[0], P[1], Q[0], Q[1]
        x=((xp-2*self.n)*(xq-2*self.n)+yp*yq)/(2*self.n)+2*self.n
        y=(yp*(xq-2*self.n)+yq*(xp-2*self.n))/(2*self.n)
        return (x, y)

    def double(self, P):
        """Doubling on B_n"""
        xp, yp = P[0], P[1]
        x=((xp-2*self.n)**2+yp**2)/(2*self.n)+2*self.n
        y=(yp*(xp-2*self.n))/self.n
        return (x, y)

    def mul(self, k, P):
        """multiplying on B_n"""
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
        "Sum S_n of cardinals on B_n "
        if not Common().is_diff(Common().pfactors(self.n))==False:
            self.leng=len(Common().pfactors(self.n))
            return self.leng/2-3*(1-3**self.leng)/4
        return Exception(str(self.n)+" does not have all prime divisors distincts. Sum of cardinals not defined")
    
    @property
    def _productp(self) ->list:
        """product of prime divisors of n that make up n"""
        _rep=[]
        for self.ki in Common().pfactors(self.n):
            for self.kj in Common().pfactors(self.n):
                if self.ki>self.kj and self.ki*self.kj==self.n:
                    if not (self.ki, self.kj) in _rep: _rep.append(((self.ki+self.kj)**2, self.ki**2-self.kj**2))
        return _rep

    @property
    def pointsZ4(self) ->list:
        """points on B_n over Z4"""
        if not len(Common().pfactors(self.n))>2:
            _points=[((self.k+2)*self.n+self.n/self.k, (self.k**2-1)*self.n/self.k) for self.k in Common().pfactors(self.n)+[1, self.n]]+self._productp
            __points=[p for p in _points if not _points.count(p)>1]
            return Common().pair_sort(__points)
        return self._points

    def negativPoints(self, l) ->list :
        """negative points on B_n"""
        self.l=l;self.rep=[]
        for P in self.l:
             if not P[0]==4*self.n: self.rep.append((P[0], -P[1])), self.rep.append((-P[0]+4*self.n, P[1])), self.rep.append((-P[0]+4*self.n, -P[1]))
        return self.rep
    
    @property
    def points(self):
        """points on B_n over different algebraic structures"""
        if self.S=="Z4":
            return self.pointsZ4
        elif self.S=="Z+":
            return [(0, 0)]+self.pointsZ4
        elif self.S=="Z":
            return Common().pair_sort(self.negativPoints(self.pointsZ4))+[(0, 0)]+self.pointsZ4
        elif self.S=="Q":
            return "B_"+str(self.n)+" has infinite solutions over "+str(self.S)
        else:
            return "Undifined structure. "+str(self.S)+" not defined."

    def morphism(self, S1, S2):
        pass


#H(2020, "Z").info()
#p=Common().pfactors(60)
#print(p)
#print(B(8, "Z").pfactors())
#print(B(60, "Z").card)
#P=(80, 40); Q=(64, 16); R=(108, 72)
#print(B(15, "Z").add(R, P))
#print(B(15, "Z").mul(100, P))
#l=[(80.0, 40.0), (108.0, 72.0), (60.0, 0.0), (256.0, 224.0), (64, 16)]
print(B(21, "Z").card)

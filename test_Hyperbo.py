import unittest
from Hyperbo import H, B
from Hyperbo import Common as Co

Co=Co(5)
class TestCommon(unittest.TestCase):
    """_summary_
    Args: None
        unittest (_type_): class implementing unit tests of Hyperbola package
    """  
    def test_is_prime(self):
        self.assertTrue(Co.is_prime(19))
    
    def test_pgcd(self):
        self.assertEqual(Co.pgcd(3, 12), 3)
        
    def test_extended_euclidean(self):
        self.assertEqual(Co.extended_euclidean(3, 5), (1, 2, -1))
        
    def test_inverse_modulo(self):
        self.assertEqual(Co.inverse_modulo(3, 5), 2)
        
    def test_pfactors(self):
        self.assertEqual(Co.pfactors(20), [2, 2, 5])
    
    def test_is_square(self):
        self.assertTrue(Co.is_square(16))
           
    def test_is_diff(self):
        self.assertFalse(Co.is_diff([1, 2, 2, 5]))
        self.assertEqual(Co.is_diff([1, 2, 3]), True)
        
    def test_is_same(self):
        self.assertTrue(Co.is_same([20, 20, 20]))
        self.assertEqual(Co.is_same([1, 2, 3]), False)
        
    def test_pair_sort(self):
        self.assertEqual(Co.pair_sort([(7, 5), (2, 3)]), [(2, 3), (7, 5)])
        self.assertEqual(Co.pair_sort([(60, 0), (108, 72), (80, 40), (64, 16), (256, 224)]), [(60, 0), (64, 16), (80, 40), (108, 72), (256, 224)])
    

class TestH(unittest.TestCase):
    """
    Args:
        unittest (_type_): _description_
        This class tests all methods of H class
    """
    def test_init(self):
        __, h = ["Z", "Q", "Z+", "F2"], {}
        for _ in range(len(__)): 
            h[_]=H(123456789987654321, __[_])
            assert h[_].n==123456789987654321 and h[_].S==__[_], f"Error, something is wrong with the H class constructor for structure {__[_]}"
    def test_is_fermat_solvable(self):
        h=H(10, "Z")
        self.assertFalse(h.is_fermat_solvable)
        
    def test_is_in_H(self):
        h=H(15, "Z+")
        self.assertTrue(h.is_in_H(4))
        
    def test_negativPoints(self):
        P=(2, 3);h=H(15, "Z+")
        self.assertListEqual(h.negativPoints(P), [(-2, -3), (2, -3), (-2, 3)])
        
    def test_card(self):
        h=H(15, "Z")
        self.assertEqual(h.card, 8)
        h=H(11, "Z")
        self.assertEqual(h.card, 4)
    
    def test_points(self):
        h=H(15, "Q")
        self.assertEqual(h.points, "Undifined. Infinite points")
    
    def test_add(self):
        h=H(15, "Z")
        P, Q = (4, 1), (1, 0)
        self.assertEqual(h.add(P, Q), (4, 1))
        
    def test_double(self):
        h=H(15, "Z")
        P=(4, 1)
        self.assertEqual(h.double(P), (17, 8))
        
    def test_mul(self):
        h=H(15, "Z")
        P=(60, 0)
        self.assertEqual(h.mul(2, P), (3600, 0))
        
    def test_plot(self):
        pass
    
class TestB(unittest.TestCase):
    """
    Args:
        unittest (_type_): _description_
        This class tests all methods of B class
    """
    def test_init(self):
        __, b = ["Z", "Q", "Z+", "F2"], {}
        for _ in range(len(__)): 
            b[_]=B(123456789987654321, __[_])
            assert b[_].n==123456789987654321 and b[_].S==__[_], f"Error, something is wrong with the H class constructor for structure {__[_]}"
    
    def test_info(self):
        pass
    
    def text_is_in_B(self):
        b=B(15, "Z")
        self.assertEqual(b.is_in_B(60), True)
    
    def test_nbr_pointsS4(self):
        b=B(15, "Z4")
        self.assertEqual(b.nbr_pointsS4, 5)
        
    def test_points(self):
        b=B(15, "Z4")
        self.assertEqual(b._points, [(60, 0), (64, 16), (80, 40), (108, 72), (256, 224)])
        
    def test_U(self):
        b=B(21, "Z")
        self.assertEqual(b.U(100), 257688760366005665518230564882810636351053761001)
        self.assertEqual(b.U(2), 5)  
        
    def test__U(self):
        b=B(21, "Z")
        self.assertEqual(b._U(2), 18)
        self.assertEqual(b._U(100), 1030755041464022662072922259531242545404215044002)
        
    def test_card(self):
        b=B(15, "Z")
        self.assertEqual(b.card, 18)
        b=B(35, "Z4")
        self.assertTrue(b.card, 5)
        
    def test_add(self):
        b=B(15, "Z")
        self.assertEqual(b.add((64, 16), (60, 0)), (64, 16))
        
    def test_double(self):
        b=B(15, "Z4")
        self.assertEqual(b.double((60, 0)), (60, 0))
        
    def test_mul(self):
        b=B(33, "Z")
        self.assertEqual(b.mul(2, (4*33, 0)), (4*33, 0))
        
    def test_card_sum(self):
        Co.n=15
        b=B(15, "Z4")
        self.assertEqual(b.card_sum, 7)
        
    def test_productp(self):
        b=B(15, "Z")
        self.assertEqual(b._productp, [(64, 16)])
        
    def test_pointsZ4(self):
        b=B(15, "Z4")
        self.assertEqual(b.pointsZ4, [(60.0, 0.0), (64, 16), (80.0, 40.0), (108.0, 72.0), (256.0, 224.0)])
    
    def test_negativePoints(self):
        b=B(15,"Z")
        self.assertEqual(b.negativPoints([(60, 0), (64, 16), (80, 40), (108, 72), (256, 224)]), [(64, -16), (-4, 16), (-4, -16), (80, -40), (-20, 40), (-20, -40), (108, -72), (-48, 72), (-48, -72), (256, -224), (-196, 224), (-196, -224)])
    
    def test_points(self):
        b=B(31, "Z4")
        self.assertEqual(b.points, [(31*4, 0)])
        
    def test_plot(self):
        pass    




if __name__=='__main__':
    unittest.main()
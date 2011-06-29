from razi.functions import functions

#from nose.tools import eq_, ok_, raises, assert_almost_equal
from nose.tools import assert_almost_equal

# test with low precision, we are not testing the correctness of the
# computed descriptors, but the sanity of the functions mapping

class TestFunctions(object):

    def setup(self):
        raise NotImplementedError("Method TestFunctions.setup  must "
                                  "be implemented in subclasses.")

    def test_mw(self):
        data = (
            ('c1ccccc1', 78.114),
            ('C1CCCCC1', 84.162),
            )
        for smiles, mw in data:
            assert_almost_equal(self.engine.scalar(functions.mw(smiles)), 
                                mw, 2)

    def test_logp(self):
        data = (
            ('c1ccccc1', 1.6866),
            ('C1CCCCC1', 2.3406),
            )
        for smiles, logp in data:
            assert_almost_equal(self.engine.scalar(functions.logp(smiles)), 
                                logp, 2)

    def test_tpsa(self):
        data = (
            ('c1ccccc1', 0.0),
            ('c1ccccc1O', 20.23),
            )
        for smiles, tpsa in data:
            assert_almost_equal(self.engine.scalar(functions.tpsa(smiles)), 
                                tpsa, 1)

from razi.expression import TxtMoleculeElement
from razi.functions import functions

#from nose.tools import eq_, ok_, raises, assert_almost_equal
from nose.tools import eq_, ok_, assert_almost_equal

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

    def test_hba(self):
        data = (
            ('c1ccccc1', 0),
            ('c1ccccc1OC', 1),
            )
        for smiles, hba in data:
            eq_(self.engine.scalar(functions.hba(smiles)), hba)

    def test_hbd(self):
        data = (
            ('c1ccccc1', 0),
            ('c1ccccc1OC', 0),
            ('c1ccccc1O', 1),
            )
        for smiles, hbd in data:
            eq_(self.engine.scalar(functions.hbd(smiles)), hbd)

    def test_num_atoms(self):
        data = (
            ('c1ccccc1', 12),
            ('c1ccccc1OC', 16),
            ('c1ccccc1O', 13),
            )
        for smiles, na in data:
            eq_(self.engine.scalar(functions.num_atoms(smiles)), na)

    def test_num_hetatoms(self):
        data = (
            ('c1ccccc1', 0),
            ('c1ccccc1OC', 1),
            ('Oc1ccccc1O', 2),
            )
        for smiles, nh in data:
            eq_(self.engine.scalar(functions.num_hetatoms(smiles)), nh)

    def test_num_hvy_atoms(self):
        data = (
            ('c1ccccc1', 6),
            ('c1ccccc1OC', 8),
            ('COc1ccccc1OC', 10),
            )
        for smiles, nh in data:
            eq_(self.engine.scalar(functions.num_hvy_atoms(smiles)), nh)

    def test_num_rings(self):
        data = (
            ('c1ccccc1', 1),
            ('c1ccccc1OC1CCCC1', 2),
            ('CCCCCCCCC', 0),
            )
        for smiles, nr in data:
            eq_(self.engine.scalar(functions.num_rings(smiles)), nr)

    def test_num_rotatable_bonds(self):
        data = (
            ('c1ccccc1', 0),
            ('c1ccccc1OC1CCCC1', 2),
            ('CCCCCCCCC', 6),
            ('CC=CC', 0),
            )
        for smiles, nr in data:
            eq_(self.engine.scalar(functions.num_rotatable_bonds(smiles)), nr)

    def test_equals(self):
        ok_(self.engine.scalar(functions.equals('c1ccccc1', 'C1=CC=CC=C1')))
        ok_(self.engine.scalar(functions.equals('c1c(O)cccc1', 'OC1=CC=CC=C1')))
        ok_(not self.engine.scalar(functions.equals('CCC', 'COC')))
        

    def test_contains(self):
        ok_(self.engine.scalar(functions.contains('c1ccccc1C', 'C1=CC=CC=C1')))
        ok_(self.engine.scalar(functions.contains('c1c(O)cccc1', 'C1=CC=CC=C1')))
        ok_(not self.engine.scalar(functions.contains('CCC', 'COC')))
        
    def test_contained_in(self):
        ok_(self.engine.scalar(functions.contained_in('c1ccccc1', 
                                                      'C1=CC=C(C)C=C1')))
        ok_(self.engine.scalar(functions.contained_in('c1c(O)cccc1', 
                                                      'C1=CC=CC=C1OC')))
        ok_(not self.engine.scalar(functions.contained_in('CCC', 'COC')))
        
    def test_match(self):
        ok_(self.engine.scalar(functions.match(TxtMoleculeElement('c1ccccc1'), 
                                               'c1acccc1')))
        ok_(not self.engine.scalar(functions.match(TxtMoleculeElement('c1ncncc1'), 
                                                   'c1acccc1')))
        # shouldn't the following also work?
        # ok_(self.engine.scalar(TxtMoleculeElement('c1ccccc1').match('c1acccc1')))
        

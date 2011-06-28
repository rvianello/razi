from database import url
from sqlalchemy import create_engine
from razi.functions import functions

#from nose.tools import eq_, ok_, raises, assert_almost_equal
from nose.tools import assert_almost_equal

# test with low precision, we are not testing the mw correctness,
# but just the sanity of the function mapping

def test_mw():
    engine = create_engine(url, echo=True)
    data = (
        ('c1ccccc1', 78.114),
        ('C1CCCCC1', 84.162),
        )
    for smiles, mw in data:
        assert_almost_equal(engine.scalar(functions.mw(smiles)), mw, 2)

def test_logp():
    engine = create_engine(url, echo=True)
    data = (
        ('c1ccccc1', 1.6866),
        ('C1CCCCC1', 2.3406),
        )
    for smiles, logp in data:
        assert_almost_equal(engine.scalar(functions.logp(smiles)), logp, 2)

def test_tpsa():
    engine = create_engine(url, echo=True)
    data = (
        ('c1ccccc1', 0.0),
        ('c1ccccc1O', 20.23),
        )
    for smiles, tpsa in data:
        assert_almost_equal(engine.scalar(functions.tpsa(smiles)), tpsa, 1)


    

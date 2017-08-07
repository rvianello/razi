import unittest

from sqlalchemy import select, func
from sqlalchemy import create_engine, MetaData
from sqlalchemy import Table, Column, Integer, Index

from rdkit.Chem import AllChem as Chem

from razi.rdkit_postgresql.types import Mol


engine = create_engine('postgresql://localhost/razi-rdkit-postgresql-test')

class MolBasicTestCase(unittest.TestCase):

    def test_mol_from_smiles(self):

        rs = engine.execute(select([ func.is_valid_smiles('c1ccccc1') ]))
        self.assertTrue(rs.fetchall()[0][0])

        rs = engine.execute(select([ func.mol_from_smiles('c1ccccc1') ]))
        self.assertIsInstance(rs.fetchall()[0][0], Chem.Mol)

        rs = engine.execute(select([ func.is_valid_smiles('c1cccc') ]))
        self.assertFalse(rs.fetchall()[0][0])

        rs = engine.execute(select([ func.mol_from_smiles('c1cccc') ]))
        self.assertIsNone(rs.fetchall()[0][0], Chem.Mol)

        rs = engine.execute(select([ func.is_valid_smiles('') ]))
        self.assertTrue(rs.fetchall()[0][0])

        rs = engine.execute(select([ func.mol_from_smiles('') ]))
        self.assertIsInstance(rs.fetchall()[0][0], Chem.Mol)

    def test_mol_from_smarts(self):

        rs = engine.execute(select([ func.is_valid_smarts('c1ccc[n,c]1') ]))
        self.assertTrue(rs.fetchall()[0][0])

        rs = engine.execute(select([ func.mol_from_smarts('c1ccc[n,c]1') ]))
        self.assertIsInstance(rs.fetchall()[0][0], Chem.Mol)

        rs = engine.execute(select([ func.is_valid_smarts('c1ccc') ]))
        self.assertFalse(rs.fetchall()[0][0])

        rs = engine.execute(select([ func.mol_from_smarts('c1ccc') ]))
        self.assertIsNone(rs.fetchall()[0][0])

    def test_mol_to_smiles(self):

        rs = engine.execute(
            select([ func.mol_to_smiles(func.mol_from_smiles('c1ccccc1')) ])
        )
        self.assertEqual(rs.fetchall()[0][0], 'c1ccccc1')

        rs = engine.execute(
            select([ func.mol_to_smiles(func.qmol('c1cccc[n,c]1')) ])
        )
        self.assertEqual(rs.fetchall()[0][0], 'c1ccncc1')

    def test_mol_to_smarts(self):

        rs = engine.execute(
            select([ func.mol_to_smarts(func.mol_from_smiles('c1ccccc1')) ])
        )
        self.assertEqual(
            rs.fetchall()[0][0],
            '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1')


metadata = MetaData(engine)

molecules = Table(
    'molecules', metadata,
    Column('id', Integer, primary_key=True),
    Column('m', Mol), Index('molecules_m', 'm', postgresql_using='gist'),
    )

class MolTestCase(unittest.TestCase):

    def setUp(self):
        metadata.create_all()

    def tearDown(self):
        metadata.drop_all()

    def test_this(self):
        pass

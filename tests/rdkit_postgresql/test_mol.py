import os
import unittest
import json

from sqlalchemy import select, func
from sqlalchemy import MetaData
from sqlalchemy import Table, Column, Integer, Index

from rdkit.Chem import AllChem as Chem

from razi.rdkit_postgresql.types import Mol

from .database import engine


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

    def test_mol_from_ctab(self):

        rs = engine.execute(select([ func.is_valid_ctab('xyz') ]))
        self.assertFalse(rs.fetchall()[0][0])

        ctab = '''chiral1.mol
  ChemDraw04200416412D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0141    0.0553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8109    0.0553    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4266    0.7697    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.0141   -0.7697    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8109   -0.1583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  1
  1  5  1  0
M  END'''

        rs = engine.execute(select([ func.is_valid_ctab(ctab) ]))
        self.assertTrue(rs.fetchall()[0][0])

        rs = engine.execute(select([ func.mol_from_ctab(ctab) ]))
        self.assertIsInstance(rs.fetchall()[0][0], Chem.Mol)

    def test_mol_to_ctab(self):

        rs = engine.execute(select([ func.mol_to_ctab(func.mol('CCC')) ]))
        self.assertEqual(rs.fetchall()[0][0], '''
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 C 2.598076 -0.000000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
''')

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
        self.assertEqual(rs.fetchall()[0][0], '*1ccccc1')

    def test_mol_to_smarts(self):

        rs = engine.execute(
            select([ func.mol_to_smarts(func.mol_from_smiles('c1ccccc1')) ])
        )
        self.assertEqual(
            rs.fetchall()[0][0],
            '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1')

    def test_mol_samestruct(self):

        rs = engine.execute(
            select([ func.mol('Cc1ccccc1') == func.mol('c1ccccc1C') ])
        )
        self.assertTrue(rs.fetchall()[0][0])

        rs = engine.execute(
            select([ func.mol('Cc1ccccc1') == func.mol('c1cnccc1C') ])
        )
        self.assertFalse(rs.fetchall()[0][0])


metadata = MetaData(engine)

molecules = Table(
    'molecules', metadata,
    Column('id', Integer, primary_key=True),
    Column('m', Mol), Index('molecules_m', 'm', postgresql_using='gist'),
    )


_test_dir = os.path.dirname(__file__)
for _ in range(2):
    _test_dir = os.path.dirname(_test_dir)
_test_dir = os.path.join(_test_dir, 'test_data')

data_filepath = os.path.join(_test_dir, 'chembl_sample.json')


class MolInsertTestCase(unittest.TestCase):

    def setUp(self):
        metadata.create_all()

    def tearDown(self):
        metadata.drop_all()

    def test_mol_insert(self):

        self.assertTrue(os.path.exists(data_filepath))

        with open(data_filepath, 'rt') as f:
            sample_data = json.load(f)

        engine.execute(
            molecules.insert(),
            [ { 'm': Chem.MolFromSmiles(s)} for _, s in sample_data ],
            )

        rs = engine.execute(select([func.count()]).select_from(molecules))

        self.assertEqual(len(sample_data), rs.fetchall()[0][0])


class MolTestCase(unittest.TestCase):

    def setUp(self):
        metadata.create_all()

        with open(data_filepath, 'rt') as f:
            sample_data = json.load(f)

        engine.execute(
            molecules.insert(),
            [ { 'm': Chem.MolFromSmiles(s)} for _, s in sample_data ],
            )

    def tearDown(self):
        metadata.drop_all()

    def test_mol_ops(self):

        molecules_count = select([func.count()]).select_from(molecules)

        rs = engine.execute(molecules_count.where(
            molecules.c.m.hassubstruct(Chem.MolFromSmiles('c1ccccc1'))
        ))
        self.assertEqual(830, rs.fetchall()[0][0])

        rs = engine.execute(molecules_count.where(
            molecules.c.m.hassubstruct(Chem.MolFromSmiles('c1ccncc1'))
        ))
        self.assertEqual(196, rs.fetchall()[0][0])

        rs = engine.execute(molecules_count.where(
            molecules.c.m.hassubstruct(func.mol_from_smarts('c1ccccc1'))
        ))
        self.assertEqual(830, rs.fetchall()[0][0])

        rs = engine.execute(molecules_count.where(
            molecules.c.m.hassubstruct(func.mol_from_smarts('c1ccncc1'))
        ))
        self.assertEqual(197, rs.fetchall()[0][0])

        rs = engine.execute(molecules_count.where(
            molecules.c.m.hassubstruct(func.mol_from_smarts('c1cc[c,n]cc1'))
        ))
        self.assertEqual(868, rs.fetchall()[0][0])

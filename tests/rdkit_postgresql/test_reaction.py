import os
import unittest

from sqlalchemy import MetaData
from sqlalchemy import select, func, bindparam
from sqlalchemy import Table, Column, Integer, Index

from rdkit.Chem import AllChem as Chem

from razi.rdkit_postgresql.types import Reaction

from .database import engine


sys_metadata = MetaData(engine)
pg_settings = Table('pg_settings', sys_metadata,
                    autoload=True, autoload_with=engine)


class ReactionBasicTestCase(unittest.TestCase):

    def test_reaction(self):

        # dummy query to ensure that GUC params are initialized
        _ = engine.execute(select([ func.rdkit_version() ])).fetchall()

        stmt = pg_settings.update()
        stmt = stmt.where(pg_settings.c.name == bindparam('param'))
        stmt = stmt.values(setting=bindparam('value'))

        engine.execute(stmt, [
            {'param': 'rdkit.ignore_reaction_agents', 'value': False},
            {'param': 'rdkit.agent_FP_bit_ratio', 'value': 0.2},
            {'param': 'rdkit.difference_FP_weight_agents', 'value': 1},
            {'param': 'rdkit.difference_FP_weight_nonagents', 'value': 10},
            {'param': 'rdkit.move_unmmapped_reactants_to_agents',
             'value': True},
            {'param': 'rdkit.threshold_unmapped_reactant_atoms',
             'value': 0.2},
            {'param': 'rdkit.init_reaction', 'value': True},
            ])

        rs = engine.execute(
            select([func.reaction_from_smiles('c1ccccc1>>c1cccnc1')])
            )
        self.assertEqual(rs.fetchall()[0][0], 'c1ccccc1>>c1ccncc1')

        rs = engine.execute(
            select([func.reaction_from_smiles('c1ccccc1>CC(=O)O>c1cccnc1')])
            )
        self.assertEqual(rs.fetchall()[0][0], 'c1ccccc1>CC(=O)O>c1ccncc1')

        rs = engine.execute(
            select([func.reaction_from_smarts('[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]')])
            )
        self.assertEqual(rs.fetchall()[0][0], 'c([cH:4][cH:3][cH:2][cH2:1])[cH2:5]>CC(=O)O>n([cH:4][cH:3][cH:2][cH2:1])[cH2:5]')

        rs = engine.execute(
            select([func.reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]')])
            )
        self.assertEqual(rs.fetchall()[0][0], 'c([cH:4][cH:3][cH:2][cH2:1])[cH2:5]>CC(=O)O.FC(F)F>n([cH:4][cH:3][cH:2][cH2:1])[cH2:5]')

        rs = engine.execute(
            select([func.reaction_from_smarts('c1ccc[n,c]c1>>c1nccnc1')])
            )
        self.assertEqual(rs.fetchall()[0][0], '*1ccccc1>>c1cnccn1')

        rs = engine.execute(
            select([func.reaction_to_smiles(func.reaction_from_smiles('c1ccccc1>>c1cccnc1'))])
            )
        self.assertEqual(rs.fetchall()[0][0], 'c1ccccc1>>c1ccncc1')

        rs = engine.execute(
            select([func.reaction_from_ctab('''$RXN

      RDKit

  1  1
$MOL

     RDKit

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  6  1  4  0
M  END
$MOL

     RDKit

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  6  1  4  0
M  END
''')])
        )
        self.assertEqual(rs.fetchall()[0][0], 'c1ccccc1>>c1ccncc1')

        rs = engine.execute(
            select([func.reaction_numreactants(func.reaction_from_smiles('[Cl].c1ccccc1>>c1cccnc1.[OH2]'))])
            )
        self.assertEqual(rs.fetchall()[0][0], 2)

        rs = engine.execute(
            select([func.reaction_numproducts(func.reaction_from_smiles('[Cl].c1ccccc1>>c1cccnc1.[OH2]'))])
            )
        self.assertEqual(rs.fetchall()[0][0], 2)

        rs = engine.execute(
            select([func.reaction_numagents(func.reaction_from_smiles('[Cl].c1ccccc1>CC(=O)O.[Na+]>c1cccnc1.[OH2]'))])
            )
        self.assertEqual(rs.fetchall()[0][0], 2)

        rs = engine.execute(
            select([func.reaction_numagents(func.reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]'))])
            )
        self.assertEqual(rs.fetchall()[0][0], 2)

        engine.execute(stmt, [
            {'param': 'rdkit.move_unmmapped_reactants_to_agents',
             'value': False},
            ])

        rs = engine.execute(
            select([func.reaction_numagents(func.reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]'))])
            )
        self.assertEqual(rs.fetchall()[0][0], 1)

        engine.execute(stmt, [
            {'param': 'rdkit.move_unmmapped_reactants_to_agents',
             'value': True},
            {'param': 'rdkit.threshold_unmapped_reactant_atoms',
             'value': 0.9},
            ])

        rs = engine.execute(
            select([func.reaction_numagents(func.reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]'))])
            )
        self.assertEqual(rs.fetchall()[0][0], 3)

        engine.execute(stmt, [
            {'param': 'rdkit.threshold_unmapped_reactant_atoms',
             'value': 0.2},
            ])

        rs = engine.execute(
            select([ func.reaction('c1ccccc1>>c1cccnc1') ==
                     func.reaction('c1ccccc1>>c1cccnc1') ])
            )
        self.assertEqual(rs.fetchall()[0][0], True)

        rs = engine.execute(
            select([ func.reaction('c1ccccc1>>c1cccnc1') ==
                     func.reaction('c1ccccc1>>c1cncnc1') ])
            )
        self.assertEqual(rs.fetchall()[0][0], False)


metadata = MetaData(engine)

reactions = Table(
    'reactions', metadata,
    Column('id', Integer, primary_key=True),
    Column('rxn', Reaction),
    Index('reactions_rxn', 'rxn', postgresql_using='gist'),
    )

reactions_unchanged = Table(
    'reactions_unchanged', metadata,
    Column('id', Integer, primary_key=True),
    Column('rxn', Reaction),
    #Index('reactions_unchanged_rxn', 'rxn', postgresql_using='gist'),
    )


_test_dir = os.path.dirname(__file__)
for _ in range(2):
    _test_dir = os.path.dirname(_test_dir)
_test_dir = os.path.join(_test_dir, 'test_data')

data_filepath = os.path.join(_test_dir, 'reaction_test_data.out.rsmi')


class ReactionTestCase(unittest.TestCase):

    def setUp(self):
        metadata.create_all()

    def tearDown(self):
        metadata.drop_all()

    def test_reaction(self):

        self.assertTrue(os.path.exists(data_filepath))

        with open(data_filepath, 'rt') as f:
            sample_data = [line.split() for line in f]

        # dummy query to ensure that GUC params are initialized
        _ = engine.execute(select([ func.rdkit_version() ])).fetchall()

        stmt = pg_settings.update()
        stmt = stmt.where(pg_settings.c.name == bindparam('param'))
        stmt = stmt.values(setting=bindparam('value'))

        engine.execute(stmt, [
            {'param': 'rdkit.ignore_reaction_agents', 'value': False},
            {'param': 'rdkit.agent_FP_bit_ratio', 'value': 0.2},
            {'param': 'rdkit.difference_FP_weight_agents', 'value': 1},
            {'param': 'rdkit.difference_FP_weight_nonagents', 'value': 10},
            {'param': 'rdkit.move_unmmapped_reactants_to_agents',
             'value': True},
            {'param': 'rdkit.threshold_unmapped_reactant_atoms',
             'value': 0.2},
            {'param': 'rdkit.init_reaction', 'value': True},
            ])

        engine.execute(
            reactions.insert(),
            [ { 'rxn': rsmiles, } for _, rsmiles in sample_data ],
            )

        rs = engine.execute(select([func.count()]).select_from(reactions))
        sz = rs.fetchall()[0][0]
        self.assertEqual(sz, len(sample_data))

        engine.execute(stmt, [
            {'param': 'rdkit.move_unmmapped_reactants_to_agents',
             'value': False},
            ])

        engine.execute(
            reactions_unchanged.insert(),
            [ { 'rxn': rsmiles, } for _, rsmiles in sample_data ],
            )

        engine.execute(stmt, [
            {'param': 'rdkit.move_unmmapped_reactants_to_agents',
             'value': True},
            ])

        rs = engine.execute(
            select([func.count()]).select_from(reactions_unchanged))
        sz = rs.fetchall()[0][0]
        self.assertEqual(sz, len(sample_data))

        rs = engine.execute(
            select([
                func.sum(func.reaction_numreactants(reactions.c.rxn))
                ])
            )
        sum_reactants = rs.fetchall()[0][0]
        self.assertEqual(sum_reactants, 1898)

        rs = engine.execute(
            select([
                func.sum(func.reaction_numreactants(reactions_unchanged.c.rxn))
                ])
            )
        sum_reactants = rs.fetchall()[0][0]
        self.assertEqual(sum_reactants, 3517)

        rs = engine.execute(
            select([
                func.sum(func.reaction_numproducts(reactions.c.rxn))
                ])
            )
        sum_products = rs.fetchall()[0][0]
        self.assertEqual(sum_products, 1157)

        rs = engine.execute(
            select([
                func.sum(func.reaction_numproducts(reactions_unchanged.c.rxn))
                ])
            )
        sum_products = rs.fetchall()[0][0]
        self.assertEqual(sum_products, 1157)

        rs = engine.execute(
            select([
                func.sum(func.reaction_numagents(reactions.c.rxn))
                ])
            )
        sum_agents = rs.fetchall()[0][0]
        self.assertEqual(sum_agents, 2528)

        rs = engine.execute(
            select([
                func.sum(func.reaction_numagents(reactions_unchanged.c.rxn))
                ])
            )
        sum_agents = rs.fetchall()[0][0]
        self.assertEqual(sum_agents, 909)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstruct('c1ccccc1>>c1ccncc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 47)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstruct('c1cnccc1>>c1ccccc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 50)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstruct('c1ccccc1>>')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 667)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstruct('c1cnccc1>>')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 83)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstruct('>>c1ccncc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 79)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstruct('>>c1ccccc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 650)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstructfp('c1ccccc1>>c1ccncc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 47)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstructfp('c1cnccc1>>c1ccccc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 50)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstructfp('c1ccccc1>>')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 667)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstructfp('c1cnccc1>>')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 83)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstructfp('>>c1ccncc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 79)

        rs = engine.execute(
            select([func.count()]).select_from(reactions).where(
                reactions.c.rxn.hassubstructfp('>>c1ccccc1')
            )
        )
        self.assertEqual(rs.fetchall()[0][0], 650)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',1),
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',1))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',1),
                    func.reaction_difference_fp('c1ncccc1>>c1ncncc1',1))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.636363636363636)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1),
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1),
                    func.reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',1))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.603448275862069)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',2),
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',2))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',2),
                    func.reaction_difference_fp('c1ncccc1>>c1ncncc1',2))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.2)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2),
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2),
                    func.reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',2))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.2)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',3),
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',3))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>>c1ccncc1',3),
                    func.reaction_difference_fp('c1ncccc1>>c1ncncc1',3))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.454545454545455)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3),
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3),
                    func.reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',3))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.444933920704846)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_structural_bfp('c1ccccc1>>c1ccncc1',1),
                    func.reaction_structural_bfp('c1ccccc1>>c1ccncc1',1))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_structural_bfp('c1ccccc1>>c1ccncc1',1),
                    func.reaction_structural_bfp('c1ncccc1>>c1ncncc1',1))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.620689655172414)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1),
                    func.reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1))
            ])
        )
        self.assertEqual(rs.fetchall()[0][0], 1.0)

        rs = engine.execute(
            select([
                func.tanimoto_sml(
                    func.reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1),
                    func.reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',1))
            ])
        )
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.514285714285714)

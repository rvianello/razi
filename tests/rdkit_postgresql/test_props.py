import unittest

from sqlalchemy import select, func

from rdkit.Chem import AllChem as Chem

from razi import rdkit_postgresql

from .database import engine


class PropsTestCase(unittest.TestCase):

    def test_props(self):

        rs = engine.execute(select([ func.mol_amw(func.mol('c1ccccc1')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 78.114)

        rs = engine.execute(select([ func.mol_logp(func.mol('c1ccccc1')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 1.6866)

        rs = engine.execute(select([ func.mol_hba(func.mol('c1ccccc1')) ]))
        self.assertEqual(rs.fetchall()[0][0], 0)

        rs = engine.execute(select([ func.mol_hbd(func.mol('c1ccccc1')) ]))
        self.assertEqual(rs.fetchall()[0][0], 0)


        rs = engine.execute(select([ func.mol_chi0n(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 3.8339648)

        rs = engine.execute(select([ func.mol_chi1n(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 2.1342905)

        rs = engine.execute(select([ func.mol_chi2n(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 1.3355491)

        rs = engine.execute(select([ func.mol_chi3n(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.7561936)

        rs = engine.execute(select([ func.mol_chi4n(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.4279941)


        rs = engine.execute(select([ func.mol_chi0v(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 3.8339648)

        rs = engine.execute(select([ func.mol_chi1v(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 2.1342905)

        rs = engine.execute(select([ func.mol_chi2v(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 1.3355491)

        rs = engine.execute(select([ func.mol_chi3v(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.7561936)

        rs = engine.execute(select([ func.mol_chi4v(func.mol('c1ccccc1O')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 0.4279941)


        rs = engine.execute(select([ func.mol_kappa1(func.mol('C12CC2C3CC13')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 2.34375)

        rs = engine.execute(select([ func.mol_kappa2(func.mol('CC(C)C1CCC(C)CCC1')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 4.132653)

        rs = engine.execute(select([ func.mol_kappa3(func.mol('CC(C)C1CCC(C)CCC1')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 2.8444445)


        rs = engine.execute(select([ func.mol_numspiroatoms(func.mol('CC1(C)CC2(C)CCC1(C)CC2')) ]))
        self.assertEqual(rs.fetchall()[0][0], 0)

        rs = engine.execute(select([ func.mol_numbridgeheadatoms(func.mol('CC1(C)CC2(C)CCC1(C)CC2')) ]))
        self.assertEqual(rs.fetchall()[0][0], 2)


        rs = engine.execute(select([ func.mol_formula(func.mol('[2H]C([3H])O')) ]))
        self.assertEqual(rs.fetchall()[0][0], 'CH4O')


        rs = engine.execute(select([ func.mol_numrotatablebonds(func.mol('CCCC')) ]))
        self.assertEqual(rs.fetchall()[0][0], 1)

        rs = engine.execute(select([ func.mol_numheavyatoms(func.mol('CCC')) ]))
        self.assertEqual(rs.fetchall()[0][0], 3)

        rs = engine.execute(select([ func.mol_numatoms(func.mol('CCC')) ]))
        self.assertEqual(rs.fetchall()[0][0], 11)

        rs = engine.execute(select([ func.mol_numheteroatoms(func.mol('CCO')) ]))
        self.assertEqual(rs.fetchall()[0][0], 1)


        rs = engine.execute(select([ func.mol_tpsa(func.mol('CCO')) ]))
        self.assertAlmostEqual(rs.fetchall()[0][0], 20.23)


        rs = engine.execute(select([ func.mol_numrings(func.mol('CCO')) ]))
        self.assertEqual(rs.fetchall()[0][0], 0)


        rs = engine.execute(select([ func.mol_murckoscaffold(func.mol('c1ccccc1CCO')) ]))
        self.assertEqual(rs.fetchall()[0][0], 'c1ccccc1')


        rs = engine.execute(select([ func.mol_to_svg(func.mol('CCO')) ]))
        self.assertEqual(rs.fetchall()[0][0], '''<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='250px' height='200px' viewBox='0 0 250 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='250.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 11.4,121.0 L 107.5,65.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 107.5,65.5 L 147.7,88.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 147.7,88.7 L 187.8,111.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 190.6 121.1
Q 190.6 114.3, 194.0 110.5
Q 197.3 106.7, 203.6 106.7
Q 209.9 106.7, 213.2 110.5
Q 216.6 114.3, 216.6 121.1
Q 216.6 128.0, 213.2 131.9
Q 209.8 135.8, 203.6 135.8
Q 197.4 135.8, 194.0 131.9
Q 190.6 128.0, 190.6 121.1
M 203.6 132.6
Q 207.9 132.6, 210.2 129.7
Q 212.6 126.8, 212.6 121.1
Q 212.6 115.6, 210.2 112.8
Q 207.9 109.9, 203.6 109.9
Q 199.3 109.9, 196.9 112.7
Q 194.6 115.5, 194.6 121.1
Q 194.6 126.8, 196.9 129.7
Q 199.3 132.6, 203.6 132.6
' fill='#FF0000'/>
<path  class='atom-2' d='M 220.0 107.0
L 223.8 107.0
L 223.8 119.1
L 238.3 119.1
L 238.3 107.0
L 242.2 107.0
L 242.2 135.4
L 238.3 135.4
L 238.3 122.3
L 223.8 122.3
L 223.8 135.4
L 220.0 135.4
L 220.0 107.0
' fill='#FF0000'/>
</svg>
''')

    maxDiff = None

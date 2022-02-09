import unittest

from sqlalchemy import select, func

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
<path class='bond-0 atom-0 atom-1' d='M 11.4,118.7 L 107.8,63.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 107.8,63.1 L 148.6,86.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 148.6,86.6 L 189.4,110.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='192.2' y='138.7' class='atom-2' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' >O</text>
<text x='219.8' y='138.7' class='atom-2' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' >H</text>
</svg>
''')

    maxDiff = None

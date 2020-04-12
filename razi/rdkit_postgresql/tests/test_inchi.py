import unittest

from sqlalchemy import select, func

from rdkit.Chem import AllChem as Chem

from razi import rdkit_postgresql

from .database import engine


class InChITestCase(unittest.TestCase):

    def test_inchi_out(self):

        rs = engine.execute(select([ func.mol_inchi(func.mol('c1ccccc1')) ]))
        self.assertEqual(rs.fetchall()[0][0], 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H')

        rs = engine.execute(select([ func.mol_inchikey(func.mol('c1ccccc1')) ]))
        self.assertEqual(rs.fetchall()[0][0], 'UHOVQNZJYSORNB-UHFFFAOYSA-N')

        rs = engine.execute(select([ func.mol_inchi(func.mol('Cc1cc(C)[n+]c(C)c1')) ]))
        self.assertEqual(rs.fetchall()[0][0], 'InChI=1S/C8H11N/c1-6-4-7(2)9-8(3)5-6/h4-5H,1-3H3/q+1')

        rs = engine.execute(select([ func.mol_inchikey(func.mol('Cc1cc(C)[n+]c(C)c1')) ]))
        self.assertEqual(rs.fetchall()[0][0], 'WDQXTFQYHUZLPX-UHFFFAOYSA-N')

import unittest

from sqlalchemy import create_engine, MetaData
from sqlalchemy import Table, Column, Integer, Index

from razi.rdkit_postgresql.types import Mol


engine = create_engine('postgresql://localhost/razi-rdkit-postgresql-test')
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

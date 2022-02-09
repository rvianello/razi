from sqlalchemy import func
from sqlalchemy.types import UserDefinedType
from sqlalchemy.dialects.postgresql.base import ischema_names

from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect

from . import comparator


class Mol(UserDefinedType):

    cache_ok = True
    comparator_factory = comparator.MolComparator

    def get_col_spec(self, **kw):
        return 'mol'

    def bind_processor(self, dialect):
        def process(value):
            # convert the Molecule instance to the value used by the
            # db driver
            if isinstance(value, Chem.Mol):
                value = memoryview(value.ToBinary())
            elif isinstance(value, str):
                value = memoryview(Chem.MolFromSmiles(value).ToBinary())
            return value
        return process

    def bind_expression(self, bindvalue):
        return func.mol_from_pkl(bindvalue)

    def column_expression(self, col):
        return func.mol_to_pkl(col, type_=self)

    def result_processor(self, dialect, coltype):
        def process(value):
            if value is None:
                return value
            elif isinstance(value, memoryview):
                return Chem.Mol(bytes(value))
            else:
                raise RuntimeError(
                    "Unexpected row value type for a Mol instance")
        return process


class QMol(UserDefinedType):

    cache_ok = True
    comparator_factory = comparator.QMolComparator

    def get_col_spec(self, **kw):
        return 'qmol'

    def bind_processor(self, dialect):
        def process(value):
            # convert the Molecule instance to the value used by the
            # db driver
            if isinstance(value, Chem.Mol):
                return Chem.MolToSmarts(value)
            elif isinstance(value, str):
                return value
            else:
                raise RuntimeError(
                    'Unexpected query value type for QMol column')

        return process


class Bfp(UserDefinedType):

    cache_ok = True
    comparator_factory = comparator.BfpComparator

    def get_col_spec(self, **kw):
        return 'bfp'

    def bind_processor(self, dialect):
        def process(value):
            if isinstance(value, ExplicitBitVect):
                value = memoryview(DataStructs.BitVectToBinaryText(value))
            return value
        return process

    def bind_expression(self, bindvalue):
        return func.bfp_from_binary_text(bindvalue)

    def column_expression(self, col):
        return func.bfp_to_binary_text(col, type_=self)

    def result_processor(self, dialect, coltype):
        def process(value):
            if value is None:
                return value
            elif isinstance(value, memoryview):
                return DataStructs.CreateFromBinaryText(bytes(value))
            else:
                raise RuntimeError(
                    "Unexpected row value type for a Bfp instance")
        return process


class Sfp(UserDefinedType):

    cache_ok = True
    comparator_factory = comparator.SfpComparator

    def get_col_spec(self, **kw):
        return 'sfp'


class Reaction(UserDefinedType):

    cache_ok = True
    comparator_factory = comparator.ReactionComparator

    def get_col_spec(self, **kw):
        return 'reaction'


# register RDKit types to PostgreSQL dialect, so that SQLAlchemy relfect infrastructure works properly
ischema_names.update({
    'mol': Mol,
    'qmol': QMol,
    'bfp': Bfp,
    'sfp': Sfp,
    'reaction': Reaction,
})

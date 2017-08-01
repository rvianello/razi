from sqlalchemy import func
from sqlalchemy.types import UserDefinedType, TypeEngine

from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect

from . import comparator

class Mol(UserDefinedType):

    comparator_factory = comparator.MolComparator

    def get_col_spec(self, **kw):
        return 'mol'

    def bind_processor(self, dialect):
        def process(value):
            # convert the Molecule instance to the value used by the
            # db driver
            if isinstance(value, Chem.Mol):
                value = memoryview(value.ToBinary())
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

    comparator_factory = comparator.QMolComparator

    def get_col_spec(self, **kw):
        return 'qmol'


class Bfp(UserDefinedType):

    comparator_factory = comparator.FpComparator

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

    def get_col_spec(self, **kw):
        return 'sfp'


class Reaction(UserDefinedType):

    def get_col_spec(self, **kw):
        return 'reaction'

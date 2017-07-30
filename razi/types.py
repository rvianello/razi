from sqlalchemy import func
from sqlalchemy.types import UserDefinedType, TypeEngine

from rdkit.Chem import AllChem as Chem


class Mol(UserDefinedType):

    def get_col_spec(self, **kw):
        return 'mol'

    def bind_expression(self, bindvalue):
        return func.mol_from_pkl(bindvalue)

    def column_expression(self, col):
        return func.mol_to_pkl(col, type_=self)

    def bind_processor(self, dialect):
        def process(value):
            # convert the Molecule instance to the value used by the
            # db driver
            if isinstance(value, Chem.Mol):
                value = memoryview(value.ToBinary())
            return value
        return process

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

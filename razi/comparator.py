from sqlalchemy import types as sqltypes
from sqlalchemy.types import UserDefinedType
from sqlalchemy.sql import operators

class MolComparator(UserDefinedType.Comparator):

    def hassubstruct(self, other):
        return self.operate(
            operators.custom_op('@>'), other, result_type=sqltypes.Boolean
            )

    def issubstruct(self, other):
        return self.operate(
            operators.custom_op('<@'), other, result_type=sqltypes.Boolean
            )

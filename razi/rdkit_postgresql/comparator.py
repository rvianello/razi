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


class QMolComparator(UserDefinedType.Comparator):

    def issubstruct(self, other):
        return self.operate(
            operators.custom_op('<@'), other, result_type=sqltypes.Boolean
            )


class FpComparator(UserDefinedType.Comparator):

    def tanimoto_sml(self, other):
        return self.operate(
            operators.custom_op('%'), other, result_type=sqltypes.Boolean
            )

    def dice_sml(self, other):
        return self.operate(
            operators.custom_op('#'), other, result_type=sqltypes.Boolean
            )


class ReactionComparator(UserDefinedType.Comparator):

    def hassubstruct(self, other):
        return self.operate(
            operators.custom_op('@>'), other, result_type=sqltypes.Boolean
            )

    def hassubstructfp(self, other):
        return self.operate(
            operators.custom_op('?>'), other, result_type=sqltypes.Boolean
            )

    def issubstruct(self, other):
        return self.operate(
            operators.custom_op('<@'), other, result_type=sqltypes.Boolean
            )

    def issubstructfp(self, other):
        return self.operate(
            operators.custom_op('<?'), other, result_type=sqltypes.Boolean
            )

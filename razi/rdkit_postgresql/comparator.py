from sqlalchemy import func
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

    def __eq__(self, other):
        return self.operate(
            operators.custom_op('@='), other, result_type=sqltypes.Boolean
            )


class QMolComparator(UserDefinedType.Comparator):

    def issubstruct(self, other):
        return self.operate(
            operators.custom_op('<@'), other, result_type=sqltypes.Boolean
            )


class BfpComparator(UserDefinedType.Comparator):

    def tanimoto_sml(self, other):
        return self.operate(
            operators.custom_op('%'), other, result_type=sqltypes.Boolean
            )

    def dice_sml(self, other):
        return self.operate(
            operators.custom_op('#'), other, result_type=sqltypes.Boolean
            )


class SfpComparator(UserDefinedType.Comparator):

    def tanimoto_sml(self, other):
        return self.operate(
            operators.custom_op('%'), other, result_type=sqltypes.Boolean
            )

    def dice_sml(self, other):
        return self.operate(
            operators.custom_op('#'), other, result_type=sqltypes.Boolean
            )

    def __add__(self, other):
        return func.add(self.expr, other)

    def __sub__(self, other):
        return func.subtract(self.expr, other)


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

    def __eq__(self, other):
        return self.operate(
            operators.custom_op('@='), other, result_type=sqltypes.Boolean
            )

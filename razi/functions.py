from sqlalchemy.sql.expression import ClauseElement
from sqlalchemy.sql import expression
try:
    from sqlalchemy.sql.functions import Function
except ImportError:
    # SQLAlchemy<0.9
    Function = expression.Function
from sqlalchemy.ext.compiler import compiles
from sqlalchemy import literal

import types


def parse_clause(clause, compiler, translate_string=None):
    """This method is used to translate a clause element (molecules,
    functions, ..).

    According to the type of the clause, a conversion to the database
    molecule type is added or the column clause (column name) or the cascaded
    clause element is returned.

    """
    from razi.chemtypes import Molecule
    from razi.expression import MoleculeElement, TxtMoleculeElement
    from razi.expression import PersistentMoleculeElement
    from razi.expression import QMoleculeElement, TxtQMoleculeElement

    if hasattr(clause, '__clause_element__'):
        # for example a column name
        return clause.__clause_element__()
    elif isinstance(clause, ClauseElement):
        # for cascaded clause elements, like other functions
        return clause
    elif isinstance(clause, MoleculeElement):
        if isinstance(clause, TxtMoleculeElement):
            return clause
        elif isinstance(clause, PersistentMoleculeElement):
            return literal(clause.desc, Molecule)
        raise TypeError
    elif isinstance(clause, QMoleculeElement):
        if isinstance(clause, TxtQMoleculeElement):
            return clause
        raise TypeError
    elif isinstance(clause, basestring) and translate_string:
        return translate_string(clause)

    # for raw parameters
    return literal(clause)


def _get_function(element, compiler, params, within_column_clause):
    """For elements of type BaseFunction, the database specific function data
    is looked up and a executable sqlalchemy.sql.expression.Function object
    is returned.
    """
    from razi.dialect import DialectManager
    chemical_dialect = DialectManager.get_chemical_dialect(compiler.dialect)
    function_data = chemical_dialect.get_function(element.__class__)

    if isinstance(function_data, list):
        """if we have a list of function names, create cascaded Function objects
        for all function names in the list::

            ['FUNC2', 'PKG1.FUNC1'] --> FUNC2(PKG1.FUNC1(..))
        """
        function = None
        for name in reversed(function_data):
            packages = name.split('.')

            if function is None:
                """for the innermost function use the passed-in parameters as argument,
                otherwise use the prior created function
                """
                args = params
            else:
                args = [function]

            function = Function(packages.pop(-1),
                        packagenames=packages,
                        *args)

        return function

    elif isinstance(function_data, types.FunctionType):
        """if we have a function, call this function with the parameters and
        return the created Function object
        """
        if hasattr(element, 'flags'):
            # when element is a BaseFunction
            flags = element.flags
        else:
            flags = {}

        return function_data(params, within_column_clause, **flags)

    else:
        packages = function_data.split('.')

        return Function(packages.pop(-1),
                        packagenames=packages,
                        *params
                        )


class BaseFunction(Function):
    """Represents a database function.

    When the function is used on a molecule column additional arguments are set
    using __call__. The column or molecule the function is called on
    is stored inside the constructor (see ChemicalComparator.__getattr__()).
    When the function is called directly all arguments are set using the
    constructor.
    """

    def __init__(self, *arguments, **kwargs):
        self.arguments = arguments
        self.flags = kwargs.copy()

        Function.__init__(self, self.__class__.__name__, **kwargs)

    def __call__(self, *arguments, **kwargs):
        if len(arguments) > 0:
            self.arguments = self.arguments + arguments

        if len(kwargs) > 0:
            self.flags.update(kwargs)

        return self


@compiles(BaseFunction)
def __compile_base_function(element, compiler, **kw):

    params = [parse_clause(argument, compiler)
              for argument in element.arguments]

    from razi.dialect import DialectManager
    chemical_dialect = DialectManager.get_chemical_dialect(compiler.dialect)

    if chemical_dialect.is_member_function(element.__class__):
        chemical = params.pop(0)
        function_name = chemical_dialect.get_function(element.__class__)

        if isinstance(function_name, str):
            """If the function is defined as String (e.g.
            "oracle_functions.dims : 'Get_Dims'"), we construct the function
            call in the query ourselves. This is because SQLAlchemy,
            at this point of the compile process, does not add parenthesis for
            functions without arguments when using Oracle.
            Otherwise let SQLAlchemy build the query. Note that this won't work
            for Oracle with functions without parameters."""

            return "%s.%s(%s)" % (
                        compiler.process(chemical),
                        function_name,
                        ", ".join([compiler.process(e) for e in params])
                        )
        else:
            function = _get_function(element, compiler, params,
                                     kw.get('within_columns_clause', False))

            return "%s.%s" % (
                compiler.process(chemical), compiler.process(function)
            )

    elif chemical_dialect.is_property(element.__class__):
        chemical = params.pop(0)
        function_name = chemical_dialect.get_function(element.__class__)

        return "%s.%s" % (compiler.process(chemical), function_name)

    else:
        function = _get_function(element, compiler, params,
                                 kw.get('within_columns_clause', False))
        return compiler.process(function)


class DialectSpecificFunction(BaseFunction):
    pass


@compiles(DialectSpecificFunction)
def __compile_dialect_specific_function(element, compiler, **kw):
    from razi.dialect import DialectManager
    chemical_dialect = DialectManager.get_chemical_dialect(compiler.dialect)
    function = chemical_dialect.get_function(element.__class__)
    args = list(element.arguments)
    return compiler.process(function(compiler, element, *args))


class functions:
    """Functions that are supported by most databases
    """

    class smiles(BaseFunction):
        """Translates a :class:`~razi.chemtypes.Molecule` into a SMILES string representation
        """
        pass

    class mw(BaseFunction):
        """Computes the molecular weight for a :class:`~razi.chemtypes.Molecule`"""
        pass

    class logp(BaseFunction):
        """Estimates the LogP for a :class:`~razi.chemtypes.Molecule`"""
        pass

    class tpsa(BaseFunction):
        """Computes the topological polar surface area for a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class hba(BaseFunction):
        """Returns the number of hydrogen-bond acceptors for a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class hbd(BaseFunction):
        """Returns the number of hydrogen-bond donors for a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class num_atoms(BaseFunction):
        """Returns the total number of atoms in a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class num_hetatoms(BaseFunction):
        """Returns the number of hetero atoms in a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class num_hvy_atoms(BaseFunction):
        """Returns the number of heavy atoms in a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class num_rings(BaseFunction):
        """Returns the number of rings in a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class num_rotatable_bonds(BaseFunction):
        """Returns the number of rotatable bonds in a :class:`~razi.chemtypes.Molecule`
        """
        pass

    class equals(DialectSpecificFunction):
        """Returns true if two :class:`~razi.chemtypes.Molecule` objects represent the same structure.

        """
        pass

    class contains(DialectSpecificFunction):
        """Returns true if the second argument represents a substructure of the :class:`~razi.chemtypes.Molecule` object. If the second argument is a string, SMILES semantic is assumed (the string is coerced to :class:`~razi.chemtypes.Molecule` type).
        """
        pass

    class contained_in(DialectSpecificFunction):
        """Superstructure test. Returns true if the first argument represents a fragment of the :class:`~razi.chemtypes.Molecule` object passed as second argument.
        """
        pass

    class match(DialectSpecificFunction):
        """Returns true if the second argument represents a substructure of the :class:`~razi.chemtypes.Molecule` object. If the second argument is a string, SMARTS semantic is assumed (the string is coerced to :class:`~razi.chemtypes.QMolecule` type).
        """
        pass

    class morgan_b(BaseFunction):
        """Returns a :class:`~razi.chemtypes.BitFingerprint` which is the bit vector Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
        """
        pass

    class morgan_feat_b(BaseFunction):
        """Returns a :class:`~razi.chemtypes.BitFingerprint` which is the bit vector Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
        """
        pass

    class atompair_b(BaseFunction):
        """Returns a :class:`~razi.chemtypes.BitFingerprint` which is the bit vector atom-pair fingerprint for a molecule.
        """
        pass

    class torsion_b(BaseFunction):
        """Returns a :class:`~razi.chemtypes.BitFingerprint` which is the bit vector topological-torsion fingerprint for a molecule.
        """
        pass

    class layered_b(BaseFunction):
        """Returns a :class:`~razi.chemtypes.BitFingerprint` which is the layered fingerprint for a molecule. This is an experimental substructure fingerprint using hashed molecular subgraphs.
        """
        pass

    class morgan_c(BaseFunction):
        """Returns a :class:`~razi.chemtypes.CntFingerprint` which is the count-based Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
        """
        pass

    class morgan_feat_c(BaseFunction):
        """Returns a :class:`~razi.chemtypes.CntFingerprint` which is the count-based Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
        """
        pass

    class atompair_c(BaseFunction):
        """Returns a :class:`~razi.chemtypes.CntFingerprint` which is the count-based atom-pair fingerprint for a molecule.
        """
        pass

    class torsion_c(BaseFunction):
        """Returns a :class:`~razi.chemtypes.CntFingerprint` which is the count-based topological-torsion fingerprint for a molecule.
        """
        pass

    class tanimoto_similarity(BaseFunction):
        """Returns the Tanimoto similarity between two fingerprints of the same type (either two :class:`~razi.chemtypes.CntFingerprint` or two :class:`~razi.chemtypes.BitFingerprint` values).
        """
        pass

    class tanimoto_similar(DialectSpecificFunction):
        """Returns true if the Tanimoto similarity between two fingerprints is higher than the currently set threshold value.
        """
        pass

    class dice_similarity(BaseFunction):
        """Returns the Dice similarity between two fingerprints of the same type (either two :class:`~razi.chemtypes.CntFingerprint` or two :class:`~razi.chemtypes.BitFingerprint` values).
        """
        pass

    class dice_similar(DialectSpecificFunction):
        """Returns true if the Dice similarity between two fingerprints is higher than the currently set threshold value.
        """
        pass

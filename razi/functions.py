from sqlalchemy import types as sqltypes
from sqlalchemy.sql import functions
from sqlalchemy.dialects import postgresql

from . import types

# GenericFunction subclasses strongly based on the geoalchemy2
# implementation

_FUNCTIONS = [
    #
    # mol constructors
    #

    ('mol_from_smiles', types.Mol, None,),
    ('mol_from_smarts', types.Mol, None,),
    ('mol_from_ctab', types.Mol, None,),
    ('mol_from_pkl', types.Mol, None,),

    ('mol_adjust_query_properties', types.Mol, None,),

    #
    # mol conversion functions
    #

    ('mol_to_smiles', sqltypes.String, None,),
    ('mol_to_smarts', sqltypes.String, None,),
    ('mol_to_ctab', sqltypes.String, None,),
    ('mol_to_svg', sqltypes.String, None,),
    ('mol_to_pkl', postgresql.BYTEA, None,),

    #
    # mol descriptors
    #

    ('mol_amw', sqltypes.Float, None,),
    ('mol_logp', sqltypes.Float, None,),
    ('mol_fractioncsp3', sqltypes.Float, None,),
    ('mol_tpsa', sqltypes.Float, None,),
    ('mol_chi0n', sqltypes.Float, None,),
    ('mol_chi1n', sqltypes.Float, None,),
    ('mol_chi2n', sqltypes.Float, None,),
    ('mol_chi3n', sqltypes.Float, None,),
    ('mol_chi4n', sqltypes.Float, None,),
    ('mol_chi0v', sqltypes.Float, None,),
    ('mol_chi1v', sqltypes.Float, None,),
    ('mol_chi2v', sqltypes.Float, None,),
    ('mol_chi3v', sqltypes.Float, None,),
    ('mol_chi4v', sqltypes.Float, None,),
    ('mol_kappa1', sqltypes.Float, None,),
    ('mol_kappa2', sqltypes.Float, None,),
    ('mol_kappa3', sqltypes.Float, None,),
]

# Iterate through _FUNCTION and create GenericFunction classes dynamically
for name, type_, doc in _FUNCTIONS:
    attributes = {'name': name}
    docs = []

    if doc is not None:
        docs.append(doc)

    if type_ is not None:
        attributes['type'] = type_

        type_str = '{0}.{1}'.format(type_.__module__, type_.__name__)
        docs.append('Return type: :class:`{0}`.'.format(type_str))

    if len(docs) != 0:
        attributes['__doc__'] = '\n\n'.join(docs)

    globals()[name] = type(name, (functions.GenericFunction,), attributes)

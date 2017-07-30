from sqlalchemy import String
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

    ('mol_to_smiles', String, None,),
    ('mol_to_smarts', String, None,),
    ('mol_to_ctab', String, None,),
    ('mol_to_svg', String, None,),
    ('mol_to_pkl', postgresql.BYTEA, None,),
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

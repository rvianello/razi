from sqlalchemy.sql import functions
from sqlalchemy.dialects import postgresql

from . import types

class mol_from_pkl(functions.GenericFunction):
    type = types.Mol

class mol_to_pkl(functions.GenericFunction):
    type = postgresql.BYTEA

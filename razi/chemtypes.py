# -*- coding: utf-8 -*-
from sqlalchemy.types import UserDefinedType
from sqlalchemy.ext.compiler import compiles

from razi.dialect import DialectManager
from razi.expression import MoleculeElement

class Molecule(UserDefinedType):
    """Molecule column type for chemical databases.
    
    Converts bind/result values to/from a generic Persistent value.
    This is used as a base class and overridden into dialect specific
    Persistent values.
    """
    
    name = 'MOLECULE'
    
    def __init__(self, chemical_index=True, **kwargs):
        self.chemical_index = chemical_index
        super(Molecule, self).__init__(**kwargs)
    
    def bind_processor(self, dialect):
        def process(value):
            if value is None:
                return value
            elif not isinstance(value, MoleculeElement):
                return value
            elif not isinstance(value.desc, MoleculeElement):
                return value.desc
            else:
                return value.desc.desc
        return process
        
    def result_processor(self, dialect, coltype=None):
        chemical_dialect = DialectManager.get_chemical_dialect(dialect)
        def process(value):
            if value is not None:
                return chemical_dialect.process_result(value, self)
            else:
                return value
        return process
    
    
@compiles(Molecule)
def _compile_molecule(type_, compiler, **kwargs):
    chemical_dialect = DialectManager.get_chemical_dialect(compiler.dialect)
    return chemical_dialect.db_column_type(type_)

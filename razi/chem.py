from sqlalchemy import Column, Table
from sqlalchemy.orm import column_property
from sqlalchemy.orm.interfaces import AttributeExtension
from sqlalchemy.sql.expression import Alias
from sqlalchemy.sql import expression
from sqlalchemy.ext.compiler import compiles

from razi.base import MoleculeBase, _to_mol, ChemicalComparator
from razi.dialect import DialectManager
from razi.functions import functions

class Molecule(MoleculeBase):
    """Molecule column type.
    
    Converts bind/result values to/from a dialect specific persistent
    molecule value.
    
    """
    
    def result_processor(self, dialect, coltype=None):
        chemical_dialect = DialectManager.get_chemical_dialect(dialect)
        def process(value):
            if value is not None:
                return chemical_dialect.process_result(value, self)
            else:
                return value
        return process

class ChemistryDDL(object):
    """A DDL extension which integrates SQLAlchemy table create/drop 
    methods with column management functions of chemical databases.
    
    Usage::
    
        sometable = Table('sometable', metadata, ...)
        
        ChemistryDDL(sometable)

        sometable.create()
    
    """
    
    def __init__(self, table):
        for event in ('before-create', 'after-create', 
                      'before-drop', 'after-drop'):
            table.ddl_listeners[event].append(self)
        self._stack = []
        
    def __call__(self, event, table, bind):
        chemical_dialect = DialectManager.get_chemical_dialect(bind.dialect)
        if event in ('before-create', 'before-drop'):
            """Remove chemistry column from column list (table._columns), 
            so that it does not show up in the create statement 
            ("create table tab (..)"). Afterwards (on event 'after-create') 
            restore the column list from self._stack.
            """
            regular_cols = [c for c in table.c 
                            if not isinstance(c.type, Molecule)]
            chem_cols = set(table.c).difference(regular_cols)
            self._stack.append(table.c)
            table._columns = expression.ColumnCollection(*regular_cols)
            
            if event == 'before-drop':
                for c in chem_cols:
                    chemical_dialect.handle_ddl_before_drop(bind, table, c)
                
        elif event == 'after-create':
            table._columns = self._stack.pop()
            
            for c in table.c:
                if isinstance(c.type, Molecule):
                    chemical_dialect.handle_ddl_after_create(bind, table, c)

        elif event == 'after-drop':
            table._columns = self._stack.pop()


class ChemicalAttribute(AttributeExtension):
    """Intercepts 'set' events on a mapped instance attribute and 
    converts the incoming value to a molecule expression.
    
    """
    
    def set(self, state, value, oldvalue, initiator):
        return _to_mol(value)
 
 
class ChemExtensionColumn(Column):
    pass
        
def ChemColumn(*args, **kw):
    """Define a declarative chemistry-aware column property.
    
    This just produces orm.column_property() with the appropriate
    extension and comparator_factory arguments.  The given arguments
    are passed through to Column.  The declarative module extracts
    the Column for inclusion in the mapped table.
    
    This method can also be used for non-declarative mappings to 
    set the properties for a chemistry column when defining the mapping.
    
    """
    if kw.has_key("comparator"):
        comparator = kw.pop("comparator")
    else:
        comparator = ChemComparator
    
    if isinstance(args[0], ChemExtensionColumn):
        # if used for non-declarative, use the column of the table definition
        column = args[0]
        args = args[1:]
    else:
        # if used for declarative, create a new column
        column = ChemExtensionColumn(*args, **kw) 
    
    return column_property(
        column, 
        extension=ChemicalAttribute(), 
        comparator_factory=comparator
    )


from sqlalchemy.orm.properties import ColumnProperty
from sqlalchemy.sql import expression
from sqlalchemy.sql.expression import ColumnClause, literal
from sqlalchemy.types import TypeEngine
from sqlalchemy.ext.compiler import compiles

from functions import functions, _get_function, BaseFunction

class MoleculeElement(object):
    """Represents a molecular structure value."""

    def __str__(self):
        return self.desc

    def __repr__(self):
        return "<%s at 0x%x; %r>" % (self.__class__.__name__, 
                                     id(self), self.desc)
    
    def __getattr__(self, name):
        return getattr(functions, name)(self)


class TxtMoleculeElement(MoleculeElement, expression.Function):
    """Represents a Molecule value expressed within application code (a SMILES).
    
    Extends expression.Function so that in a SQL expression context the value 
    is interpreted as 'txt_to_mol(value)' or as the equivalent function in 
    the currently used database.
    """
    
    def __init__(self, desc, chemical_type='MOLECULE'):
        assert isinstance(desc, basestring)
        self.desc = desc
        self.chemical_type = chemical_type
        expression.Function.__init__(self, "")


@compiles(TxtMoleculeElement)
def __compile_txtmoleculeelement(element, compiler, **kw):
    function = _get_function(element, compiler, [element.desc], 
                             kw.get('within_columns_clause', False))
    
    return compiler.process(function)


class PersistentMoleculeElement(MoleculeElement):
    """Represents a Molecule value loaded from the database."""
    
    def __init__(self, desc):
        self.desc = desc

    
class MoleculeBase(TypeEngine):
    """Base Molecule column type for all spatial databases.
    
    Converts bind/result values to/from a generic Persistent value.
    This is used as a base class and overridden into dialect specific
    Persistent values.
    """
    
    name = 'MOLECULE'
    
    def __init__(self, chemical_index=True, **kwargs):
        self.chemical_index = chemical_index
        super(MoleculeBase, self).__init__(**kwargs)
    
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
        def process(value):
            if value is not None:
                return PersistentMoleculeElement(value)
            else:
                return value
        return process
    

# ORM integration

def _to_mol(value):
    """Interpret a value as a Molecule-compatible construct."""

    if hasattr(value, '__clause_element__'):
        return value.__clause_element__()
    elif isinstance(value, (expression.ClauseElement, MoleculeElement)):
        return value
    elif isinstance(value, basestring):
        return TxtMoleculeElement(value)
    elif value is None:
        return None
    else:
        raise Exception("Invalid type")
    
class ChemicalComparator(ColumnProperty.ColumnComparator):
    """Intercepts standard Column operators on mapped class attributes
        and overrides their behavior.
    """
    
    def __getattr__(self, name):
        return getattr(functions, name)(self)
        
    # override the __eq__() operator (allows to use '==' on molecules)
    def __eq__(self, other): 
        return functions.equals(self, other)
    

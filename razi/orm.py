from sqlalchemy import Column
from sqlalchemy.orm import column_property
from sqlalchemy.orm.interfaces import AttributeExtension
from sqlalchemy.orm.properties import ColumnProperty
#from sqlalchemy.sql import expression

from razi.expression import _to_mol
#from razi.dialect import DialectManager
from razi.functions import functions

class ChemicalAttribute(AttributeExtension):
    """Intercepts 'set' events on a mapped instance attribute and 
    converts the incoming value to a molecule expression.
    
    **** DEPRECATED ****
    should be replaced by ORM Attribute event listeners
    """
    def set(self, state, value, oldvalue, initiator):
        return _to_mol(value)
 

class ChemComparator(ColumnProperty.ColumnComparator):
    """Intercepts standard Column operators on mapped class attributes
        and overrides their behavior.
    """
    
    def __getattr__(self, name):
        return getattr(functions, name)(self)
        
    # override the __eq__() operator (allows to use '==' on molecules)
    def __eq__(self, other): 
        return functions.equals(self, other)

 
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
        extension=ChemicalAttribute(), #    **** DEPRECATED ****
        comparator_factory=comparator
    )
    


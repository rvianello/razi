from sqlalchemy.sql import expression
from sqlalchemy.ext.compiler import compiles

from razi.functions import functions, _get_function

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


    

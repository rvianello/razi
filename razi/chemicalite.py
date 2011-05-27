from sqlalchemy import select, func

from razi.chem import ChemComparator
from razi.molecule import Molecule, \
    PersistentMoleculeElement, TxtMoleculeElement
from razi.dialect import ChemicalDialect 
from razi.functions import functions, BaseFunction

def _configure_connection(connection):
    # FIXME - name of extension modules shouldn't be hardcoded
    signtree = 'libsigntree.so'
    chemicalite = 'libchemicalite.so'
    connection.enable_load_extension(True)
    connection.load_extension(signtree)
    connection.load_extension(chemicalite)


class ChemicaLiteComparator(ChemComparator):
    """Comparator class used for ChemicaLite
    """
    def __getattr__(self, name):
        try:
            return ChemComparator.__getattr__(self, name)
        except AttributeError:
            return getattr(chemicalite_functions, name)(self)


class ChemicaLitePersistentSpatialElement(PersistentMoleculeElement):
    """Represents a Molecule value as loaded from the database."""
    
    def __init__(self, desc):
        self.desc = desc
        
    def __getattr__(self, name):
        try:
            return PersistentMoleculeElement.__getattr__(self, name)
        except AttributeError:
            return getattr(chemicalite_functions, name)(self)


class chemicalite_functions(functions):
    """Functions only supported by ChemicaLite
    """
    
    #class func(BaseFunction):
    #    """Func(m)"""
    #    pass
    
    pass



class ChemicaLiteDialect(ChemicalDialect):
    """Implementation of ChemicalDialect for ChemicaLite."""
    
    __functions = { 
        TxtMoleculeElement: 'mol',
        #functions.func1: 'Func1',
        #chemicalite_functions.func2 : 'Func2',
        }
    
    def _get_function_mapping(self):
        return ChemicaLiteDialect.__functions
    
    def process_result(self, value, type):
        return ChemicaLitePersistentSpatialElement(value)
    
    def db_column_type(self, type_):
        if isinstance(type_, Molecule):
            return 'MOLECULE'
        raise NotImplementedError("DB column for type %s not supported by the "
                                  "current chemical dialect " % type(type_))

    def _ddl_after_create(self, table, column, bind):
        """This method is called on chemical columns after the mapped table 
        was created in the database by SQLAlchemy.
        """
        if column.type.chemical_index:
            # call chemicalite function that creates the index virtual table 
            # and triggers
            q = select([func.mol_structural_index(table.name, column.name)])
            bind.execute(q)
    
    def _ddl_before_drop(self, table, column, bind):
        """This method is called on chemical columns before the mapped table 
        is dropped from the database by SQLAlchemy.
        """
        import pdb; pdb.set_trace()
        if column.type.chemical_index:
            # call chemicalite function that drops the triggers
            # bind.execute(select([func.DisableChemicalIndex(table.name, column.name)]).execution_options(autocommit=True))
            # drop the index virtual table
            # bind.execute("DROP TABLE idx_%s_%s" % (table.name, column.name));
            pass


    

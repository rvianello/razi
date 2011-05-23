from sqlalchemy import select, func

from razi.base import ChemicalComparator, PersistentMoleculeElement, \
    TxtMoleculeElement
from razi.dialect import ChemicalDialect 
from razi.functions import functions, BaseFunction

from geoalchemy.geometry import GeometryExtensionColumn


class ChemicaLiteComparator(ChemicalComparator):
    """Comparator class used for ChemicaLite
    """
    def __getattr__(self, name):
        try:
            return ChemicalComparator.__getattr__(self, name)
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
        #functions.func1: 'Func1',
        #chemicalite_functions.func2 : 'Func2',
        }
    
    def _get_function_mapping(self):
        return ChemicaLiteDialect.__functions
    
    def process_result(self, value, type):
        return ChemicaLitePersistentSpatialElement(TxtMoleculeElement(value))
    
    def handle_ddl_before_drop(self, bind, table, column):
        if column.type.chemical_index:
            # call chemicalite function that drops the triggers
            # bind.execute(select([func.DisableChemicalIndex(table.name, column.name)]).execution_options(autocommit=True))
            # drop the index virtual table
            # bind.execute("DROP TABLE idx_%s_%s" % (table.name, column.name));
    
    def handle_ddl_after_create(self, bind, table, column):
        if column.type.spatial_index:
            # call chemicalite function that creates the index virtual table and triggers
            # bind.execute("SELECT CreateChemicalIndex('%s', '%s')" % (table.name, column.name))
            # bind.execute("VACUUM %s" % table.name) ?
    

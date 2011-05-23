# -*- coding: utf-8 -*-
from sqlalchemy import select, func, and_

from razi.base import ChemicalComparator, PersistentMoleculeElement, \
    TxtMoleculeElement
from razi.dialect import ChemicalDialect 
from razi.functions import functions, BaseFunction

class PostgresRDKitComparator(ChemicalComparator):
    """Comparator class used for PostgreSQL+RDKit
    """
    def __getattr__(self, name):
        try:
            return ChemicalComparator.__getattr__(self, name)
        except AttributeError:
            return getattr(pgrdkit_functions, name)(self)


class PostgresRDKitPersistentMoleculeElement(PersistentMoleculeElement):
    """Represents a Molecule value as loaded from the database."""

    def __init__(self, desc):
        self.desc = desc
        
    def __getattr__(self, name):
        try:
            return PersistentMoleculeElement.__getattr__(self, name)
        except AttributeError:
            return getattr(pgrdkit_functions, name)(self)


class pgrdkit_functions(functions):
    """Functions only supported by PostgreSQL+RDKit
    """
    #class func(BaseFunction):
    #    """Func(m)"""
    #    pass
    
    pass

class PostgresRDKitDialect(ChemicalDialect):
    """Implementation of ChemicalDialect for PostgreSQL+RDKit."""
    
    __functions = {
                   #functions.func1: 'Func1',
                   #pgrdkit_functions.func2 : 'Func2',
                  }
    
    def _get_function_mapping(self):
        return PostgresRDKitDialect.__functions
    
    def process_result(self, value, type):
        return PostgresRDKitPersistentMoleculeElement(TxtMoleculeElement(value))
    
    def handle_ddl_before_drop(self, bind, table, column):
        # drop the chemical column
        pass
    
    def handle_ddl_after_create(self, bind, table, column):    
        # add the chemical column

        # eventually also add a structural index and not NULL constraint 
        if column.type.chemical_index:
            bind.execute("CREATE INDEX \"idx_%s_%s\" ON \"%s\".\"%s\" USING GIST (%s)" % 
                         (table.name, column.name, (table.schema or 'public'), table.name, column.name))
            
        if not column.nullable:
            bind.execute("ALTER TABLE \"%s\".\"%s\" ALTER COLUMN \"%s\" SET not null" % 
                            ((table.schema or 'public'), table.name, column.name))
            

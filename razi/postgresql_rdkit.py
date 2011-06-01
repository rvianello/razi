# -*- coding: utf-8 -*-
from sqlalchemy import select, func, and_

from razi.chem import ChemComparator
from razi.chemtypes import Molecule
from razi.expression import PersistentMoleculeElement, \
    TxtMoleculeElement, TxtPatternElement
from razi.dialect import ChemicalDialect 
from razi.functions import functions, BaseFunction

class PostgresRDKitComparator(ChemComparator):
    """Comparator class used for PostgreSQL+RDKit
    """
    def __getattr__(self, name):
        try:
            return ChemComparator.__getattr__(self, name)
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
        #TxtMoleculeElement: '',
        #TxtPatternElement: '',
        
        #functions.smiles: '',
        #functions.mw: '',
        #functions.logp: '',
        #functions.tpsa: '',
        #functions.hba: '',
        #functions.hbd: '',
        #functions.num_atoms: '',
        #functions.num_hetatoms: '',
        #functions.num_hvy_atoms: '',
        #functions.num_rings: '',
        #functions.num_rotatable_bonds: '',
        
        #functions.equals: pgrdkit_functions._equals,
        #functions.contains: pgrdkit_functions._contains,
        #functions.contained_in: pgrdkit_functions._contained_in, 
        #functions.matches: pgrdkit_functions._matches,                    #functions.func1: 'Func1',
        
        #pgrdkit_functions.func2 : 'Func2',
        }
    
    def _get_function_mapping(self):
        return PostgresRDKitDialect.__functions
    
    def process_result(self, value, type):
        return PostgresRDKitPersistentMoleculeElement(value)
    
    def db_column_type(self, type_):
        if isinstance(type_, Molecule):
            return 'mol'
        raise NotImplementedError("DB column for type %s not supported by the "
                                  "current chemical dialect " % type(type_))

    def _ddl_before_drop(self, table, column, bind):
        # nothing to do?
        pass
    
    def _ddl_after_create(self, table, column, bind):    
        # eventually add a structural index and not NULL constraint 
        if column.type.chemical_index:
            bind.execute("CREATE INDEX \"idx_%s_%s\" ON \"%s\".\"%s\" USING GIST (%s)" % 
                         (table.name, column.name, (table.schema or 'public'), 
                          table.name, column.name))
            

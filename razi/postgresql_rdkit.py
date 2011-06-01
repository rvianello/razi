# -*- coding: utf-8 -*-
#from sqlalchemy.sql import select, func, and_
from sqlalchemy.ext.compiler import compiles

from razi.chem import ChemComparator
from razi.chemtypes import Molecule
from razi.expression import PersistentMoleculeElement, \
    TxtMoleculeElement, TxtPatternElement
from razi.dialect import ChemicalDialect 
from razi.functions import functions, parse_clause, BaseFunction

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
    
    class _Cast(BaseFunction):
        def __init__(self, totype, *args, **kw):
            self.totype = totype
            BaseFunction.__init__(self, *args, **kw)
            
    @staticmethod
    def mol(params, within_column_clause, **flags):
        # TODO: check carefully 'within_column_clause'
        return pgrdkit_functions._Cast('mol', *params, **flags)

    @staticmethod
    def qmol(params, within_column_clause, **flags):
        # TODO: check carefully 'within_column_clause'
        return pgrdkit_functions._Cast('qmol', *params, **flags)


@compiles(pgrdkit_functions._Cast)
def __compile_pgrdkit_cast(element, compiler, **kw):
    chemical = parse_clause(element.arguments[0], compiler)
    return "%s::%s" % (compiler.process(chemical), element.totype)
         
        
class PostgresRDKitDialect(ChemicalDialect):
    """Implementation of ChemicalDialect for PostgreSQL+RDKit."""
    
    __functions = {
        TxtMoleculeElement: pgrdkit_functions.mol,
        TxtPatternElement: pgrdkit_functions.qmol,
        
        #functions.smiles: '',
        functions.mw: 'mol_amw',
        functions.logp: 'mol_logp',
        functions.tpsa: 'mol_tpsa',
        functions.hba: 'mol_hba',
        functions.hbd: 'mol_hbd',
        functions.num_atoms: 'mol_numatoms',
        functions.num_hetatoms: 'mol_numheteroatoms',
        functions.num_hvy_atoms: 'mol_numheavyatoms',
        functions.num_rings: 'mol_numrings',
        functions.num_rotatable_bonds: 'mol_numrotatablebonds',
        
        #functions.equals: pgrdkit_functions._equals,
        #functions.contains: pgrdkit_functions._contains,
        #functions.contained_in: pgrdkit_functions._contained_in, 
        #functions.matches: pgrdkit_functions._matches,                    
        
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
            

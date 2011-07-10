# -*- coding: utf-8 -*-
from sqlalchemy import Column
from sqlalchemy.sql import func, cast, text
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.sql.expression import ClauseElement, Executable

from razi.chem import ChemComparator
from razi.chemtypes import Molecule, QMolecule, BitFingerprint, CntFingerprint
from razi.expression import \
    PersistentMoleculeElement, PersistentBitFingerprintElement, \
    PersistentCntFingerprintElement, \
    TxtMoleculeElement, TxtQMoleculeElement
from razi.dialect import ChemicalDialect 
from razi.functions import functions, parse_clause #, BaseFunction

class GUC(Executable, ClauseElement):
    
    def __init__(self, variable, *args, **kwargs):
        self.variable = variable
        #Executable.__init__(self, self.__class__.__name__, **kwargs)
        
    def set(self, value):
        xpr = text('SET %s=:value' % self.variable)
        return xpr.execution_options(autocommit=True).params(value=value)
        
    def get(self):
        return text('SHOW %s' % self.variable)

@compiles(GUC)
def __compile_guc(element, compiler, **kwargs):
    return compiler.process(element.get())
    
tanimoto_threshold = GUC('rdkit.tanimoto_threshold')
dice_threshold = GUC('rdkit.dice_threshold')

class PostgresRDKitComparator(ChemComparator):
    """Comparator class used for PostgreSQL+RDKit
    """
    def __getattr__(self, name):
        try:
            return ChemComparator.__getattr__(self, name)
        except AttributeError:
            return getattr(pgrdkit_functions, name)(self)


class PostgresRDKitPersistentMoleculeElement(PersistentMoleculeElement):

    def __init__(self, desc):
        self.desc = desc
        
    def __getattr__(self, name):
        try:
            return PersistentMoleculeElement.__getattr__(self, name)
        except AttributeError:
            return getattr(pgrdkit_functions, name)(self)


class PostgresRDKitPersistentBitFingerprintElement(PersistentBitFingerprintElement):

    def __init__(self, desc):
        self.desc = desc
        
    def __getattr__(self, name):
        try:
            return PersistentBitFingerprintElement.__getattr__(self, name)
        except AttributeError:
            return getattr(pgrdkit_functions, name)(self)


class PostgresRDKitPersistentCntFingerprintElement(PersistentCntFingerprintElement):

    def __init__(self, desc):
        self.desc = desc
        
    def __getattr__(self, name):
        try:
            return PersistentCntFingerprintElement.__getattr__(self, name)
        except AttributeError:
            return getattr(pgrdkit_functions, name)(self)


class pgrdkit_functions(functions):
    """Functions only supported by PostgreSQL+RDKit
    """
    
    @staticmethod
    def _equals(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler, TxtMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtMoleculeElement)
        if isinstance(m1, Column) and\
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
            return m1.op('@=')(m2)
        elif isinstance(m2, Column) and\
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
            return m2.op('@=')(m1)
        else:
            return func.mol_eq(m1, m2)
                                      
    @staticmethod
    def _contains(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler, TxtMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtMoleculeElement)
        if isinstance(m1, Column) and\
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
            return m1.op('@>')(m2)
        elif isinstance(m2, Column) and\
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
            return m2.op('<@')(m1)
        else:
            return func.substruct(m1, m2)
                                      
    @staticmethod
    def _contained_in(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler, TxtMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtMoleculeElement)
        if isinstance(m1, Column) and\
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
            return m1.op('<@')(m2)
        elif isinstance(m2, Column) and\
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
            return m2.op('@>')(m1)
        else:
            return func.rsubstruct(m1, m2)
                                      
    @staticmethod
    def _match(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler, TxtQMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtQMoleculeElement)
        if isinstance(m1, Column) and\
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
            return m1.op('@>')(m2)
        elif isinstance(m2, Column) and\
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
            return m2.op('<@')(m1)
        else:
            return func.substruct(m1, m2)
                                      
    @staticmethod
    def _tanimoto(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler)
        m2 = parse_clause(arg2, compiler)
        if (isinstance(m1, Column) and \
            isinstance(m1.type, BitFingerprint) and \
            m1.type.chemical_index) or \
           (isinstance(m2, Column) and \
            isinstance(m2.type, BitFingerprint) and \
            m2.type.chemical_index):
            # apparently, we can just use the mod operator to generate the 
            # desired SQL expression
            return m1 % m2
        else:
            return func.tanimoto_sml_op(m1, m2)
                                      
    @staticmethod
    def _dice(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler)
        m2 = parse_clause(arg2, compiler)
        if (isinstance(m1, Column) and \
            isinstance(m1.type, BitFingerprint) and \
            m1.type.chemical_index) or \
           (isinstance(m2, Column) and \
            isinstance(m2.type, BitFingerprint) and \
            m2.type.chemical_index):
            return m1.op('#')(m2)
        else:
            return func.dice_sml_op(m1, m2)
                                      
    @staticmethod
    def mol(params, within_column_clause, **flags):
        # TODO: check carefully 'within_column_clause'
        (param,) = params
        return cast(param, Molecule)
        
    @staticmethod
    def qmol(params, within_column_clause, **flags):
        # TODO: check carefully 'within_column_clause'
        (param,) = params
        return cast(param, QMolecule)


class PostgresRDKitDialect(ChemicalDialect):
    """Implementation of ChemicalDialect for PostgreSQL+RDKit."""
    
    __functions = {
        TxtMoleculeElement: pgrdkit_functions.mol,
        TxtQMoleculeElement: pgrdkit_functions.qmol,
        
        # molecular descriptors
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
        
        # molecule comparison ops
        functions.equals: pgrdkit_functions._equals,
        functions.contains: pgrdkit_functions._contains,
        functions.contained_in: pgrdkit_functions._contained_in, 
        functions.match: pgrdkit_functions._match,                    
        
        # bitstring fingerprint generation
        functions.morgan_b: 'morganbv_fp',
        functions.morgan_feat_b: 'featmorganbv_fp',
        functions.atompair_b: 'atompairbv_fp',
        functions.torsion_b: 'torsionbv_fp',
        functions.layered_b: 'layered_fp',
                
        # count vector fingerprint generation
        functions.morgan_c: 'morgan_fp',
        functions.morgan_feat_c: 'featmorgan_fp',
        functions.atompair_c: 'atompair_fp',
        functions.torsion_c: 'torsion_fp',
                
        # fingerprint similarity
        functions.tanimoto_similarity : 'tanimoto_sml',
        functions.tanimoto_similar: pgrdkit_functions._tanimoto,
        functions.dice_similarity: 'dice_sml',
        functions.dice_similar: pgrdkit_functions._dice,
        
        #pgrdkit_functions.func2 : 'Func2',
        }
    
    def _get_function_mapping(self):
        return PostgresRDKitDialect.__functions
    
    def process_result(self, value, type_):
        if isinstance(type_, Molecule):
            return PostgresRDKitPersistentMoleculeElement(value)
        elif isinstance(type_, BitFingerprint):
            return PostgresRDKitPersistentBitFingerprintElement(value)
        elif isinstance(type_, CntFingerprint):
            return PostgresRDKitPersistentCntFingerprintElement(value)
        raise NotImplementedError("DB column for type %s not supported by the "
                                  "current chemical dialect " % type(type_))
    
    def db_column_type(self, type_):
        if isinstance(type_, Molecule):
            return 'mol'
        elif isinstance(type_, QMolecule):
            return 'qmol'
        elif isinstance(type_, BitFingerprint):
            return 'bfp'
        elif isinstance(type_, CntFingerprint):
            return 'sfp'
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
            

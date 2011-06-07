from sqlalchemy import Column, select, func
from sqlalchemy.sql.expression import and_, table, column

from razi.chem import ChemComparator
from razi.chemtypes import Molecule
from razi.expression import PersistentMoleculeElement, TxtMoleculeElement
from razi.expression import TxtQMoleculeElement
from razi.dialect import ChemicalDialect 
from razi.functions import functions, parse_clause #, BaseFunction

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

    @staticmethod
    def _str_idx_constraint(m1, m2, op):
        idx_tbl = table("str_idx_%s_%s" % (m1.table.fullname, m1.key), 
                        column("s"), column("id"))
        return table(m1.table.fullname, column("rowid")).c.rowid.in_(
            select([idx_tbl.c.id])
                .where(idx_tbl.c.s.op(op)(func.mol_signature(m2)))
                )
                            
    @staticmethod
    def _equals(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler, TxtMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtMoleculeElement)
        q1 = func.mol_same(m1, m2)
        if isinstance(m1, Column) and \
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m1, m2, '='))
        elif isinstance(m2, Column) and \
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m2, m1, '='))
        else:
            return q1
    
    @staticmethod
    def _contains(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler, TxtMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtMoleculeElement)
        q1 = func.mol_is_substruct(m1, m2)
        if isinstance(m1, Column) and \
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m1, m2, '>='))
        elif isinstance(m2, Column) and \
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m2, m1, '<='))
        else:
            return q1
    
    @staticmethod
    def _contained_in(compiler, element, arg1, arg2): 
        m1 = parse_clause(arg1, compiler, TxtMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtMoleculeElement)
        q1 = func.mol_substruct_of(m1, m2)
        if isinstance(m1, Column) and \
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m1, m2, '<='))
        elif isinstance(m2, Column) and \
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m2, m1, '>='))
        else:
            return q1
    
    @staticmethod
    def _matches(compiler, element, arg1, arg2):
        m1 = parse_clause(arg1, compiler, TxtQMoleculeElement)
        m2 = parse_clause(arg2, compiler, TxtQMoleculeElement)
        q1 = func.mol_is_substruct(m1, m2)
        if isinstance(m1, Column) and \
            isinstance(m1.type, Molecule) and \
            m1.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m1, m2, '>='))
        elif isinstance(m2, Column) and \
            isinstance(m2.type, Molecule) and \
            m2.type.chemical_index:
                return and_(
                    q1,
                    chemicalite_functions._str_idx_constraint(m2, m1, '<='))
        else:
            return q1
    

class ChemicaLiteDialect(ChemicalDialect):
    """Implementation of ChemicalDialect for ChemicaLite."""
    
    __functions = { 
        TxtMoleculeElement: 'mol',
        TxtQMoleculeElement: 'qmol',
        
        functions.smiles: 'mol_smiles',
        functions.mw: 'mol_mw',
        functions.logp: 'mol_logp',
        functions.tpsa: 'mol_tpsa',
        functions.hba: 'mol_hba',
        functions.hbd: 'mol_hbd',
        functions.num_atoms: 'mol_num_atms',
        functions.num_hetatoms: 'mol_num_hetatms',
        functions.num_hvy_atoms: 'mol_num_hvyatms',
        functions.num_rings: 'mol_num_rings',
        functions.num_rotatable_bonds: 'mol_num_rotatable_bonds',
        
        functions.equals: chemicalite_functions._equals,
        functions.contains: chemicalite_functions._contains,
        functions.contained_in: chemicalite_functions._contained_in, 
        functions.matches: chemicalite_functions._matches, 
            
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


    

from sqlalchemy import schema as _schema
from sqlalchemy import event as _event
from sqlalchemy import pool as _pool
from sqlalchemy.dialects.postgresql.base import PGDialect
from sqlalchemy.dialects.sqlite.base import SQLiteDialect
#from sqlalchemy.dialects.mysql.base import MySQLDialect
#from sqlalchemy.dialects.oracle.base import OracleDialect
#from sqlalchemy.dialects.mssql.base import MSDialect

#from sqlalchemy import func

#from razi.functions import functions
#from razi.molecule import TxtMoleculeElement

class ChemicalDialect(object):
    """This class bundles all required classes and methods to support 
    a database dialect. It is supposed to be subclassed.
    The child classes must be added to DialectManager.__initialize_dialects(),
    so that they can be used.
    
    """
    
    __functions = {
                   #functions.func1: 'Func',
                  }
    
    def get_function(self, function_class):
        """This method is called to translate a generic function into a 
        database dialect specific function.
        
        It either returns a string, a list or a method that returns a Function 
        object (see also functions._get_function)
        
            - String: A function name, e.g.::
                    
                    functions.func: 'PACKAGE.FUNCTION'
            
            - List: A list of function names that are called cascaded, e.g.::
            
                    functions.func: 
                      ['TO_CHAR', 'PACKAGE.FUNCTION'] is compiled as:
                      TO_CHAR(PACKAGE.FUNCTION(..))
                    
            - Method: A method that accepts a list/set of arguments and returns
                      a Function object, e.g.::
            
                    functions.equals : lambda params, within_column_clause : (func.F_EQUAL(*params) == 'TRUE')
    
        """
        if self._get_function_mapping() is not None:
            if function_class in self._get_function_mapping():
                if self._get_function_mapping()[function_class] is None:
                    raise NotImplementedError("Operation '%s' is not supported for '%s'" 
                                                % (function_class.__name__, self.__class__.__name__)) 
                else:
                    return self._get_function_mapping()[function_class]
        
        return ChemicalDialect.__functions[function_class]
    
    def is_member_function(self, function_class):   
        """Returns True if the passed-in function should be called as member 
        function, e.g.::
        
            Point.the_geom.dims is compiled as:
            points.the_geom.Get_Dims()
            
        """
        return False
    
    def is_property(self, function_class):   
        """Returns True if the passed-in function should be called as property, 
        e.g.::
        
            Point.the_geom.x is compiled as (for MS Server):
            points.the_geom.x
            
        """
        return False
    
    def _get_function_mapping(self):
        """This method can be overridden in subclasses, to set a database 
        specific name for a function, or to add new functions, that are only 
        supported by the database.
        """
        return None
    
    def process_result(self, value, type):
        """This method is called when a molecule value from the database is
        transformed into a MoleculeElement object. 
        
        """
        raise NotImplementedError("Method ChemicalDialect.process_result must "
                                  "be implemented in subclasses.")
    
    def db_column_type(self, type_):
        raise NotImplementedError("Method ChemicalDialect.db_column_type must "
                                  "be implemented in subclasses.")
        
    def handle_ddl_after_create(self, table, connection):
        """This method is called after the mapped table was created in the 
        database by SQLAlchemy.
        """
        from razi.chemtypes import Molecule
        chemical_cols = [c for c in table.c if isinstance(c.type, Molecule)]
        for c in chemical_cols:
            self._ddl_after_create(table, c, connection.engine)
            
    def _ddl_after_create(self, table, column, bind):
        """This method is called on chemical columns after the mapped table 
        was created in the database by SQLAlchemy.
        """
        raise NotImplementedError("Method "
                                  "ChemicalDialect._ddl_after_create must "
                                  "be implemented in subclasses.")
    
    def handle_ddl_before_drop(self, table, connection):
        """This method is called before the mapped table is dropped from the
        database by SQLAlchemy.
        """
        from razi.chemtypes import Molecule
        chemical_cols = [c for c in table.c if isinstance(c.type, Molecule)]
        for c in chemical_cols:
            self._ddl_before_drop(table, c, connection.engine)

    def _ddl_before_drop(self, table, column, bind):
        """This method is called on chemical columns before the mapped table 
        is dropped from the database by SQLAlchemy.
        """
        raise NotImplementedError("Method "
                                  "ChemicalDialect._ddl_before drop must "
                                  "be implemented in subclasses.")
    

class DialectManager(object):
    """This class is responsible for finding a chemical dialect (e.g. 
    PGRDKitDialect or ChemicaLiteDialect) for a SQLAlchemy database dialect.
    
    It can be used by calling "DialectManager.get_chemical_dialect(dialect)", 
    which returns the corresponding chemical dialect.
    The chemical dialect has to be listed in __initialize_dialects().
    
    """
    
    # all available chemical dialects 
    # {(SQLAlchemy dialect class: chemical dialect class)}    
    __dialects_mapping = None
    
    # all instantiated dialects 
    # {(chemical dialect class: chemical dialect instance)}    
    __chemical_dialect_instances = {}     
    
    @staticmethod
    def __initialize_dialects():
        #further spatial dialects can be added here
        from razi.postgresql_rdkit import PostgresRDKitDialect
        from razi.chemicalite import ChemicaLiteDialect
        #from razi.mysql import ?
        #from razi.oracle import ?
        #from razi.mssql import ?
            
        DialectManager.__dialects_mapping = {
                PGDialect: PostgresRDKitDialect,
                SQLiteDialect:  ChemicaLiteDialect,
                }

    @staticmethod
    def __dialects():
        if DialectManager.__dialects_mapping is None:
            DialectManager.__initialize_dialects()
            
        return DialectManager.__dialects_mapping
    
    @staticmethod
    def get_chemical_dialect(dialect):
        """This method returns a chemical dialect instance for a given 
        SQLAlchemy dialect.
        The instances are cached, so that for every chemical dialect exists 
        only one instance.
        
        """
        possible_chemical_dialects = [chem_dialect 
                                      for (main_dialect, chem_dialect) 
                                      in DialectManager.__dialects().items() 
                                      if isinstance(dialect, main_dialect)]
        
        if possible_chemical_dialects:
            # take the first possible spatial dialect
            chem_dialect = possible_chemical_dialects[0]
            if chem_dialect not in DialectManager.__chemical_dialect_instances:
                # if there is no instance for the given dialect yet, create one
                chem_dialect_instance = chem_dialect()
                DialectManager.__chemical_dialect_instances[chem_dialect] \
                    = chem_dialect_instance
                
            return DialectManager.__chemical_dialect_instances[chem_dialect]
        else:
            raise NotImplementedError('Dialect "%s" is not supported by '
                                      'Razi' % (dialect.name))
        
        
# custom configuration of db connections

@_event.listens_for(_pool.Pool, "connect")
def _configure_connection(dbapi_conn, conn_record):
    module_name = dbapi_conn.__class__.__module__.split('.')[0].lower()
    if module_name.find('sqlite') >= 0: 
        # sqlite3 or pysqlite2
        from razi import chemicalite
        chemicalite._configure_connection(dbapi_conn)
       
# custom ddl code

@_event.listens_for(_schema.Table, "after_create")
def _handle_ddl_after_create(table, connection, **kwargs):
    chemical_dialect = DialectManager.get_chemical_dialect(connection.dialect)
    chemical_dialect.handle_ddl_after_create(table, connection)
    
@_event.listens_for(_schema.Table, "before_drop")
def _handle_ddl_before_drop(table, connection, **kwargs):
    chemical_dialect = DialectManager.get_chemical_dialect(connection.dialect)
    chemical_dialect.handle_ddl_before_drop(table, connection)
    

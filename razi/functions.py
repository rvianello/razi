from sqlalchemy.sql.expression import Function, ClauseElement
from sqlalchemy.ext.compiler import compiles
from sqlalchemy import literal
from sqlalchemy.types import NullType, TypeDecorator

import types

def parse_clause(clause, compiler):
    """This method is used to translate a clause element (molecules, 
    functions, ..).
    According to the type of the clause, a conversion to the database 
    molecule type is added or the column clause (column name) or the cascaded 
    clause element is returned.
        
    """
    from razi.base import MoleculeElement, TxtMoleculeElement
    
    if hasattr(clause, '__clause_element__'):
        # for example a column name
        return clause.__clause_element__()
    elif isinstance(clause, ClauseElement):
        # for cascaded clause elements, like other functions
        return clause
    elif isinstance(clause, MoleculeElement):
        if isinstance(clause, TxtMoleculeElement):
            return clause
        return clause.desc
    elif isinstance(clause, basestring):
        return TxtMoleculeElement(clause)
    
    # for raw parameters    
    return literal(clause)


def _get_function(element, compiler, params, within_column_clause):
    """For elements of type BaseFunction, the database specific function data 
    is looked up and a executable sqlalchemy.sql.expression.Function object 
    is returned.
    """
    from razi.dialect import DialectManager 
    chemical_dialect = DialectManager.get_chemical_dialect(compiler.dialect)
    function_data = chemical_dialect.get_function(element.__class__)
    
    if isinstance(function_data, list):
        """if we have a list of function names, create cascaded Function objects
        for all function names in the list::
        
            ['FUNC2', 'PKG1.FUNC1'] --> FUNC2(PKG1.FUNC1(..))
        """
        function = None
        for name in reversed(function_data):
            packages = name.split('.')
            
            if function is None:
                """for the innermost function use the passed-in parameters as argument,
                otherwise use the prior created function
                """
                args = params
            else:
                args = [function]
                
            function = Function(packages.pop(-1), 
                        packagenames=packages,
                        *args
                        )
        
        return function
    
    elif isinstance(function_data, types.FunctionType):
        """if we have a function, call this function with the parameters and 
        return the created Function object
        """
        if hasattr(element, 'flags'):
            # when element is a BaseFunction
            flags = element.flags
        else:
            flags = {}
            
        return function_data(params, within_column_clause, **flags)
    
    else:
        packages = function_data.split('.')
        
        return Function(packages.pop(-1), 
                        packagenames=packages, 
                        *params
                        )

    
class BaseFunction(Function):
    """Represents a database function.
    
    When the function is used on a molecule column additional arguments are set
    using __call__. The column or molecule the function is called on 
    is stored inside the constructor (see ChemicalComparator.__getattr__()).
    When the function is called directly all arguments are set using the 
    constructor.
    """
    
    def __init__(self, *arguments, **kwargs):
        self.arguments = arguments
        self.flags = kwargs.copy()
        
        Function.__init__(self, self.__class__.__name__, **kwargs)
        
    def __call__(self, *arguments, **kwargs):
        if len(arguments) > 0:
            self.arguments =  self.arguments + arguments
        
        if len(kwargs) > 0:
            self.flags.update(kwargs)
        
        return self

@compiles(BaseFunction)
def __compile_base_function(element, compiler, **kw):
    
    params = [parse_clause(argument, compiler) 
              for argument in element.arguments]
    
    from razi.dialect import DialectManager 
    chemical_dialect = DialectManager.get_chemical_dialect(compiler.dialect)
    
    if chemical_dialect.is_member_function(element.__class__):
        chemical = params.pop(0)
        function_name = chemical_dialect.get_function(element.__class__)
        
        if isinstance(function_name, str):
            """If the function is defined as String (e.g. 
            "oracle_functions.dims : 'Get_Dims'"), we construct the function 
            call in the query ourselves. This is because SQLAlchemy, 
            at this point of the compile process, does not add parenthesis for 
            functions without arguments when using Oracle.
            Otherwise let SQLAlchemy build the query. Note that this won't work
            for Oracle with functions without parameters."""
            
            return "%s.%s(%s)" % (
                        compiler.process(chemical),
                        function_name,
                        ", ".join([compiler.process(e) for e in params]) 
                        )
        else:
            function = _get_function(element, compiler, params,
                                     kw.get('within_columns_clause', False))
            
            return "%s.%s" % (
                compiler.process(chemical), compiler.process(function)      
                )
            
    elif database_dialect.is_property(element.__class__):
        chemical = params.pop(0)
        function_name = chemical_dialect.get_function(element.__class__)
        
        return "%s.%s" % (compiler.process(chemical), function_name)
        
    else:
        function = _get_function(element, compiler, params, 
                                 kw.get('within_columns_clause', False))
        return compiler.process(function)

def BooleanFunction(function, compare_value = 'TRUE'):
    """Wrapper for database function that return 'Boolean-like' values.
    
    This function adds the necessary comparison to ensure
    that in the WHERE clause the function is compared to `compare_value`
    while in the SELECT clause the function result is returned.
    
    :param function: The function that needs boolean wrapping
    :param compare_value: The value that the function should be compare to
    """
    def function_handler(params, within_column_clause):
        return check_comparison(function(*params), within_column_clause, True, 
                                compare_value)
    
    return function_handler
  
def check_comparison(function, 
                     within_column_clause, returns_boolean, compare_value):
    """Because Oracle and MS SQL do not want to know Boolean and functions 
    return 0/1 or the string 'TRUE', we manually have to add a comparison, but 
    only if the function is called inside the where-clause, not in the select 
    clause.
    
    For example:
    select .. from .. where SDO_EQUAL(.., ..) = 'TRUE'
    select SDO_EQUAL(.., ..) from ..
    
    """
    if returns_boolean and not within_column_clause:
        return (function == compare_value)
    else: 
        return function

class functions:
    """Functions that are supported by most databases
    """
    
    #class func1(BaseFunction):
    #    """Func(m)"""
    #    pass


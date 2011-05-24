import sqlalchemy

def create_engine(*args, **kwargs):

    signtree = kwargs.pop('signtree', 'libsigntree.so')
    chemicalite = kwargs.pop('chemicalite', 'libchemicalite.so')
    
    engine = sqlalchemy.create_engine(*args, **kwargs)

    if engine.name == 'sqlite':
        connection = engine.raw_connection().connection
        connection.enable_load_extension(True)
        connection.load_extension(signtree)
        connection.load_extension(chemicalite)

    return engine

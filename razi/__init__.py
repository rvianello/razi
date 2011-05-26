import inspect

from sqlalchemy import pool
from sqlalchemy import event

import chemicalite

@event.listens_for(pool.Pool, "connect")
def _configure_connection(dbapi_conn, conn_record):
    module_name = dbapi_conn.__class__.__module__.split('.')[0].lower()
    if module_name.find('sqlite') >= 0: 
        # sqlite3 or pysqlite2
        chemicalite._configure_connection(dbapi_conn)

__all__ = sorted(name for name, obj in locals().items()
                 if not (name.startswith('_') or inspect.ismodule(obj)))

__version__ = '0.0.0'


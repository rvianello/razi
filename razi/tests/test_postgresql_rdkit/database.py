import os
try:
    url = os.environ['POSTGRESQL_RDKIT_TESTDB']
except KeyError:
    url = 'postgresql://db_user:db_password@hostname:1234/razi_test'


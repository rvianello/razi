#!/bin/bash
set -e

psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" <<-EOSQL
    CREATE DATABASE $RAZI_DB;
    CREATE USER "$RAZI_USER";
    ALTER ROLE "$RAZI_USER" WITH PASSWORD '$RAZI_PASSWORD';
    GRANT ALL PRIVILEGES ON DATABASE $RAZI_DB TO "$RAZI_USER";
EOSQL

# depending on the context (e.g. docker-compose up) pg_restore may be problematic to execute.
# it may be better to defer it to a later stage than cluster initialization. 
#pg_restore --username "$POSTGRES_USER" --dbname "$CHEMBL_DB" "/var/lib/postgresql/data/chembl_${CHEMBL_RELEASE}_postgresql.dmp"

# but also install the rdkit cartridge
psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$RAZI_DB" <<-EOSQL
    CREATE EXTENSION rdkit;
EOSQL



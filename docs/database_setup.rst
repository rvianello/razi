Database setup
==============

Depending on the database backend, some configuration steps may be required before a chemical database is created.

postgresql-rdkit
----------------

The RDKit extension for PostgreSQL must be installed. Please refer to the `RDKit wiki <http://code.google.com/p/rdkit/wiki/BuildingTheCartridge>`_ for detailed instructions.

creation of a chemical database template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Creating the template database.
    $ createdb -E UTF8 template_rdkit
    # Allows non-superusers the ability to create from this template
    $ psql -d postgres -c "UPDATE pg_database SET datistemplate='true' WHERE datname='template_rdkit';"
    # Loading the RDKit cartridge
    $ psql -d template_rdkit -f `pg_config --sharedir`/contrib/rdkit.sql

creation of the database
^^^^^^^^^^^^^^^^^^^^^^^^

::

    $ createdb -T template_rdkit my_chem_db


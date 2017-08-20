Creating a chemical database
============================

This tutorial is based on a similar document that is part of the `RDKit official documentation <http://www.rdkit.org/docs/Cartridge.html#creating-databases>`_ and it illustrates how to build a chemical database and perform some simple search queries using Razi. The Python bindings for the RDKit libraries will be used in some data pre-processing steps, so you'll need to have them available on your system.

Some basic understanding of SQLAlchemy is assumed. For additional details, please refer to the excellent `SQLALchemy ORM tutorial <http://www.sqlalchemy.org/docs/orm/tutorial.html>`_.

Download the ChEMBL compounds data
----------------------------------

Download the `ChEMBLdb <https://www.ebi.ac.uk/chembl/>`_ database and decompress it::

    $ wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_23/chembl_23_chemreps.txt.gz
    $ gunzip chembl_23_chemreps.txt.gz


Database creation
-----------------

Then, create a new database and configure it to use the RDKit extension::

    $ createdb -Udb_user razi_tutorial
    $ psql razi_tutorial db_user

    razi_tutorial=# create extension rdkit;
    CREATE EXTENSION
    razi_tutorial=# \q

*The database username and password, together with the name of the database itself, are most likely to differ in your case. Please replace them with the values appropriate to your work environment.*

Connection to the database
--------------------------

Start your python interpreter and configure a database connection::

    from sqlalchemy import create_engine
    engine = create_engine('postgresql://db_user:db_password@host:1234/razi_tutorial')

also, define the database session factory object::

    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)


Schema definition
-----------------

For this tutorial we only need one single database entity, mapped to a python class. Following the SQLAlchemy's *declarative* style, a base class is first defined::

    from sqlalchemy.ext.declarative import declarative_base
    Base = declarative_base(bind=engine)

then, the definition of the mapped entity follows::

    from sqlalchemy import Column, Index, Integer, String
    from razi.rdkit_postgresql.types import Mol

    class Compound(Base):
        __tablename__ = 'compounds'

        id = Column(Integer, primary_key=True)
        name = Column(String)
        structure = Column(Mol)

        __table_args__ = (
            Index('compounds_structure', 'structure',
                  postgresql_using='gist'),
            )

        def __init__(self, name, structure):
            self.name = name
            self.structure = structure

        def __repr__(self):
            return '(%s) < %s >' % (self.name, self.structure)

to actually create the database schema, a method of the ``Base`` class ``metadata`` attribute is called::

    Base.metadata.create_all()

In the present case this last command creates a table named ``compounds`` with columns ``id``, ``name`` and ``structure``, and it also includes the creation of a structural index on the column with type :class:`~razi.rdkit_postgresql.types.Mol`.

Inserting data
--------------

The data from ChEMBL is formatted as a csv file, tab-separated and beginning with a header line::

    $ head -n3  chembl_23_chemreps.txt
    chembl_id	canonical_smiles	standard_inchi	standard_inchi_key
    CHEMBL153534	Cc1cc(cn1C)c2csc(N=C(N)N)n2	InChI=1S/C10H13N5S/c1-6-3-7(4-15(6)2)8-5 <snip>
    CHEM1BL440060	CC[C@H](C)[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@@H](N)CCSC) <snip>

To semplify the processing of these data records we define a namedtuple matching the input format::

    from collections import namedtuple
    Record = namedtuple('Record', 'chembl_id, smiles, inchi, inchi_key')

During this tutorial only a portion of the whole database is actually imported, and at the same time we want to make sure that troublesome SMILES strings are skipped (SMILES containing errors, compounds that are too big and/or other strings that the RDKit cartridge won't be able to process). File parsing and data filtering can therefore be performed with a function similar to the following::

    import csv
    from rdkit import Chem

    def read_chembldb(filepath, limit=0):

        with open(filepath, 'rt') as inputfile:
            reader = csv.reader(inputfile, delimiter='\t', skipinitialspace=True)
            next(reader) # skip header

            for count, record in enumerate(map(Record._make, reader), 1):

                smiles = record.smiles

                # skip problematic smiles
                if len(smiles) > 300: continue
                smiles = smiles.replace('=N#N','=[N+]=[N-]')
                smiles = smiles.replace('N#N=','[N-]=[N+]=')
                if not Chem.MolFromSmiles(smiles):
                    continue

                yield count, record.chembl_id, smiles
                if count == limit:
                    break

The ``read_chembldb`` function above is a python generator, producing for each valid record a python tuple containing the record counter and the ``chembl_id`` and ``smiles`` strings.

With this function importing the compounds into the database reduces to a simple loop *(please note that depending on the available hardware resources importing the whole database may require a few hours; to keep this tutorial short we'll limit the processing to the first 25K compounds, a dataset size the usually corresponds to a few minutes)*::

    session = Session()
    for count, chembl_id, smiles in read_chembldb('chembl_08_chemreps.txt', 25000):
        compound = Compound(chembl_id, smiles)
	session.add(compound)
    session.commit()

Querying the database
---------------------

Finally, we can perform some queries. We can for example verify the number of compounds actually loaded into the database::

    >>> print session.query(Compound).count()
    24956

or display the first 5 compounds::

    >>> for compound in session.query(Compound)[:5]:
    ...     print compound
    ...
    (CHEMBL6582) < NC(=O)c1cc(-c2ccc(Cl)cc2)nnc1Cl >
    (CHEMBL6583) < Cc1cnc(NS(c2cccc3c(N(C)C)cccc23)(=O)=O)cn1 >
    (CHEMBL6584) < CN(C)/C=N/c1nc(/N=C\N(C)C)c2c(ncc(Sc3cc(Cl)c(Cl)cc3)n2)n1 >
    (CHEMBL6585) < CC12C(C[C@@H](I)[C@@H]1O)C1C(c3ccc(O)cc3CC1)CC2 >
    (CHEMBL6637) < C/C(=C\Cn1oc(=O)[nH]c1=O)c1ccc(OCCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)cc1 >


Finally (and hopefully more interestingly), here's a first example of a more chemistry-aware query, searching the database for a given substructure::

    >>> # which compounds contain 'c1cccc2c1nncc2' as a substructure?
    ...
    >>> subset = session.query(Compound)
    >>> subset = subset.filter(Compound.structure.contains('c1cccc2c1nncc2'))
    >>> print subset.count()
    2
    >>> for compound in subset: print compound
    ...
    (CHEMBL12112) < CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1 >
    (CHEMBL26025) < Cc1cccc(NC(=O)Nc2ccc3nnccc3c2)c1 >


Please notice how the SQLAlchemy's ORM API allows the incremental specification of the filtering clause (or clauses) associated to the main selection query and how the ``subset`` instance is actually used twice, in two distinct queries (to compute the number of record matching the clause and to retrieve the actual records). In addition to this, the returned records can also serve as the basis for further queries, also using the chemical functions provided by the database backend:

    >>> for compound in subset:
    ...     # run a query to compute the molecular weight for this compound
    ...     print session.scalar(compound.structure.mw)
    ...
    488.701
    278.315

Using Razi with the SQLAlchemy SQL Expression Language
======================================================

This tutorial is based on a similar document available from the `RDKit wiki <http://code.google.com/p/rdkit/wiki/DatabaseCreation>`_ and it illustrates how to build a chemical database and perform structural searches using Razi's extensions to the SQLAlchemy SQL Expression Language. The Python bindings for the RDKit libraries will be used in some data pre-processing steps, so you'll need to have them available on your system.

Download data
-------------

Download the `ChEMBLdb <https://www.ebi.ac.uk/chembldb/index.php>`_ database and decompress it::

    $ wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_08/chembl_08_chemreps.txt.gz
    $ gunzip chembl_08_chemreps.txt.gz


Database creation
-----------------

Then, a new chemical database is created::
 
    $ createdb -Udb_user -Ttemplate_rdkit chembldb_tutorial

This initialization step includes the installation of the RDKit PostgreSQL extension (for details about PostgreSQL+RDKit setup see the dedicated section of this document).

Connection to the database
--------------------------

Start your python interpreter and configure a database connection::

    >>> from sqlalchemy import create_engine
    >>> engine = create_engine('postgresql://db_user:db_password@host:1234/chembldb_tutorial')


Schema definition
-----------------

For this tutorial we will only need one single table in the database schema::

    >>> from sqlalchemy import Table, Column, Integer, String, MetaData
    >>> from razi.chemtypes import Molecule
    >>> metadata = MetaData()
    >>> compounds = Table('compounds', metadata,
    ...     Column('id', Integer, primary_key=True),
    ...     Column('chembl_id', String),
    ...     Column('structure', Molecule),
    ... )
    >>>  

and::

    >>> metadata.create_all(engine)

This last command creates the `compounds` table and also implicitly includes the creation of a structural index on the column with type `Molecule`.

The data from ChEMBL is formatted as a csv file, tab-separated and beginning with a header line::

    $ head -n3  chembl_08_chemreps.txt 
    chembl_id	chebi_id	molweight	canonical_smiles	inchi	inchi_key
     CHEMBL6582	100733	268.09877	NC(=O)c1cc(nnc1Cl)c2ccc(Cl)cc2	InChI=1/C11H7Cl2N3 <snip>
     CHEMBL6583	100744	342.41542	CN(C)c1cccc2c(cccc12)S(=O)(=O)Nc3cnc(C)cn3	In <snip>

Define a namedtuple matching the format of the csv records in the input file::

    >>> from collections import namedtuple
    >>> Record = namedtuple('Record', 'chembl_id, chebi_id, mw, smiles, inchi, inchi_key')
    >>> 

In this tutorial only a subset of this data will be used, and at the same time we want to make sure that troublesome SMILES strings are skipped (SMILES containing errors, compounds that are too big and/or other strings that the RDKit cartridge won't be able to process). File parsing and data filtering can be performed with a function similar to the following::

    import csv 
    from rdkit import Chem

    def read_chembldb(filepath, limit=0):
    
        inputfile = open(filepath, 'rb')
        reader = csv.reader(inputfile, delimiter='\t', skipinitialspace=True)
        reader.next() # skip header line
    
        for count, record in enumerate(map(Record._make, reader), 1):
    
	    smiles = record.smiles

            # skip problematic compounds
            if len(smiles) > 300: continue
            smiles = smiles.replace('=N#N','=[N+]=[N-]')
            smiles = smiles.replace('N#N=','[N-]=[N+]=')
            if not Chem.MolFromSmiles(smiles): continue
    
            yield {'chembl_id': record.chembl_id, 
	           'structure': smiles}
    
            if count == limit: 
	        break

The ``read_chembldb`` function above is a python generator, returning for each record a dictionary containing the `chembl_id` and `smiles` strings. With this function importing the compounds into the database reduces to a simple loop *(please note that depending on the available hardware resources importing the whole database may require a few hours; to keep this tutorial short we'll limit processing to the first 25K compounds)*::

    >>> conn = engine.connect()
    >>> for kw in read_chembldb('chembl_08_chemreps.txt', 25000):
    >>>     conn.execute(compounds.insert(), **kw)

and we can finally perform some queries. We can for example verify the number of compounds in the database::

    >>> engine.scalar(compounds.count())
    >>> 25000

or display the first 5 compounds::

    (CHEMBL6582) NC(=O)c1cc(-c2ccc(Cl)cc2)nnc1Cl
    (CHEMBL6583) Cc1cnc(NS(c2cccc3c(N(C)C)cccc23)(=O)=O)cn1
    (CHEMBL6584) CN(C)/C=N/c1nc(/N=C\N(C)C)c2c(ncc(Sc3cc(Cl)c(Cl)cc3)n2)n1
    (CHEMBL6585) CC12C(C[C@@H](I)[C@@H]1O)C1C(c3ccc(O)cc3CC1)CC2
    (CHEMBL6637) C/C(=C\Cn1oc(=O)[nH]c1=O)c1ccc(OCCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)cc1

Finally (and hopefully more interestingly), here's a first example of a chemical query, searching the database for a given substructure::

    In [12]: # which compounds contain 'c1cccc2c1nncc2' as a substructure?
    
    In [13]:
    (CHEMBL12112) CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1
    (CHEMBL26025) Cc1cccc(NC(=O)Nc2ccc3nnccc3c2)c1

Creating a chemical database
============================

This tutorial is based on a similar document available from the `RDKit wiki <http://code.google.com/p/rdkit/wiki/DatabaseCreation>`_ and it illustrates how to build a chemical database and perform structural searches using Razi. The Python bindings for the RDKit libraries will be used in some data pre-processing steps, so you'll need to have them available on your system.

Some basic understanding of SQLAlchemy is assumed. For additional details, please refer to the excellent `SQLALchemy ORM tutorial <http://www.sqlalchemy.org/docs/orm/tutorial.html>`_.

Download sample data
--------------------

Download the `ChEMBLdb <https://www.ebi.ac.uk/chembldb/index.php>`_ database and decompress it::

    $ wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_08/chembl_08_chemreps.txt.gz
    $ gunzip chembl_08_chemreps.txt.gz


Database creation
-----------------

Then, a new chemical database is created::
 
    $ createdb -Udb_user -Ttemplate_rdkit razi_tutorial

This initialization step includes the installation of the RDKit PostgreSQL extension (for details about PostgreSQL+RDKit setup see the dedicated section in this document).

*The database user name and password, together with the name of the database itself, are most likely to change in your case. Just replace them consistently with the values appropriate to your work environment.* 

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

    from sqlalchemy import Column, Integer, String
    from razi.orm import ChemColumn
    from razi.chemtypes import Molecule
    
    class Compound(Base):
        __tablename__='compounds'
        
        id = Column(Integer, primary_key=True)
        name = Column(String)
        structure = ChemColumn(Molecule)
        
        def __init__(self, name, structure):
            self.name = name
            self.structure = structure
            
        def __repr__(self):
            return '(%s) < %s >' % (self.name, self.structure)

to actually create the database schema, a method of the ``Base`` class ``metadata`` attribute is called::

    Base.metadata.create_all()

In the present case this last command creates the ``compounds`` table and also implicitly includes the creation of a structural index on the column with type ``Molecule``.

Inserting data
--------------

The data from ChEMBL is formatted as a csv file, tab-separated and beginning with a header line::

    $ head -n3  chembl_08_chemreps.txt 
    chembl_id	chebi_id	molweight	canonical_smiles	inchi	inchi_key
     CHEMBL6582	100733	268.09877	NC(=O)c1cc(nnc1Cl)c2ccc(Cl)cc2	InChI=1/C11H7Cl2N3 <snip>
     CHEMBL6583	100744	342.41542	CN(C)c1cccc2c(cccc12)S(=O)(=O)Nc3cnc(C)cn3	In <snip>

To semplify the processing of these data records we define a namedtuple matching the input format::

    from collections import namedtuple
    Record = namedtuple('Record', 'chembl_id, chebi_id, mw, smiles, inchi, inchi_key')

During this tutorial only a subset of this data is actually imported into the database, and at the same time we want to make sure that troublesome SMILES strings are skipped (SMILES containing errors, compounds that are too big and/or other strings that the RDKit cartridge won't be able to process). File parsing and data filtering can be performed with a function similar to the following::

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
    
            yield count, record.chembl_id, smiles
    
            if count == limit: 
	        break

The ``read_chembldb`` function above is a python generator, producing for each valid record a python tuple containing the record counter and the ``chembl_id`` and ``smiles`` strings.

With this function importing the compounds into the database reduces to a simple loop *(please note that depending on the available hardware resources importing the whole database may require a few hours; to keep this tutorial short we'll limit the processing to the first 25K compounds, a dataset size the usually corresponds to a few minutes)*::

    session = Session()
    for count, chembl_id, smiles in read_chembldb('chembl_08_chemreps.txt', 25000):
        compound = Compound(chembl_id, smiles)
	session.add(compound)
	if not count % 1000:
	    session.commit()

Querying the database
---------------------

Finally, we can perform some queries. We can for example verify the number of compounds actually loaded into the database::

    >>> print session.query(Compound).count()
    24956
    >>> 

or display the first 5 compounds::

    >>> for compound in session.query(Compound)[:5]:
    ...     print compound
    ... 
    (CHEMBL6582) < NC(=O)c1cc(-c2ccc(Cl)cc2)nnc1Cl >
    (CHEMBL6583) < Cc1cnc(NS(c2cccc3c(N(C)C)cccc23)(=O)=O)cn1 >
    (CHEMBL6584) < CN(C)/C=N/c1nc(/N=C\N(C)C)c2c(ncc(Sc3cc(Cl)c(Cl)cc3)n2)n1 >
    (CHEMBL6585) < CC12C(C[C@@H](I)[C@@H]1O)C1C(c3ccc(O)cc3CC1)CC2 >
    (CHEMBL6637) < C/C(=C\Cn1oc(=O)[nH]c1=O)c1ccc(OCCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)cc1 >
    >>> 

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
    >>> 

Please notice how the SQLAlchemy's ORM API allows the incremental specification of the filtering clause (or clauses) associated to the main selection query and how the ``subset`` instance is actually used twice (to compute the number of record matching the query and to retrieve the actual records). In addition to this, the returned records can also be used as the basis for further queries, also using the chemical functions provided by the database backend:

    >>> for compound in subset: 
    ...     # a query returning the computed molecular weight for each compound
    ...     print session.scalar(compound.structure.mw)
    ... 
    488.701
    278.315
    >>> 

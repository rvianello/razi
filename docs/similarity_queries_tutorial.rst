.. _similarity-queries-tutorial:

Example similarity queries
==========================

This tutorial is again based on a similar document available from the `RDKit wiki <http://code.google.com/p/rdkit/wiki/ExampleSimilarityQueries>`_ and it illustrates how to use Razi to perform structure similarity queries on a chemical database.


Database creation
-----------------

Create a new database for this tutorial::
 
    $ createdb -Udb_user -Ttemplate_rdkit razi_tutorial

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

The schema is composed of one single database entity, mapped to a python class. As described in the first tutorial, a base class is first defined::

    from sqlalchemy.ext.declarative import declarative_base
    Base = declarative_base(bind=engine)

then, the definition of the mapped entity follows::

    from sqlalchemy import Column, Integer, String
    from razi.orm import ChemColumn
    from razi.chemtypes import Molecule, BitFingerprint
    
    class Compound(Base):
        __tablename__='compounds'
        
        id = Column(Integer, primary_key=True)
        name = Column(String)
        structure = ChemColumn(Molecule)
        atompair = ChemColumn(BitFingerprint)
        torsion = ChemColumn(BitFingerprint)
        morgan = ChemColumn(BitFingerprint)
    
        def __init__(self, name, structure):
            self.name = name
            self.structure = structure
            self.atompair = self.structure.atompair_b()
            self.torsion = self.structure.torsion_b()
            self.morgan = self.structure.morgan_b(2)
        
        def __repr__(self):
            return '(%s) < %s >' % (self.name, self.structure)


and the database schema is created::

    Base.metadata.create_all()

In the present case this last command creates the ``compounds`` table and also implicitly includes the creation of indices on the columns with types :class:`~razi.chemtypes.Molecule` and  :class:`~razi.chemtypes.BitFingerprint`. Please notice how in the constructor the fingerprint fields are initialized by database backend expressions invoked on the ``structure`` column.

Inserting data
--------------

To populate the database the same data and code we used in the first tutorial is used again (this time we import a more extended dataset, but still small enough to keep the processing time acceptably short for a tutorial. Feel free to modify the number of compounds imported, just be aware that some results may change with the imported dataset)::

    session = Session()
    for count, chembl_id, smiles in read_chembldb('chembl_08_chemreps.txt', 50000):
        compound = Compound(chembl_id, smiles)
	session.add(compound)
    session.commit()


Querying the database
---------------------

We are now ready to reproduce the similarity queries executed in the original RDKit tutorial. 

Similarity search using Morgan fingerprints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tanimoto similarity using the current similarity cutoff::

    >>> # search compounds similar to the following (known to be in the db):
    >>> query_cmpnd = 'CN(C)/C=N/c1nc(/N=C\N(C)C)c2c(ncc(Sc3cc(Cl)c(Cl)cc3)n2)n1'
    >>> # compute a binary fingerprint for the query compound 
    >>> # (actually, it's the corresponding SQL expression whose value is to be computed by the db)
    >>> from razi.functions import functions
    >>> from razi.expression import TxtMoleculeElement
    >>> query_bfp = functions.morgan_b(TxtMoleculeElement(query_cmpnd), 2)
    >>> # determine the number of compunds with Tanimoto similarity above
    >>> # the current threshold value:
    >>> print session.query(Compound).filter(Compound.morgan.tanimoto_similar(query_bfp)).count()
    2
    
Or using the Dice similarity::

    >>> print session.query(Compound).filter(Compound.morgan.dice_similar(query_bfp)).count()
    6

Including the similarity values in the search results::

    >>> constraint = Compound.morgan.dice_similar(query_bfp)
    >>> dice_sml = Compound.morgan.dice_similarity(query_bfp).label('dice')
    >>> from sqlalchemy import desc
    >>> results = session.query(Compound, dice_sml).filter(constraint).order_by(desc(dice_sml))
    >>> for row in results: print row.Compound, row.dice
    (CHEMBL6584) < CN(C)/C=N/c1nc(/N=C\N(C)C)c2c(ncc(Sc3cc(Cl)c(Cl)cc3)n2)n1 > 1.0
    (CHEMBL6544) < Nc1nc(N)c2c(ncc(Sc3cc(Cl)c(Cl)cc3)n2)n1 > 0.666666666667
    (CHEMBL6618) < Nc1nc(N)c2c(ncc(Sc3cc4c(cccc4)cc3)n2)n1 > 0.52380952381
    (CHEMBL6465) < Nc1nc(N)c2c(ncc(Sc3cc(Cl)c(Cl)cc3Cl)n2)n1 > 0.506024096386
    (CHEMBL6631) < COc1ccc(Sc2cnc3c(c(N)nc(N)n3)n2)cc1 > 0.5
    (CHEMBL57035) < CCN(CC)CCCNc1ncc2cc(-c3c(Cl)cccc3Cl)c(/N=C\N(C)C)nc2n1 > 0.5


Similarity search using other fingerprints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At this point using the other fingerprint types basically only requires redefining the ``query_bfp`` fingerprint and the query constraint. For example, Tanimoto similarity between topological torsion fingerprints using the current similarity cutoff::

    >>> query_bfp = functions.torsion_b(TxtMoleculeElement(query_cmpnd))
    >>> constraint = Compound.torsion.tanimoto_similar(query_bfp)
    >>> tanimoto_sml = Compound.torsion.tanimoto_similarity(query_bfp).label('tanimoto')
    >>> results = session.query(Compound, tanimoto_sml).filter(constraint).order_by(desc(tanimoto_sml))

and Tanimoto similarity between atom-pair fingerprints using the current similarity cutoff is almost identical:: 

    >>> query_bfp = functions.atompair_b(TxtMoleculeElement(query_cmpnd))
    >>> constraint = Compound.atompair.tanimoto_similar(query_bfp)
    >>> tanimoto_sml = Compound.atompair.tanimoto_similarity(query_bfp).label('tanimoto')
    >>> results = session.query(Compound, tanimoto_sml).filter(constraint).order_by(desc(tanimoto_sml))


Changing the similarity cutoff values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The threshold values used by the Tanimoto and Dice filter operators are mapped to two expressions defined in module :py:mod:`razi.postgresql_rdkit`::

    >>> from razi.postgresql_rdkit import tanimoto_threshold, dice_threshold
    >>> session.scalar(tanimoto_threshold), session.scalar(dice_threshold)
    (u'0.5', u'0.5')
 
The same expressions provide a mechanism to set a different cutoff::

    >>> session.execute(tanimoto_threshold.set(0.65))
    <sqlalchemy.engine.base.ResultProxy object at 0x1bbc5a10>
    >>> session.scalar(tanimoto_threshold), session.scalar(dice_threshold)
    (u'0.65', u'0.5')


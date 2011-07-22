.. _structure-queries-tutorial:

Example structure queries
=========================

This tutorial is based on a similar document available from the `RDKit wiki <http://code.google.com/p/rdkit/wiki/ExampleStructureQueries>`_ and it illustrates how to use Razi to perform substructure and superstructure queries on a chemical database.

No dedicated database is created for this tutorial. The same database used in the ":ref:`similarity-queries-tutorial`" tutorial can be used instead. If you are no longer connected to the database, the connection configuration and schema definition must be entered again. If you are still connected just skip to section ":ref:`structure-queries-tutorial-querying`".

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

Then define the mapping to the database table::

    from sqlalchemy.ext.declarative import declarative_base
    Base = declarative_base(bind=engine)
    
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
        
        
        def __init__(self, name, structure):
            self.name = name
            self.structure = structure
            
        def __repr__(self):
            return '(%s) < %s >' % (self.name, self.structure)

.. _structure-queries-tutorial-querying:

Querying the database
---------------------

Substructure queries
^^^^^^^^^^^^^^^^^^^^

Retrieve the number of molecules containing a triazine:

    >>> constraint = Compound.structure.contains('c1ncncn1')
    >>> print session.query(Compound).filter(constraint).count()
    56
    
Retrieve the number of molecules containing a coumarin::

    >>> constraint = Compound.structure.contains('O=C1OC2=CC=CC=C2C=C1')
    >>> print session.query(Compound).filter(constraint).count()
    166

Get the first 10 of those::

    >>> for c in session.query(Compound).filter(constraint)[:10]: print c
    (CHEMBL58793) < OC(=O)CCCCc1cc(=O)oc2c1ccc(O)c2CN1CCCC1 >
    (CHEMBL56784) < [Na+].COc1ccc(-c2c3n(c4c(=O)oc5cc(OS([O-])(=O)=O)c(OC)cc5c42)CCc2cc(OC)c(OC)cc2-3)cc1O >
    (CHEMBL54909) < COc1cc2ccc(=O)oc2c(O)c1O >
    (CHEMBL50150) < COc1ccc(CCn2cc(-c3ccc(OC)c(OC)c3)c3c4c(oc(=O)c23)cc(OC)c(OC)c4)cc1OC >
    (CHEMBL50201) < CC(C)CCc1c(O)ccc2c1oc(=O)cc2 >
    (CHEMBL59509) < OC(=O)CCCCc1cc(=O)oc2c1ccc(O)c2CNc1ccccc1 >
    (CHEMBL57330) < CCCN(C1COc2cccc(OC)c2C1)CCCCNC(=O)c1c2c(oc(=O)c1)c1c3c(c2)CCCN3CCC1 >
    (CHEMBL57173) < C/C(CC/C=C(\C)C1=CC(=O)C(C)(C)O1)=C\COc1cc2oc(=O)ccc2cc1 >
    (CHEMBL57138) < COc1ccc(-c2c3n(c4c(=O)oc5cc(O)c(OC)cc5c42)CCc2c(OC)c(OC)c(OC)cc2-3)cc1O >
    (CHEMBL56918) < C/C(=C\COc1ccc2c(oc(=O)cc2)c1)C1=CC(=O)C(C)(C)O1 >


Including property filters
^^^^^^^^^^^^^^^^^^^^^^^^^^

Differently from the original RDKit tutorial, chemical descriptor were not introduced into the current database schema. Filtering based on chemical properties can still be introduced, with the difference that these properties are in this case computed while processing the query::

    >>> mw = Compound.structure.mw.label('mw')
    >>> logp = Compound.structure.logp.label('logp')
    >>> # compounds containing coumarin as substructure, with molecular weight
    >>> # not above 200, ordered by ascending estimated logp
    >>> subset = session.query(Compound, mw, logp).filter(constraint).filter(mw <= 200).order_by(logp)
    >>> for row in subset: print row.Compound.name, row.mw, row.logp
    CHEMBL32810 178.143 1.2042
    CHEMBL51628 162.144 1.4986
    CHEMBL12252 192.17 1.51262
    CHEMBL6466 146.145 1.793
    CHEMBL49732 176.171 1.8016
    CHEMBL12626 176.171 1.80702
    CHEMBL12208 176.171 1.80702
    CHEMBL12279 160.172 2.10142
    CHEMBL12636 190.198 2.11002
    CHEMBL19240 190.198 2.11544
    CHEMBL53569 186.166 2.5392
    CHEMBL6355 196.205 2.9462


Other kinds of structural searches
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Superstructure queries
~~~~~~~~~~~~~~~~~~~~~~

Look for molecules in the database that are substructures of a query (i.e. where the query is a superstructure of the database molecule)::

    >>> constraint = Compound.structure.contained_in('c1ccc(C(COC(c2c(=O)oc3c(ccc(O)c3)c2)=O)=O)cc1')
    >>> for c in session.query(Compound).filter(constraint)[:10]: print c
    (CHEMBL51628) < O=c1oc2cc(O)ccc2cc1 >
    (CHEMBL44857) < CCCOC(=O)C >
    (CHEMBL44215) < CCOC=O >
    (CHEMBL545) < CCO >
    (CHEMBL14688) < CO >
    (CHEMBL17564) < C >
    (CHEMBL15972) < O=Cc1ccccc1 >
    (CHEMBL14687) < CCCO >
    (CHEMBL16264) < CCOCC >
    (CHEMBL14079) < COC(=O)C >


SMARTS-based Queries
~~~~~~~~~~~~~~~~~~~~

``contains`` substructure queries are by default executed using SMILES semantics. In order to do SMARTS-based queries, one may use ``match``, as this example shows:  

    >>> constraint = Compound.structure.match('cc(c)NC(=O)N')
    >>> for c in session.query(Compound).filter(constraint)[:10]: print c
    (CHEMBL6997) < CSCC[C@H](NC(Nc1cc(C)ccc1)=O)C(=O)N[C@@H](CC(C)C)C(N[C@@H](Cc1ccccc1)C(O)=O)=O >
    (CHEMBL6500) < CCOC(c1ccc(NC(=O)Nc2c(C)cc3c(c2)C(C)(C)CC(C)(C)S3)cc1)=O >
    (CHEMBL6218) < COc1cc2c(c(N)nc(N3CCN(C(=O)Nc4ccccc4)CC3)n2)cc1OC >
    (CHEMBL7610) < COc1ccc(C[C@H](NC(Nc2cc3n(Cc4c(Cl)cccc4Cl)cc(CN4CCCC4)c3cc2)=O)C(N[C@@H](CCCNC(=N)N)C(NCc2ccccc2)=O)=O)cc1 >
    (CHEMBL7667) < CCCCNS(=NC(=O)Nc1ccc(Cl)cc1)(=O)c1ccc(C)cc1 >
    (CHEMBL7955) < CCNS(=NC(=O)Nc1ccc(Cl)cc1)(=O)c1ccc(C)cc1 >
    (CHEMBL7851) < Cc1c(Cl)c(C)cc(S(N)(=NC(=O)Nc2ccc(Cl)cc2)=O)c1 >
    (CHEMBL7627) < COc1ccc(C[C@H](NC(Nc2cc3n(Cc4ccc(F)cc4)cc(CNC4CCCC4)c3cc2)=O)C(N[C@@H](CCCN=C(N)N)C(NCc2ccccc2)=O)=O)cc1 >
    (CHEMBL7346) < CCOC(c1ccc(NC(=O)Nc2cc3c(cc2)N(C)C(C)(C)C=C3C)cc1)=O >
    (CHEMBL7520) < CSCC[C@H](NC(Nc1ccccc1)=O)C(N[C@@H](CC(C)C)C(N[C@@H](Cc1ccccc1)C(O)=O)=O)=O >


Exact match queries
~~~~~~~~~~~~~~~~~~~

Matching full structures is supported by using :class:`~razi.functions.functions.equals`::

    >>> print session.query(Compound).filter(Compound.structure.equals('c1ncncn1')).count()

or by just using the equality operator ``==``::

    >>> print session.query(Compound).filter(Compound.structure == 'c1ncncn1').count()


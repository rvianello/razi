from sqlalchemy import types as sqltypes
from sqlalchemy.sql import functions
from sqlalchemy.dialects import postgresql

from . import types

# GenericFunction subclasses strongly based on the geoalchemy2
# implementation

_FUNCTIONS = [
    #
    # mol constructors
    #

    ('mol', types.Mol, None,),

    ('mol_from_smiles', types.Mol, None,),
    ('mol_from_smarts', types.Mol, None,),
    ('mol_from_ctab', types.Mol, None,),
    ('mol_from_pkl', types.Mol, None,),

    ('mol_adjust_query_properties', types.Mol, None,),

    #
    # mol conversion functions
    #

    ('mol_to_smiles', sqltypes.String, None,), # also supports qmol input
    ('mol_to_cxsmiles', sqltypes.String, None,),
    ('mol_to_smarts', sqltypes.String, None,), # also supports qmol input
    ('mol_to_ctab', sqltypes.String, None,),
    ('mol_to_svg', sqltypes.String, None,),
    ('mol_to_pkl', postgresql.BYTEA, None,),

    #
    # mol descriptors
    #

    ('mol_amw', sqltypes.Float, None,),
    ('mol_logp', sqltypes.Float, None,),
    ('mol_fractioncsp3', sqltypes.Float, None,),
    ('mol_tpsa', sqltypes.Float, None,),
    ('mol_chi0n', sqltypes.Float, None,),
    ('mol_chi1n', sqltypes.Float, None,),
    ('mol_chi2n', sqltypes.Float, None,),
    ('mol_chi3n', sqltypes.Float, None,),
    ('mol_chi4n', sqltypes.Float, None,),
    ('mol_chi0v', sqltypes.Float, None,),
    ('mol_chi1v', sqltypes.Float, None,),
    ('mol_chi2v', sqltypes.Float, None,),
    ('mol_chi3v', sqltypes.Float, None,),
    ('mol_chi4v', sqltypes.Float, None,),
    ('mol_kappa1', sqltypes.Float, None,),
    ('mol_kappa2', sqltypes.Float, None,),
    ('mol_kappa3', sqltypes.Float, None,),

    ('mol_hba', sqltypes.Integer, None,),
    ('mol_hbd', sqltypes.Integer, None,),
    ('mol_numrotatablebonds', sqltypes.Integer, None,),
    ('mol_numatoms', sqltypes.Integer, None,),
    ('mol_numheavyatoms', sqltypes.Integer, None,),
    ('mol_numheteroatoms', sqltypes.Integer, None,),
    ('mol_numrings', sqltypes.Integer, None,),
    ('mol_numaromaticrings', sqltypes.Integer, None,),
    ('mol_numaliphaticrings', sqltypes.Integer, None,),
    ('mol_numsaturatedrings', sqltypes.Integer, None,),
    ('mol_numaromaticheterocycles', sqltypes.Integer, None,),
    ('mol_numaliphaticheterocycles', sqltypes.Integer, None,),
    ('mol_numsaturatedheterocycles', sqltypes.Integer, None,),
    ('mol_numaromaticcarbocycles', sqltypes.Integer, None,),
    ('mol_numaliphaticcarbocycles', sqltypes.Integer, None,),
    ('mol_numsaturatedcarbocycles', sqltypes.Integer, None,),
    ('mol_numheterocycles', sqltypes.Integer, None,),
    ('mol_numspiroatoms', sqltypes.Integer, None,),
    ('mol_numbridgeheadatoms', sqltypes.Integer, None,),

    ('mol_formula', sqltypes.String, None,),
    ('mol_inchi', sqltypes.String, None,),
    ('mol_inchikey', sqltypes.String, None,),
    ('mol_murckoscaffold', sqltypes.String, None,),
    ('mol_nm_hash', sqltypes.String, None,),

    #
    # qmol constructors
    #
    ('qmol_from_smiles', types.QMol, None,),
    ('qmol_from_ctab', types.QMol, None,),

    #
    # bfp constructors
    #

    ('bfp_from_binary_text', types.Bfp, None,),
    ('avalon_fp', types.Bfp, None,),
    ('layered_fp', types.Bfp, None,),
    ('maccs_fp', types.Bfp, None,),
    ('rdkit_fp', types.Bfp, None,),
    ('atompairbv_fp', types.Bfp, None,),
    ('featmorganbv_fp', types.Bfp, None,),
    ('morganbv_fp', types.Bfp, None,),
    ('torsionbv_fp', types.Bfp, None,),
    ('reaction_structural_bfp', types.Bfp, None,),

    #
    # bfp conversion functions
    #

    ('bfp_to_binary_text', postgresql.BYTEA, None,),

    #
    # bfp integer descriptors
    #

    ('bfp_size', sqltypes.Integer, None,),

    #
    # sfp constructors
    #

    ('atompair_fp', types.Sfp, None,),
    ('featmorgan_fp', types.Sfp, None,),
    ('morgan_fp', types.Sfp, None,),
    ('torsion_fp', types.Sfp, None,),
    ('reaction_difference_fp', types.Sfp, None,),

    #
    # sfp processing functions
    #

    ('add', types.Sfp, None,),
    ('subtract', types.Sfp, None,),
    ('all_values_gt', sqltypes.Boolean, None,),
    ('all_values_lt', sqltypes.Boolean, None,),

    #
    # reaction constructors
    #

    ('reaction', types.Reaction, None,),

    ('reaction_from_smiles', types.Reaction, None,),
    ('reaction_from_smarts', types.Reaction, None,),
    ('reaction_from_ctab', types.Reaction, None,),

    #
    # reaction conversion functions
    #

    ('reaction_to_smiles', sqltypes.String, None,),
    ('reaction_to_smarts', sqltypes.String, None,),
    ('reaction_to_ctab', sqltypes.String, None,),
    ('reaction_to_svg', sqltypes.String, None,),

    #
    # reaction descriptors
    #

    ('reaction_numreactants', sqltypes.Integer, None,),
    ('reaction_numproducts', sqltypes.Integer, None,),
    ('reaction_numagents', sqltypes.Integer, None,),

    #
    # similarity functions
    #

    ('tanimoto_sml', sqltypes.Float, None,),
    ('dice_sml', sqltypes.Float, None,),
    ('tversky_sml', sqltypes.Float, None,),

    ('tanimoto_dist', sqltypes.Float, None,),
    ('dice_dist', sqltypes.Float, None,),

    #
    # utility functions
    #

    ('rdkit_version', sqltypes.String, None,),
    ('is_valid_smiles', sqltypes.Boolean, None,),
    ('is_valid_smarts', sqltypes.Boolean, None,),
    ('is_valid_ctab', sqltypes.Boolean, None,),
    ('is_valid_mol_pkl', sqltypes.Boolean, None,),

    #
    # substructure search
    #
    ('substruct_count', sqltypes.Integer, None,),

]

# Iterate through _FUNCTIONS and create GenericFunction classes dynamically
for name, type_, doc in _FUNCTIONS:
    attributes = {'name': name, 'inherit_cache': True}
    docs = []

    if doc is not None:
        docs.append(doc)

    if type_ is not None:
        attributes['type'] = type_

        type_str = '{0}.{1}'.format(type_.__module__, type_.__name__)
        docs.append('Return type: :class:`{0}`.'.format(type_str))

    if len(docs) != 0:
        attributes['__doc__'] = '\n\n'.join(docs)

    globals()[name] = type(name, (functions.GenericFunction,), attributes)

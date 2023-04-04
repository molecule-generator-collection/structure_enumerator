import argparse
import itertools
import os
import random
import sys
import yaml

from rdkit import Chem, Geometry
from rdkit.Chem import rdDepictor, AllChem
rdDepictor.SetPreferCoordGen(True)


def get_parser():
    parser = argparse.ArgumentParser(
        description="",
        usage=f"python {os.path.basename(__file__)} -c CONFIG_FILE"
    )
    parser.add_argument(
        "-c", "--config", type=str, required=True,
        help="path to a config file"
    )
    return parser.parse_args()


def enumerate_all_products(core, rgroup_dict, output_format, randomOrder=False):
    """
    The original function was written by G. Landrum in the following RDKit blog post.
    Some modifications were done by S. Ishida.
    URL: https://greglandrum.github.io/rdkit-blog/posts/2022-03-14-rgd-and-molzip.html
    """
    order = itertools.product(*rgroup_dict.values())
    if randomOrder:
        order = list(order)
        random.shuffle(order)

    if output_format == "sdf":
        conformer = core.GetConformer()
        corePos = {i: Geometry.Point2D(conformer.GetAtomPosition(i))
                   for i in range(conformer.GetNumAtoms())}
        
    for tpl in order:
        tm = Chem.RWMol(core)
        for r in tpl:
            tm.InsertMol(r)
        prod = Chem.molzip(tm)
        if prod is not None:
            if output_format == 'sdf':
                rdDepictor.Compute2DCoords(prod, canonOrient=False, coordMap=corePos)
            yield prod

def main():
    args = get_parser()
    with open(args.config, 'r') as f:
        conf = yaml.load(f, Loader=yaml.SafeLoader)

    if conf['output_format'] not in {'smi', 'sdf'}:
        print('[ERROR] Please specify the output format from [smi, sdf].', 
              f'Your specified format: {conf["output_format"]}')
        sys.exit(1)

    suppl = Chem.SmilesMolSupplier(conf['core'], titleLine=False, nameColumn=-1)
    core = next(suppl)
    AllChem.Compute2DCoords(core)

    r_group_mol_dict = {}
    for r, fname in conf['rgroup'].items():
        with open(fname, 'r') as f:
            smiles_list = [l.strip('\n') for l in f.readlines()]
        rnum_mapped_mol_list = []
        for smi in smiles_list:
            map_num = r.replace('r', '')
            mol = Chem.MolFromSmiles(smi.replace('*', f'[*:{map_num}]'))
            rnum_mapped_mol_list.append(mol)
        r_group_mol_dict[r] = rnum_mapped_mol_list 

    product_generator = enumerate_all_products(core, r_group_mol_dict, conf['output_format'], randomOrder=False)
    mol_products = []
    seen_products = set()
    for prod in product_generator:
        Chem.SanitizeMol(prod)
        smi = Chem.MolToSmiles(prod)
        if smi not in seen_products:
            mol_products.append(prod)
            seen_products.add(smi)
        
    os.makedirs(os.path.dirname(conf['output']), exist_ok=True)
    if conf['output_format'] == 'smi':
        with open(conf['output'], 'w') as f:
            f.writelines('\n'.join(seen_products))
    elif conf['output_format'] == 'sdf':
        with Chem.SDWriter(conf['output']) as w:
            for m in mol_products:
                w.write(m)


if __name__ == "__main__":
    main()

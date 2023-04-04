# structure_enumerator

## Requirements

- Python: 3.10
- RDKit: 2022.9.5

## How to install

```python
pip install rdkit==2022.9.5
```

## How to use

Prepare a config file, and then just run the following command:

```bash
python run.py -c config.yaml
```

### How to prepare configuration file

Specify a core structure and R-group list in a configuration file as follows.

```yaml
core: data/core_structure.smi
rgroup:
  r1: data/r-group01.smi
  r2: data/r-group02.smi
  r3: data/r-group03.smi
output: result/enumerated_structures.sdf
output_format: sdf
```

> **Note**
> Only SMILES format is supported as an input format.

> **Note**
> Number of `rgroup` key (r`1`, r`2`, etc.) must corresponds to atom-map number of dummy atoms in the core structure.

You can specify `smi` or `sdf` as an `output_format`.

- If you want to enumerate all structures quickly, specify `smi`.
- If you want to preserve the positions of the core atoms, specify `sdf`.

import esm
import esm.inverse_folding
import torch

from utils import deep_mutational_scan


if __name__ == '__main__':
    pdb_files = [
        'data/1RLC.pdb',
        'data/4RUB_AS.pdb',
    ]
    pdb_chain_names = [
        {
            'NtRbcL': 'L',
            'NtRbcS': 'S',
        },
        {
            'NtRbcL': 'A',
            'NtRbcS': 'S',
        },
    ]

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()
    if torch.cuda.is_available():
        model = model.cuda()

    for pdb_file, chains in zip(pdb_files, pdb_chain_names):
        pdb_prefix = pdb_file.split('/')[-1].split('.')[0]
    
        structure = esm.inverse_folding.util.load_structure(pdb_file)
        coords, native_seqs = esm.inverse_folding.multichain_util.extract_coords_from_complex(structure)
    
        with open(f'target/results_esmif1_complex_all_{pdb_prefix}.txt', 'w') as f_all, \
             open(f'target/results_esmif1_complex_better_{pdb_prefix}.txt', 'w') as f_better:
            f_all.write('chain\tmutation\tlog_likelihood\tbetter\n')
    
            for name, chain_id in chains.items():
                f_better.write(f'{name}:\n')
    
                native_seq = native_seqs[chain_id]
                wildtype_ll, _ = esm.inverse_folding.multichain_util.score_sequence_in_complex(
                    model, alphabet, coords, chain_id, native_seq
                )
    
                for pos, wt, mt in deep_mutational_scan(native_seq):
                    mutant_seq = native_seq[:pos] + mt + native_seq[(pos + 1):]
                    mutant_ll, _ = esm.inverse_folding.multichain_util.score_sequence_in_complex(
                        model, alphabet, coords, chain_id, mutant_seq
                    )
                    f_all.write(f'{name}\t{wt}{pos + 1}{mt}\t{mutant_ll}\t{mutant_ll > wildtype_ll}\n')
                    if mutant_ll > wildtype_ll:
                        f_better.write(f'\t{wt}{pos + 1}{mt}\n')

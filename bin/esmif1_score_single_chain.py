import esm
import esm.inverse_folding
import torch

from utils import deep_mutational_scan


if __name__ == '__main__':
    pdb_file = 'data/1RLC.pdb'

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()
    if torch.cuda.is_available():
        model = model.cuda()

    chains = {
        'NtRbcL': 'L',
        'NtRbcS': 'S',
    }
    with open('target/results_esmif1_single_chain_all.txt', 'w') as f_all, \
         open('target/results_esmif1_single_chain_better.txt', 'w') as f_better:
        f_all.write('chain\tmutation\tlog_likelihood\tbetter\n')

        for name, chain_id in chains.items():
            f_better.write(f'{name}:\n')

            coords, native_seq = esm.inverse_folding.util.load_coords(pdb_file, chain_id)

            wildtype_ll, _ = esm.inverse_folding.util.score_sequence(
                model, alphabet, coords, native_seq
            )

            for pos, wt, mt in deep_mutational_scan(native_seq):
                mutant_seq = native_seq[:pos] + mt + native_seq[(pos + 1):]
                mutant_ll, _ = esm.inverse_folding.util.score_sequence(
                    model, alphabet, coords, mutant_seq
                )
                f_all.write(f'{name}\t{wt}{pos + 1}{mt}\t{mutant_ll}\t{mutant_ll > wildtype_ll}\n')
                if mutant_ll > wildtype_ll:
                    f_better.write(f'\t{wt}{pos + 1}{mt}\n')

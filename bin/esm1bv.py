import numpy as np
import scipy.special
import sys
import torch


def get_model_name(name):
    if name == 'esm1b':
        from fb_model import FBModel
        model = FBModel(
            'esm1b_t33_650M_UR50S',
            repr_layer=[-1],
        )
    elif name.startswith('esm1v'):
        from fb_model import FBModel
        model = FBModel(
            'esm1v_t33_650M_UR90S_' + name[-1],
            repr_layer=[-1],
        )
    else:
        raise ValueError('Model {} not supported'.format(name))

    return model


def encode(seq, model):
    return model.encode(seq)


def decode(embedding, model, exclude=set()):
    if exclude == 'unnatural':
        exclude = set([
            'B', 'J', 'O', 'U', 'X', 'Z', '-', '.',
        ])
    
    logits = model.decode(embedding)

    assert(logits.shape[0] == embedding.shape[0])
    assert(logits.shape[1] == len(model.alphabet_.all_toks))

    valid_idx = [
        idx for idx, tok in enumerate(model.alphabet_.all_toks)
    ]
    logits = logits[:, valid_idx]

    argmax = ''.join([
        model.alphabet_.all_toks[valid_idx[tok_idx]]
        if ('<' not in model.alphabet_.all_toks[valid_idx[tok_idx]] and
            model.alphabet_.all_toks[valid_idx[tok_idx]] not in exclude) else '.'
        for tok_idx in np.argmax(logits, 1)
    ])

    return argmax


def reconstruct(seq, model, encode_kwargs={}, decode_kwargs={}):
    return decode(
        encode(seq, model, **encode_kwargs),
        model, **decode_kwargs
    )


def diff(seq_old, seq_new, start=0, end=None):
    different_muts = []
    for idx, (ch_old, ch_new) in enumerate(
            zip(seq_old, seq_new)
    ):
        if idx < start:
            continue
        if end is not None and idx >= end:
            continue
        if ch_new == '.':
            continue
        if ch_old != ch_new:
            different_muts.append((idx, ch_old, ch_new))
    return different_muts


def reconstruct_multi_models(
        wt_seq,
        model_names=[
            'esm1b',
            'esm1v1',
            'esm1v2',
            'esm1v3',
            'esm1v4',
            'esm1v5',
        ],
        alpha=None,
        return_names=False,
):
    mutations_models, mutations_model_names = {}, {}
    for model_name in model_names:
        model = get_model_name(model_name)
        if alpha is None:
            wt_new = reconstruct(
                wt_seq, model, decode_kwargs={ 'exclude': 'unnatural' }
            )
            mutations_model = diff(wt_seq, wt_new)
        else:
            mutations_model = soft_reconstruct(
                wt_seq, model, alpha=alpha,
            )
        for mutation in mutations_model:
            if mutation not in mutations_models:
                mutations_models[mutation] = 0
                mutations_model_names[mutation] = []
            mutations_models[mutation] += 1
            mutations_model_names[mutation].append(model.name_)
        del model

    if return_names:
        return mutations_models, mutations_model_names

    return mutations_models


if __name__ == '__main__':
    sequences = {
        'NtRbcL': 'MSPQTETKASVGFKAGVKEYKLTYYTPEYQTKDTDILAAFRVTPQPGVPPEEAGAAVAAESSTGTWTTVWTDGLTSLDRYKGRCYRIERVVGEKDQYIAYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPPAYVKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYECLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEALYKAQAETGEIKGHYLNATAGTCEEMIKRAVFARELGVPIVMHDYLTGGFTANTSLAHYCRDNGLLLHIHRAMHAVIDRQKNHGIHFRVLAKALRMSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFVEQDRSRGIYFTQDWASLPGVLPVASGGIHVWHMPALTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVALEACVKARNEGRDLAQEGNEIIREACKWSPELAAACEVWKEIVFNFAAVDVLDK',
        'NtRbcS': 'MQVWPPINKKKYETLSYLPDLSQEQLLSEVEYLLKNGWVPCLEFETEHGFVYRENNKSPGYYDGRYWTMWKLPMFGCTDATQVLAEVEEAKKAYPQAWIRIIGFDNVRQVQCISFIAYKPEGY',
    }
    
    for name, seq in sequences.items():
        print(f'{name}:')
        mutations_models = reconstruct_multi_models(seq)
        for k, v in sorted(mutations_models.items(), key=lambda item: -item[1]):
            mut_str = f'{k[1]}{k[0] + 1}{k[2]}'
            print(f'\t{mut_str}\t{v}')

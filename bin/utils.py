AAs = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
]

def deep_mutational_scan(sequence, exclude_noop=True):
    for pos, wt in enumerate(sequence):
        for mt in AAs:
            if exclude_noop and wt == mt:
                continue
            yield (pos, wt, mt)

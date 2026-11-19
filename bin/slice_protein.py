#!/usr/bin/env python3
"""sliceProtein: the pieces of a protein (or of a slice of one) left after removing annotated regions.

Shared by extract_unannotated_parquet_slices.py (Pfam regions) and mask_mgnifam_regions.py (MGnifam regions).
"""


def sliceProtein(mgyp, sequence, metadata, min_sequence_length, start=1):
    """Returns FASTA records of the gaps between the merged regions of metadata ("p" and "m" lists whose last two
    values are a 1-based inclusive start and end in sequence coordinates) that are at least min_sequence_length long.

    start is the protein coordinate of sequence[0], so headers ">{mgyp}_{s}_{e}" stay in protein coordinates when
    sequence is itself a slice.
    """
    pfams = metadata.get("p", [])
    mgnifams = metadata.get("m", [])

    sorted_regions = sorted((region[-2], region[-1]) for region in pfams + mgnifams)
    merged_regions = []

    for region_start, end in sorted_regions:
        if not merged_regions or region_start > merged_regions[-1][1]:
            merged_regions.append([region_start, end])
        else:
            merged_regions[-1][1] = max(merged_regions[-1][1], end)

    offset = start - 1
    sliced_sequences = []
    current_start = 1
    for region_start, end in merged_regions:
        if region_start - current_start >= min_sequence_length:
            sliced_sequences.append(
                f">{mgyp}_{current_start + offset}_{region_start - 1 + offset}\n"
                f"{sequence[current_start - 1 : region_start - 1]}\n"
            )
        current_start = end + 1

    if len(sequence) - current_start + 1 >= min_sequence_length:
        sliced_sequences.append(
            f">{mgyp}_{current_start + offset}_{len(sequence) + offset}\n{sequence[current_start - 1 :]}\n"
        )

    return sliced_sequences

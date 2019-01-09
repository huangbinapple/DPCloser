import sys
import os
import getopt

import cmap_file
import xmap_file


def print_help():
    body = '<contig_cmap> <contig_key> <ref_cmap> <align_dir> <output_file>'
    print("python3 {} {}".format(__file__, body))

def compare_refs(refs_, refs, align_option=True):
    transform_arrays = {}
    for ref_id, cmap_ in refs_.items():
        if len(cmap_.positions) != len(refs[ref_id].positions):
            transform_array = []
            transform_arrays[ref_id] = transform_array
            index = 0
            positions = refs[ref_id].positions
            for position_ in cmap_.positions:
                if position_ == positions[index]:
                    transform_array.append(index)
                else:
                    if align_option:
                        transform_array.append(index)
                    else:
                        transform_array.append(index + 1)
                    index += 1
                index += 1
            assert len(transform_array) == len(cmap_.positions)
        else:
            transform_arrays[ref_id] = list(range(len(cmap_.positions)))
    return transform_arrays

def retrieve_fragments(fragments, indexs):
    result = []
    if indexs[1] > indexs[0]:
        result.append(fragments[
            indexs[0] :\
            indexs[0] + 1
        ])
        print(indexs[0], indexs[0] + 1)
        for i in range(1, len(indexs)):
            result.append(fragments[
                indexs[i - 1] + 1 :\
                indexs[i] + 1
            ])
            print(indexs[i - 1] + 1, indexs[i] + 1)
        result.append(fragments[
            indexs[i] + 1 : \
            indexs[i] + 2
            ])
        print(indexs[i] + 1, indexs[i] + 2)
    else:
        result.append(fragments[
            indexs[0] + 1 :
            indexs[0] : -1
        ])
        for i in range(0, len(indexs) - 1):
            result.append(fragments[
                indexs[i] : \
                indexs[i + 1] : -1
            ])
        result.append(fragments[
            0 : indexs[-1] + 1
        ][::-1])
    return result

def main():
    interface = 'h'
    options, args = getopt.getopt(sys.argv[1:], interface)
    for option, value in options:
        if option == '-h':
            print_help()
            sys.exit()
        else:
            print_help()
            sys.exit()
    contig_cmap, contig_key, ref_cmap, align_dir, output_file = args

    # read key file.
    contig_id2node_name = {}
    with open(contig_key) as f:
        for i in range(4):
            f.readline()
        for line in f:
            tokens = line.split('\t')
            contig_id2node_name[tokens[0]] = tokens[1]

    # read contig_camp file.
    contigs = cmap_file.read_file(contig_cmap)
    print(contigs['363'].fragments)
    print(len(contigs['363'].fragments))

    # read ref cmap file.
    refs = cmap_file.read_file(ref_cmap)

    # Get the three alignment file.
    refs_, alignment_file = None, None
    for file_name in os.listdir(align_dir):
        if file_name.endswith('_r.cmap'):
            refs_ = cmap_file.read_file(align_dir + '/' + file_name)
        elif file_name.endswith('xmap'):
            alignment_file = align_dir + '/' + file_name
    assert refs_, alignment_file

    # Compare refs and refs_.
    transform_arrays = compare_refs(refs_, refs)
    transform_arrays_rev = compare_refs(refs_, refs, False)

    with open(output_file, 'w') as fout:
        for id_query, id_ref, oritation, alignment_info in \
                xmap_file.read_file(alignment_file):
            node_name = contig_id2node_name[id_query]
            if node_name not in ('node252292'):
                pass
                # continue
            ref_indexs, query_indexs = list(zip(*alignment_info))
            ref_indexs = list(ref_indexs)
            query_indexs = list(query_indexs)
            # fout.write('\t'.join(map(str, ref_indexs)))
            # fout.write('\n')
            # fout.write('\t'.join(map(str, query_indexs)))
            # fout.write('\n')
            # From 1-indexed to 0-indexed.
            for indexs in ref_indexs, query_indexs:
                for i in range(len(indexs)):
                    indexs[i] -= 1
            # fout.write('\t'.join(map(str, ref_indexs)))
            # fout.write('\n')
            # fout.write('\t'.join(map(str, query_indexs)))
            # fout.write('\n')
            if oritation == '+':
                ref_indexs = [transform_arrays[id_ref][i] \
                    for i in ref_indexs]
            else:
                ref_indexs = [transform_arrays_rev[id_ref][i] \
                    for i in ref_indexs]
            
            ref_fragments = refs[id_ref].fragments
            query_fragments = contigs[id_query].fragments
            # Write first line.
            node_length = int(contigs[id_query].length)
            direction = '1' if oritation == '+' else '0'
            fout.write(' '.join(map(str, (node_name, node_length,
                direction, ref_indexs[0], ref_indexs[-1] + 1))))
            fout.write('\n')
            # Write second line.
            fout.write('XXX\n')
            # Write thrid line (subject).
            fragment_groups = retrieve_fragments(ref_fragments, 
                ref_indexs)
            # fout.write(' '.join(map(str, fragment_groups)))
            fout.write('; '.join(
                (' '.join(map(lambda x: str(x) + ',400', ele)) for \
                    ele in fragment_groups)
            ))
            fout.write('\n')
            # Write fourth line (query).
            fragment_groups = retrieve_fragments(query_fragments,
                query_indexs)
            fout.write('; '.join(
                (' '.join(map(str, ele)) for ele in fragment_groups)
            ))
            # fout.write(' '.join(map(str, fragment_groups)))
            fout.write('\n')

            # Write fifth line.
            fout.write('XXX\n')


if __name__ == "__main__":
    main()


    
# python3 getTax.py <accession_file> <accession_col> <acc2taxid> <ranklineage> <outacc2tax.txt>
# OUTPUT: outacc2tax.txt
# accession|taxid|lineage_information
import logging
import sys
import tlog


def overlap(a, b, min_overlap_len=1):
    if a[1] < b[0] + min_overlap_len or a[0] > b[1] - min_overlap_len:
        return 0
    else:
        return 1


def gettax(input_acc_file, acc2taxid_file, ranklineage_file, acctax_output, acc_col=2):
    load_acc_n = 0
    # used acc add
    used_accession = dict()
    with open(input_acc_file) as f:
        for line in f:
            temp = line.rstrip().split()
            acc = temp[acc_col - 1].split('.')[0]
            if acc not in used_accession:
                # ['taxid','linkage']
                used_accession[acc] = ['', '']
                load_acc_n += 1

    # acc to taxid
    load_taxid_n = 0
    success_acc2taxid_n = 0
    used_tax = dict()
    with open(acc2taxid_file) as f:
        taxidx = 2
        if acc2taxid_file.endswith('FULL'):  # for nr
            taxidx = 1
        for line in f:
            temp = line.rstrip().split()
            acc = temp[0].split('.')[0]
            if acc in used_accession:
                used_accession[acc][0] = temp[taxidx]
                if temp[taxidx] not in used_tax:
                    used_tax[temp[taxidx]] = dict()
                    load_taxid_n += 1
                used_tax[temp[taxidx]][acc] = 1
                success_acc2taxid_n += 1

    # taxid to lineage
    with open(ranklineage_file) as f:
        for line in f:
            temp = [x.strip() for x in line.split('|')]
            if temp[0] in used_tax:
                for acc in used_tax[temp[0]]:
                    used_accession[acc][1] = '|'.join(temp[1:])

    with open(acctax_output, 'w') as fout:
        for key in used_accession:
            fout.write(key+"|"+used_accession[key][0]+"|"+used_accession[key][1]+'\n')
    return load_acc_n, load_taxid_n, success_acc2taxid_n


def detectcontamination(blastout, tax_accession_out, idt, lth, whitelist,
                        fasta_seq,
                        white_seq_path,
                        black_seq_path,
                        black_table,
                        force_remove):

    # load acctax
    acc2tax = dict()
    blackacc = dict()
    whiteacc_n = 0
    with open(tax_accession_out) as fta:
        for line in fta:
            temp = line.rstrip().split('|')
            acc2tax[temp[0]] = temp[1:]
            if whitelist not in temp[2:]:
                blackacc[temp[0]] = 1
            else:
                whiteacc_n += 1
    if whiteacc_n == 0:
        # warning
        logging.warning("# No hit of white list, "
                        "which means only contamination exists,"
                        "please check the input of '-wl,--white_list' !")

    # load blastout as filter
    blackquery = dict()  # {(query: {(qstart,qend):acc)}}
    whitequery = dict()  # {(query: {(qstart,qend):acc)}}
    with open(blastout) as fbo:
        for line in fbo:
            temp = line.rstrip().split('\t')
            idt_idx = 9
            lth_idx = 4
            qstart_idx = 5
            qend_idx = 6
            # high similar record
            if float(temp[idt_idx]) >= idt and float(temp[lth_idx]) >= lth:
                acc = temp[1].split('.')[0]
                qstart = float(temp[qstart_idx])
                qend = float(temp[qend_idx])
                query = temp[0]
                if acc in blackacc:
                    if query not in blackquery:
                        blackquery[query] = dict()
                    blackquery[query][(qstart, qend)] = acc
                else:
                    if query not in whitequery:
                        whitequery[query] = dict()
                    whitequery[query][(qstart, qend)] = acc

        blackquery_flag = dict()
        for query in blackquery:
            blackquery_flag[query] = 1
            # remove similar records with region both in white/black query
            if force_remove == 0:
                if query in whitequery:
                    for bregion in blackquery[query]:
                        if blackquery_flag[query] == 0:
                            break
                        for wregion in whitequery[query]:
                            # at least 10bp overlap
                            if overlap(bregion, wregion, 10):
                                blackquery_flag[query] = 0
                                break

        # write
        with open(fasta_seq) as f:
            with open(white_seq_path, 'w') as fw:
                with open(black_seq_path, 'w') as fb:
                    with open(black_table, 'w') as ft:
                        for line in f:
                            if line.startswith('>'):
                                tag = line[1:].rstrip()
                                name = tag.split(' ')[0]
                                toblackseq = 0
                                if name in blackquery:
                                    if blackquery_flag[name] == 1:
                                        toblackseq = 1
                                        accs = set([blackquery[name][region] for region in blackquery[name]])
                                        lineages = ['|'.join(acc2tax[acc]) for acc in accs if acc in acc2tax]
                                        ft.write('\t'.join([name, ','.join(accs), ','.join(lineages)]) + '\n')
                            # sequences with no blast hits will be writen in while list file
                            if toblackseq == 1:
                                fb.write(line)
                            else:
                                fw.write(line)


if __name__ == "__main__":
    try:
        input_acc = sys.argv[1]
        acccol = int(sys.argv[2])
        acc2taxid = sys.argv[3]
        ranklineage = sys.argv[4]
        output = sys.argv[5]
    except Exception as e:
        print(e)
        print('''
        python3 getTax.py <accession_file> <accession_col> <acc2taxid> <ranklineage> <output_file>
    
        <accession_file> - output such as blast
        <accession_col>  - the colnum of accession (e.g. blast: 2)
        <acc2taxid> - accesiontotax file
        <ranklineage> - rank lingeage file
        ''')
        exit()
    gettax(input_acc, acc2taxid, ranklineage, output, acccol)

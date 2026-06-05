import sys
import re

def print_help():
    print("Make fasta from matching result.\nusage:\n\tpython make_fa.py row.fa result.txt out.fa\n"
          "\trow.fa: row assembly fasta\n\tresult.txt: matching result from program matching\n\tout.fa: output fasta")

# ????????????
def complement(seq):
    """
    ????DNA???????????
    """
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return "".join([comp_dict.get(base, base) for base in seq])

def parse_fasta(file_path):
    """
    ????FASTA??,??????,
    ????ID,????????
    """
    sequences = {}
    current_id = None
    current_seq_list = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = "".join(current_seq_list)
                current_id = line[1:].split()[0]  # ??????,?????????ID
                current_seq_list = []
            else:
                current_seq_list.append(line)
        
        # ?????????????
        if current_id:
            sequences[current_id] = "".join(current_seq_list)
            
    return sequences

try:
    fain = sys.argv[1]
    if fain == "help":
        print_help()
        exit(0)

    oin = sys.argv[2]
    orderin = open(oin)
    faout_f = sys.argv[3]
    faout = open(faout_f, "w")
    
    # ?????????? SeqIO
    record_dict = parse_fasta(fain)
    h_appends = []

    # ???????????????? t[0:-1] ???? t[0:-1]
    tmp_set = set()

    # ??? orderin,???????????ID(t[0:-1])??????
    # ???????,??????????
    orderin_lines = []
    for line in orderin.readlines():
        tmp = re.split("\t", line.strip())
        first_t = tmp[0][0:-1]  # ????? t[0:-1]
        if first_t not in tmp_set:
            tmp_set.add(first_t)
            orderin_lines.append(line)
        else:
            print(f"Skipping duplicate entry for {first_t}")

    # ?? orderin ???,????
    orderin.seek(0)
    # ?????? orderin_lines
    for line in orderin_lines:
        if line.startswith("iter") or line.startswith("self"):
            continue
        seq = ""
        tmp = re.split("\t", line.strip())
        for t in tmp:
            try:
                # ?????????????
                tmp_seq = record_dict[t[0:-1]]
                if t[-1] == '-':
                    # ????????:???,???
                    tmp_seq = complement(tmp_seq)[::-1]
                seq += tmp_seq
                h_appends.append(t[0:-1])
            except KeyError:
                print(f"Sequence {t[0:-1]} not found in the record_dict.")
        faout.write(">" + "_".join(tmp) + "\n")
        faout.write(seq + "\n")

    # ????? record_dict ????
    for l in record_dict.keys():
        if l in h_appends:
            continue
        faout.write(">" + l + "\n")
        faout.write(record_dict[l] + "\n")

    faout.close()

except IndexError:
    print_help()
except FileNotFoundError:
    print("Error: One or more input files not found. Please check the file paths.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")


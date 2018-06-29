import io_util as io
import methyl_noe_parser

def parseNOEData(noe_files):
    """
    New NOE scoring function
    :param noe_files:
    :return:
    """
    print noe_files
    noe_lines = io.readFile(noe_files)
    noe_data = []
    total_noe_count = 0.0
    for noe in noe_lines:
        if len(noe) <= 1:
            pass
        else:
            res1, atm1, res2, atom2, noe, tol = noe.split()
            total_noe_count += 1.0
            if float(tol) == 0.0:
                tol = 0.1
            noe_data.append([int(res1), atm1, int(res2), atom2, float(noe), float(tol)])
    return noe_data, total_noe_count

def parseBMRBblockMR(block_file, aa_seq, seq_correction):

    return methyl_noe_parser.methylNOEParser(block_file, aa_seq, seq_correction)

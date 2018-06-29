import math
import numpy
import copy
import utility.io_util as io
import bbRMSD
import os


def getHCoorMatrix(ss_list, smotif):
    # change this to match the format from bbrmsd
    noe_matrix = {}
    count = 0
    for i in range(1, len(smotif), 5):
        coo1 = [[smotif[i][3]], [smotif[i][4]], [smotif[i][5]], ['H']]
        noe_matrix[ss_list[count]] = coo1
        count += 1
    return noe_matrix


def processBBC(bbc):
    tlen = len(bbc)
    x, y, z = [None] * tlen, [None] * tlen, [None] * tlen
    for i in range(0, tlen):
        x[i] = bbc[i][3]
        y[i] = bbc[i][4]
        z[i] = bbc[i][5]
    return [x, y, z]


def getBackboneCoors(ss_list, smotif):
    backbone = {}
    count = 0
    j = 0
    for i in range(0, len(ss_list)):
        backbone[ss_list[count]] = smotif[j:j + 5]
        j = j + 5
        count += 1
    return backbone


def getILVARotamers(res_type, bbc, spin):
    """

    :param res_type:
    :param bbc:
    :param spin:
    :return:
    """

    cwd = (os.path.dirname(os.path.realpath(__file__)))
    file_name = cwd + '/sidechainRotamers/' + res_type + "_sc.pickle"
    rotamers = io.readPickle(file_name)
    rmsd_cutoff = 0.1
    bbc = processBBC(bbc)
    spin_coors, cluster_coors = bbRMSD.bbrmsd(bbc, rotamers, rmsd_cutoff, spin, res_type)

    return spin_coors, cluster_coors


def checkNoe(atom1_coor, atom2_coor, noedef):
    noe_bool = False
    dist = 999.999
    error_array = []
    dist_array = []
    for i in range(0, len(atom1_coor[0])):
        x1, y1, z1 = atom1_coor[0][i], atom1_coor[1][i], atom1_coor[2][i]
        for j in range(0, len(atom2_coor[0])):
            x2, y2, z2 = atom2_coor[0][j], atom2_coor[1][j], atom2_coor[2][j]
            dist = math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1))
            dist = round(dist, 2)

            if dist <= 1.8:
                dist_array.append(999.999)
                error_array.append(999.999)

            elif (dist <= (noedef[4] + noedef[6])) and (dist >= (noedef[4] - noedef[5])):
                noe_bool = True
                return noe_bool, dist, 0.0, 0.0
            else:
                if (dist > (noedef[4] + noedef[6])):
                    error = (dist - (noedef[4] + noedef[6])) * (dist - (noedef[4] + noedef[6]))
                    error_array.append(error)
                    dist_array.append(dist)

                else:
                    error = (dist - (noedef[4] - noedef[5])) * (dist - (noedef[4] - noedef[5]))
                    error_array.append(error)
                    dist_array.append(dist)
    dist_array.sort()
    error_array.sort()
    return noe_bool, dist, dist_array[0], error_array[0]


def getAtomCoors(noedef, coorH_matrix, bb_matrix, cluster_protons, cluster_sidechains, resi):
    """

    :param noedef:
    :param coorH_matrix:
    :param bb_matrix:
    :param cluster_protons:
    :param cluster_sidechains:
    :param resi:
    :return:
    """
    # noedef = [31, 'H', 47, 'HG1', 4.2, 2.4, 0.63, 'E', 'V']
    methyls = {'I': ['HG2', 'HD1'], 'L': ['HD1', 'HD2', 'HD'], 'V': ['HG1', 'HG2', 'HG'], 'A': ['HB']}
    cluster_proton_entries = cluster_protons.keys()
    atom1_coor, atom2_coor = [], []

    if noedef[1] == 'H': # if just protons
        atom1_coor = coorH_matrix[noedef[0]]

    elif noedef[1] in methyls[noedef[7]]: # having methyl entries
        if (noedef[0], noedef[1]) in cluster_proton_entries:
            atom1_coor = cluster_protons[(noedef[0], noedef[1])]
        else:
            bb1_coors = bb_matrix[noedef[0]]
            atom1_coor, sidechains1 = getILVARotamers(noedef[7], bb1_coors, noedef[1])
            if atom1_coor:
                cluster_protons[(noedef[0], noedef[1])] = atom1_coor
                cluster_sidechains[(noedef[0], noedef[1])] = sidechains1
    else:
        print "Not a proton spin or missing residue number"

    if noedef[3] == 'H':
        atom2_coor = coorH_matrix[noedef[2]]
    elif noedef[3] in methyls[noedef[8]]:
        if (noedef[2], noedef[3]) in cluster_proton_entries:
            atom2_coor = cluster_protons[(noedef[2], noedef[3])]
        else:
            bb2_coors = bb_matrix[noedef[2]]
            atom2_coor, sidechains2 = getILVARotamers(noedef[8], bb2_coors, noedef[3])
            if atom2_coor:
                cluster_protons[(noedef[2], noedef[3])] = atom2_coor
                cluster_sidechains[(noedef[2], noedef[3])] = sidechains2
    else:
        print "Not a proton spin or missing residue number"

    return atom1_coor, atom2_coor, cluster_protons, cluster_sidechains


def noeinresi(noe, resi):
    if len(noe) == 9:
        if (noe[0] in resi) and (noe[2] in resi):
            return True
        else:
            return False
    else:
        tnoe = noe[0]
        if (tnoe[0] in resi) and (tnoe[2] in resi):
            return True
        else:
            return False


def noeinresi2(noe, resi):
    if len(noe) == 9:
        if (noe[0] in resi[0]) and (noe[2] in resi[1]):
            return True
        elif (noe[0] in resi[1]) and (noe[2] in resi[0]):
            return True
        else:
            return False
    else:
        tnoe = noe[0]
        if (tnoe[0] in resi[0]) and (tnoe[2] in resi[1]):
            return True
        elif (tnoe[0] in resi[1]) and (tnoe[2] in resi[0]):
            return True
        else:
            return False


def extractUnSatisfiedNoes(noes_found, impossible_noes, data):
    unsatisfied = []
    for entry in data:
        if (entry in noes_found) or (entry in impossible_noes):
            pass
        else:
            unsatisfied.append(entry)
    return unsatisfied


def getNOEenergy(error_array, total_noes):
    noe_energy = numpy.sum(error_array)
    noe_energy = noe_energy / float(total_noes)
    return noe_energy / math.pow(total_noes, 1 / 3.0)


def s1ILVApdf(s1_def, s2_def, smotif, exp_data, stage):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param exp_data:
    :return:
    """

    # methyls = {'I': ['HG2', 'HD1'], 'L': ['HD1', 'HD2', 'HD'], 'V': ['HG1', 'HG2', 'HG'], 'A': ['HB']}
    noe_data = exp_data['ilva_noes']
    cluster_protons = {}      # What does it carry
    cluster_sidechains = {}   # what does it carry

    satisfied_noes = []
    unsatisfied_noes = []
    impossible_noes = []
    error_array = []

    max_noe_limit = exp_data['max_noe_dist']
    max_violations = exp_data['max_violations']
    noe_energy_cutoff = exp_data['noe_energy_cutoff']

    noes_found = 0.0
    total_noes = 0.0
    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)
    resi = ss1_list + ss2_list
    resi.sort()

    coorH_matrix = getHCoorMatrix(ss1_list, smotif[1])
    coor2_matrix = getHCoorMatrix(ss2_list, smotif[2])
    coorH_matrix.update(coor2_matrix)

    bb_matrix = getBackboneCoors(ss1_list, smotif[1])
    bb2_matrix = getBackboneCoors(ss2_list, smotif[2])
    bb_matrix.update(bb2_matrix)

    smotif_noe_data = []
    for entry in noe_data:
        if noeinresi(entry, resi):
            smotif_noe_data.append(entry)

    count = 0.0
    tol_noe_count = 0
    for noedef in smotif_noe_data:

        if len(noedef) == 9:
            # noedef = [31, 'H', 47, 'HG1', 4.2, 2.4, 0.63, 'E', 'V']
            atom1_coor, atom2_coor, cluster_protons, cluster_sidechains = getAtomCoors(noedef, coorH_matrix, bb_matrix,
                                                                                       cluster_protons,
                                                                                       cluster_sidechains, resi)
            if atom1_coor and atom2_coor:
                noe_bool, dist, lowest_dist, error = checkNoe(atom1_coor, atom2_coor, noedef)
                error_array.append(error)
                satisfied_noes.append(noedef)
                total_noes += 1.0
                noes_found += 1.0
                if lowest_dist > max_noe_limit:
                    tol_noe_count += 1
                    if tol_noe_count > max_violations:
                        return 0.0, noes_found, 0.00, [satisfied_noes,
                                                       unsatisfied_noes], cluster_protons, cluster_sidechains
            else:
                impossible_noes.append(noedef)
                total_noes += 1
        else:
            dist = 999.999
            lowest_dist = 0.0
            error = 0.0
            noe_bool = False
            for noe in noedef:
                atom1_coor, atom2_coor, cluster_protons, cluster_sidechains = getAtomCoors(noe, coorH_matrix, bb_matrix,
                                                                                           cluster_protons,
                                                                                           cluster_sidechains, resi)
                if atom1_coor and atom2_coor:
                    noe_bool, dist, lowest_dist, error = checkNoe(atom1_coor, atom2_coor, noe)
                    if noe_bool:
                        break

            if noe_bool:
                error_array.append(error)
                satisfied_noes.append(noedef)
                noes_found += 1
                total_noes += 1
            elif error or lowest_dist:
                error_array.append(error)
                satisfied_noes.append(noedef)
                noes_found += 1
                total_noes += 1
                if lowest_dist > max_noe_limit:
                    tol_noe_count += 1
                    if tol_noe_count > max_violations:
                        return 0.0, noes_found, 0.00, [satisfied_noes,
                                                       unsatisfied_noes], cluster_protons, cluster_sidechains
            elif dist == 999.999:
                impossible_noes.append(noedef)
                total_noes += 1
            else:
                print "1:WTH did i miss"

    noe_energy = getNOEenergy(error_array, total_noes)
    probability = noes_found / total_noes
    threshold = exp_data['expected_noe_prob'][stage - 1]

    if probability < threshold:
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains
    elif noe_energy > noe_energy_cutoff:
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains
    elif (noe_energy <= noe_energy_cutoff) and (probability >= threshold):
        unsatisfied_noes = extractUnSatisfiedNoes(satisfied_noes, impossible_noes, noe_data)
        return probability, total_noes, noe_energy, [satisfied_noes, impossible_noes, unsatisfied_noes,
                                                     error_array], cluster_protons, cluster_sidechains
    else:
        print "Fatal error in logic in s1ILVApdf"
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains

def getSxCoorMatrix(coor_array, native_sse):
    resi = {}
    native_sse_range = range(native_sse[4], native_sse[5] + 1)
    count_array = 0
    for i in range(1, len(coor_array[0]), 5):
        resi[native_sse_range[count_array]] = [[coor_array[0][i]], [coor_array[1][i]], [coor_array[2][i]], ['H']]
        count_array += 1
    return resi


def getSxbbCoorMatrix(coor_array, native_sse):
    """
    Double checked this def works
    :param coor_array:
    :param native_sse:
    :return:
    """
    resi = {}
    native_sse_range = range(native_sse[4], native_sse[5] + 1)
    count_array = 0

    for i in range(0, len(coor_array[0]), 5):
        x, y, z = [None] * 5, [None] * 5, [None] * 5
        for j in range(0, 5):
            x[j] = coor_array[0][i + j]
            y[j] = coor_array[1][i + j]
            z[j] = coor_array[2][i + j]
        resi[native_sse_range[count_array]] = [x, y, z]
        count_array += 1
    return resi


def getSxILVARotamers(res_type, bbc, spin):
    import glob, os
    cwd = (os.path.dirname(os.path.realpath(__file__)))
    file_name = cwd + '/sidechainRotamers/' + res_type + "_sc.pickle"
    rotamers = io.readPickle(file_name)
    rmsd_cutoff = 0.1
    spin_coors, cluster_coors = bbRMSD.bbrmsd(bbc, rotamers, rmsd_cutoff, spin, res_type)
    return spin_coors, cluster_coors


def getSxAtomCoors(noedef, coorH_matrix, bb_matrix, cluster_protons, cluster_sidechains, resi):
    methyls = {'I': ['HG2', 'HD1'], 'L': ['HD1', 'HD2', 'HD'], 'V': ['HG1', 'HG2', 'HG'], 'A': ['HB']}

    cluster_proton_entries = cluster_protons.keys()
    atom1_coor, atom2_coor = [], []

    if noedef[1] == 'H':
        atom1_coor = coorH_matrix[noedef[0]]
    elif noedef[1] in methyls[noedef[7]]:
        if (noedef[0], noedef[1]) in cluster_proton_entries:
            atom1_coor = cluster_protons[(noedef[0], noedef[1])]
        else:
            bb1_coors = bb_matrix[noedef[0]]
            atom1_coor, sidechains1 = getSxILVARotamers(noedef[7], bb1_coors, noedef[1])
            if atom1_coor:
                cluster_protons[(noedef[0], noedef[1])] = atom1_coor
                cluster_sidechains[(noedef[0], noedef[1])] = sidechains1
    else:
        print "Not a proton spin or missing residue number"

    if noedef[3] == 'H':
        atom2_coor = coorH_matrix[noedef[2]]
    elif noedef[3] in methyls[noedef[8]]:
        if (noedef[2], noedef[3]) in cluster_proton_entries:
            atom2_coor = cluster_protons[(noedef[2], noedef[3])]
        else:
            bb2_coors = bb_matrix[noedef[2]]
            atom2_coor, sidechains2 = getSxILVARotamers(noedef[8], bb2_coors, noedef[3])
            if atom2_coor:
                cluster_protons[(noedef[2], noedef[3])] = atom2_coor
                cluster_sidechains[(noedef[2], noedef[3])] = sidechains2
    else:
        print "Not a proton spin or missing residue number"

    return atom1_coor, atom2_coor, cluster_protons, cluster_sidechains


def sX2ILVApdf(transformed_coors, native_sse_order, sorted_noe_data, cluster_protons, cluster_sidechains,
               exp_data, stage):

    sse_coors = copy.deepcopy(transformed_coors)
    satisfied_noes = copy.deepcopy(sorted_noe_data[0])
    impossible_noes = copy.deepcopy(sorted_noe_data[1])
    noe_data = copy.deepcopy(sorted_noe_data[2])  # this array only contains the NOEs that are never computed before!
    error_array = copy.deepcopy(sorted_noe_data[3])

    max_noe_limit = exp_data['max_noe_dist']
    max_violations = exp_data['max_violations']
    noe_energy_cutoff = exp_data['noe_energy_cutoff']
    unsatisfied_noes = []

    noes_found = len(satisfied_noes)
    total_noes = len(satisfied_noes) + len(impossible_noes)

    coorH_matrix = {}
    for i in range(0, len(sse_coors)):
        tcoor = getSxCoorMatrix(sse_coors[i], native_sse_order[i])
        coorH_matrix.update(tcoor)

    resi = coorH_matrix.keys()

    bb_matrix = {}
    for i in range(0, len(sse_coors)):
        tcoor = getSxbbCoorMatrix(sse_coors[i], native_sse_order[i])
        bb_matrix.update(tcoor)

    smotif_noe_data = []
    for entry in noe_data:  # this only extracts NOEs that have resi's for the current set
        if noeinresi(entry, resi):
            smotif_noe_data.append(entry)

    tol_noe_count = 0

    for noedef in smotif_noe_data:

        if len(noedef) == 9:
            atom1_coor, atom2_coor, cluster_protons, cluster_sidechains = getSxAtomCoors(noedef, coorH_matrix,
                                                                                         bb_matrix, cluster_protons,
                                                                                         cluster_sidechains, resi)
            if atom1_coor and atom2_coor:
                noe_bool, dist, lowest_dist, error = checkNoe(atom1_coor, atom2_coor, noedef)
                error_array.append(error)
                satisfied_noes.append(noedef)
                total_noes += 1.0
                noes_found += 1.0
                if lowest_dist > max_noe_limit:
                    tol_noe_count += 1
                    if tol_noe_count > max_violations:
                        return 0.0, noes_found, 0.00, [satisfied_noes,
                                                       unsatisfied_noes], cluster_protons, cluster_sidechains
            else:
                impossible_noes.append(noedef)
                total_noes += 1.0
        else:
            dist = 999.999
            lowest_dist = 0.0
            error = 0.0
            noe_bool = False
            for noe in noedef:
                atom1_coor, atom2_coor, cluster_protons, cluster_sidechains = getSxAtomCoors(noe, coorH_matrix,
                                                                                             bb_matrix, cluster_protons,
                                                                                             cluster_sidechains, resi)
                if atom1_coor and atom2_coor:
                    noe_bool, dist, lowest_dist, error = checkNoe(atom1_coor, atom2_coor, noe)
                    if noe_bool:
                        break
            if noe_bool:
                error_array.append(error)
                satisfied_noes.append(noedef)
                noes_found += 1.0
                total_noes += 1.0
            elif error or lowest_dist:
                error_array.append(error)
                satisfied_noes.append(noedef)
                noes_found += 1.0
                total_noes += 1.0
                if lowest_dist > max_noe_limit:
                    tol_noe_count += 1
                    if tol_noe_count > max_violations:
                        return 0.0, noes_found, 0.00, [satisfied_noes, unsatisfied_noes,
                                                       error_array], cluster_protons, cluster_sidechains
            elif dist == 999.999:
                impossible_noes.append(noedef)
                total_noes += 1
            else:
                print "2+:WTH did i miss"

    noe_energy = getNOEenergy(error_array, total_noes)
    probability = noes_found / total_noes
    threshold = exp_data['expected_noe_prob'][stage - 1]
    if probability < threshold:
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains
    elif noe_energy > noe_energy_cutoff:
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains
    elif (noe_energy <= noe_energy_cutoff) and (probability >= threshold):
        unsatisfied_noes = extractUnSatisfiedNoes(satisfied_noes, impossible_noes, noe_data)
        return probability, total_noes, noe_energy, [satisfied_noes, impossible_noes, unsatisfied_noes,
                                                     error_array], cluster_protons, cluster_sidechains
    else:
        print "Fatal error in logic in sXILVApdf"
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains


def noe_in_pairOLD(tsse_ordered, exp_data, pair):
    resi = []
    sse_ordered = [tsse_ordered[pair[0]], tsse_ordered[pair[1]]]
    for i in range(0, 2):
        native_sse_range = range(sse_ordered[i][4], sse_ordered[i][5] + 1)
        resi = resi + native_sse_range
    noe_data = exp_data['ilva_noes']
    smotif_noe_data = []
    print "noe in pair:", resi
    for entry in noe_data:
        if noeinresi(entry, resi):
            smotif_noe_data.append(entry)
            print entry
    if len(smotif_noe_data) > 0:
        return True
    else:
        return False


def noe_in_pair(tsse_ordered, exp_data, pair):
    resi = []
    sse_ordered = [tsse_ordered[pair[0]], tsse_ordered[pair[1]]]
    for i in range(0, 2):
        native_sse_range = range(sse_ordered[i][4], sse_ordered[i][5] + 1)
        resi.append(native_sse_range)
    noe_data = exp_data['ilva_noes']
    smotif_noe_data = []
    for entry in noe_data:
        if noeinresi2(entry, resi):
            smotif_noe_data.append(entry)
    if len(smotif_noe_data) > 0:
        return True
    else:
        return False


def refineILVA(transformed_coors, sse_ordered, exp_data, old_noe_energy):

    #sse_coors = copy.deepcopy(transformed_coors)
    sse_coors = (transformed_coors)[:]
    cluster_protons = {}
    cluster_sidechains = {}
    satisfied_noes = []
    impossible_noes = []
    error_array = []
    noe_data = exp_data['ilva_noes']
    max_noe_limit = exp_data['max_noe_dist']
    max_violations = exp_data['max_violations']
    unsatisfied_noes = []

    noes_found = len(satisfied_noes)
    total_noes = len(satisfied_noes) + len(impossible_noes)

    coorH_matrix = {}
    try:
        for i in range(0, len(sse_coors)):
            tcoor = getSxCoorMatrix(sse_coors[i], sse_ordered[i])
            coorH_matrix.update(tcoor)

        resi = coorH_matrix.keys()

    except IndexError:
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains

    bb_matrix = {}
    for i in range(0, len(sse_coors)):
        tcoor = getSxbbCoorMatrix(sse_coors[i], sse_ordered[i])
        bb_matrix.update(tcoor)

    smotif_noe_data = []
    for entry in noe_data:
        if noeinresi(entry, resi):
            smotif_noe_data.append(entry)
    tol_noe_count = 0

    for noedef in smotif_noe_data:

        if len(noedef) == 9:
            atom1_coor, atom2_coor, cluster_protons, cluster_sidechains = getSxAtomCoors(noedef, coorH_matrix,
                                                                                         bb_matrix, cluster_protons,
                                                                                         cluster_sidechains, resi)
            if atom1_coor and atom2_coor:
                noe_bool, dist, lowest_dist, error = checkNoe(atom1_coor, atom2_coor, noedef)
                error_array.append(error)
                satisfied_noes.append(noedef)
                total_noes += 1.0
                noes_found += 1.0
                if lowest_dist > max_noe_limit:
                    tol_noe_count += 1
                    if tol_noe_count > max_violations:
                        return 0.0, noes_found, 0.00, [satisfied_noes,
                                                       unsatisfied_noes], cluster_protons, cluster_sidechains
            else:
                impossible_noes.append(noedef)
                total_noes += 1.0
        else:
            dist = 999.999
            lowest_dist = 0.0
            error = 0.0
            noe_bool = False
            for noe in noedef:
                atom1_coor, atom2_coor, cluster_protons, cluster_sidechains = getSxAtomCoors(noe, coorH_matrix,
                                                                                             bb_matrix, cluster_protons,
                                                                                             cluster_sidechains, resi)
                if atom1_coor and atom2_coor:
                    noe_bool, dist, lowest_dist, error = checkNoe(atom1_coor, atom2_coor, noe)
                    if noe_bool:
                        break
            if noe_bool:
                error_array.append(error)
                satisfied_noes.append(noedef)
                noes_found += 1.0
                total_noes += 1.0
            elif error or lowest_dist:
                error_array.append(error)
                satisfied_noes.append(noedef)
                noes_found += 1.0
                total_noes += 1.0
                if lowest_dist > max_noe_limit:
                    tol_noe_count += 1
                    if tol_noe_count > max_violations:
                        return 0.0, noes_found, 0.00, [satisfied_noes, unsatisfied_noes,
                                                       error_array], cluster_protons, cluster_sidechains
            elif dist == 999.999:
                impossible_noes.append(noedef)
                total_noes += 1
            else:
                print "2+:WTH did i miss"

    noe_energy = numpy.sum(error_array)
    noe_energy = noe_energy / float(total_noes)
    noe_energy = noe_energy / math.pow(total_noes, 1 / 3.0)

    if noe_energy > old_noe_energy:
        return 0.001, noes_found, 0.00, [satisfied_noes, unsatisfied_noes], cluster_protons, cluster_sidechains
    else:
        unsatisfied_noes = extractUnSatisfiedNoes(satisfied_noes, impossible_noes, noe_data)
        return (noes_found / total_noes), total_noes, noe_energy, [satisfied_noes, impossible_noes, unsatisfied_noes,
                                                                   error_array], cluster_protons, cluster_sidechains

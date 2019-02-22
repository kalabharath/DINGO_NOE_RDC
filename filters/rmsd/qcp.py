#!/usr/bin/env python

"""
Project_Name: rmsd/qcp, File_name: qcp.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 8/05/15 , Time:4:27 AM
"""

import copy

import qcprot


def dumpPDBCoo(coo_array):
    """

    :param coo_array:
    :return:
    """
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]

        pdb_line = "%-6s%5d %2s%1s%s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" \
                   % ('ATOM', i + 1, 'CA', "", 'ALA', 'A', i + 1, " ", x, y, z, 1.0, 30.0, ' ', ' ')
        print pdb_line
    print 'TER'
    return True


def dumpPDBCoo2(coo_array):
    """

    :param coo_array:
    :return:
    """
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        atom = coo_array[3][i]
        res_no = coo_array[4][i]
        res = coo_array[5][i]

        pdb_line = "%-6s%5d  %-2s%5s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" \
                   % ('ATOM', i + 1, atom, res, 'A', res_no, " ", x, y, z, 1.0, 30.0, ' ', ' ')
        print pdb_line
    print 'TER'
    return True


def getCAcoo(frag):

    """
    :param frag:
    :return:
    """
    # print frag
    factor = len(frag[0]) / 5
    x, y, z = [None] * factor, [None] * factor, [None] * factor
    index = 0
    for i in range(2, len(frag[0]), 5):
        x[index] = (frag[0][i])
        y[index] = (frag[1][i])
        z[index] = (frag[2][i])
        index += 1
    return [x, y, z]


def getXcoo(frag, atom_num):
    """
    :param frag:
    :return:
    """
    # print frag
    factor = len(frag[0]) / 5
    x, y, z = [None] * factor, [None] * factor, [None] * factor
    index = 0
    for i in range(atom_num, len(frag[0]), 5):
        x[index] = (frag[0][i])
        y[index] = (frag[1][i])
        z[index] = (frag[2][i])
        index += 1
    return [x, y, z]


def getCcoo(frag):
    """
    :param frag:
    :return:
    """
    # print frag
    factor = len(frag[0]) / 5
    x, y, z = [None] * factor, [None] * factor, [None] * factor
    index = 0
    for i in range(3, len(frag[0]), 5):
        x[index] = (frag[0][i])
        y[index] = (frag[1][i])
        z[index] = (frag[2][i])
        index += 1
    return [x, y, z]


def getcoo(frag):
    """

    :param frag:
    :return:
    """
    # [41, 'ASP', 'N', 28.117, -17.694, 7.215]

    x, y, z, atom_type = [None] * len(frag), [None] * len(frag), [None] * len(frag), [None] * len(frag)
    res_no, res = [None] * len(frag), [None] * len(frag)
    for i in range(0, len(frag)):
        x[i] = frag[i][3]
        y[i] = frag[i][4]
        z[i] = frag[i][5]
        atom_type[i] = frag[i][2]
        res_no[i] = frag[i][0]
        res[i] = frag[i][1]
    return [x, y, z, atom_type, res_no, res]


def centerCoo(coo_array):
    """

    :param coo_array:
    :return:
    """

    xsum, ysum, zsum = 0, 0, 0

    for i in range(0, len(coo_array[0])):
        xsum += coo_array[0][i]
        ysum += coo_array[1][i]
        zsum += coo_array[2][i]

    xsum /= len(coo_array[0])
    ysum /= len(coo_array[0])
    zsum /= len(coo_array[0])

    for i in range(0, len(coo_array[0])):
        coo_array[0][i] -= xsum
        coo_array[1][i] -= ysum
        coo_array[2][i] -= zsum

    return coo_array, [xsum, ysum, zsum]


def applyTranslation(frag, cen_mass):
    """

    :param frag:
    :param cen_mass:
    :return:
    """

    for i in range(0, len(frag[0])):
        frag[0][i] += cen_mass[0]
        frag[1][i] += cen_mass[1]
        frag[2][i] += cen_mass[2]
    return frag


def translateCM(coo_array, cm):
    """

    :param coo_array:
    :param cm:
    :return:
    """
    xsum, ysum, zsum = cm[0], cm[1], cm[2]
    for i in range(0, len(coo_array[0])):
        coo_array[0][i] -= xsum
        coo_array[1][i] -= ysum
        coo_array[2][i] -= zsum
    return coo_array


def applyRot(frag, rotmat):
    """

    :param frag:
    :param rotmat:
    :return:
    """
    for i in range(0, len(frag[0])):
        x = rotmat[0] * frag[0][i] + rotmat[1] * frag[1][i] + rotmat[2] * frag[2][i]
        y = rotmat[3] * frag[0][i] + rotmat[4] * frag[1][i] + rotmat[5] * frag[2][i]
        z = rotmat[6] * frag[0][i] + rotmat[7] * frag[1][i] + rotmat[8] * frag[2][i]
        frag[0][i] = x
        frag[1][i] = y
        frag[2][i] = z
    return frag


def get_dist(r1, r2):
    import math
    x1, y1, z1 = r1[0], r1[1], r1[2]
    x2, y2, z2 = r2[0], r2[1], r2[2]
    return math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1))


def rmsdQCP(psmotif, csmotif, direction, cutoff, previous_sse_index):

    """
    use rmsdQCP
    :param psmotif:
    :param csmotif:
    :param direction:
    :param cutoff:
    :param previous_sse_index:
    :return:
    """

    psmotif = (psmotif[1])[:]
    psmotif_index = previous_sse_index[-1]

    if psmotif_index == 0:
        indexa =1
        indexb = 2
    elif psmotif_index == 1:
        indexa = 2
        indexb = 1
    else:
        print "Error"

    # Extra careful about this if-else loop, the trick is finding the the correct arrays
    if direction == 'left':
        native_fraga = getcoo(psmotif[indexa])
        frag_a = getcoo(psmotif[indexa])
        frag_b = getcoo(csmotif[2])
        native_fragb_2ndsse = (csmotif[1])[:]
        native_fraga_2ndsse = getcoo(psmotif[indexb])
    else:
        native_fraga = getcoo(psmotif[indexa])
        frag_a = getcoo(psmotif[indexa])
        frag_b = getcoo(csmotif[1])
        native_fragb_2ndsse = (csmotif[2])[:]
        native_fraga_2ndsse = getcoo(psmotif[indexb])

    frag_a, a_cen = centerCoo(frag_a)
    frag_b, b_cen = centerCoo(frag_b)

    frag_aca = getCAcoo(frag_a)
    frag_bca = getCAcoo(frag_b)

    fraglen = len(frag_aca[0])

    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, frag_aca[0][i])
        qcprot.SetDArray(1, i, xyz1, frag_aca[1][i])
        qcprot.SetDArray(2, i, xyz1, frag_aca[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, frag_bca[0][i])
        qcprot.SetDArray(1, i, xyz2, frag_bca[1][i])
        qcprot.SetDArray(2, i, xyz2, frag_bca[2][i])

    rot = qcprot.MakeDvector(9)

    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********

    if rmsd > cutoff:
        # exit without further computation
        # free memory
        qcprot.FreeDMatrix(xyz1)
        qcprot.FreeDMatrix(xyz2)
        qcprot.FreeDArray(rot)
        return rmsd, []

    rotmat = [None] * 9
    for i in range(0, 9):
        rotmat[i] = qcprot.GetDvector(i, rot)

    # translate the other SSE of the current smotif
    sse_2nd_coos = getcoo(native_fragb_2ndsse)

    cm_sse2nd = translateCM(sse_2nd_coos, b_cen)
    rot_sse_2nd = applyRot(cm_sse2nd, rotmat)
    trans_sse2nd = applyTranslation(rot_sse_2nd, a_cen)

    # return 3 arrays of coordinates
    # The order problem ?

    if direction == 'left':
        if psmotif_index == 0:
            transformed_coor = [trans_sse2nd, native_fraga, native_fraga_2ndsse]
        else:
            transformed_coor = [trans_sse2nd, native_fraga_2ndsse, native_fraga]
    else:
        if psmotif_index == 1: # double check this
            transformed_coor = [native_fraga_2ndsse, native_fraga, trans_sse2nd]
        else:
            transformed_coor = [native_fraga, native_fraga_2ndsse, trans_sse2nd]

    # free memory
    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    return rmsd, transformed_coor


def rmsdQCP3(previous_smotif, csmotif, direction, cutoff, previous_sse_index):

    """
    :param previous_smotif:
    :param csmotif:
    :param direction:
    :param cutoff:
    :param previous_sse_index:
    :return:
    """
    psmotif = (previous_smotif[2][1])[:]
    psmotif_index = previous_sse_index[-1]

    if direction == 'left':
        frag_b = getcoo(csmotif[2])
        native_fragb_2ndsse = (csmotif[1])[:]
        frag_a = copy.deepcopy(psmotif[psmotif_index])
    else:
        frag_a = copy.deepcopy(psmotif[psmotif_index])
        frag_b = getcoo(csmotif[1])
        native_fragb_2ndsse = (csmotif[2])[:]

    frag_a, a_cen = centerCoo(frag_a)
    frag_b, b_cen = centerCoo(frag_b)

    frag_aca = getCAcoo(frag_a)
    frag_bca = getCAcoo(frag_b)

    fraglen = len(frag_aca[0])

    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, frag_aca[0][i])
        qcprot.SetDArray(1, i, xyz1, frag_aca[1][i])
        qcprot.SetDArray(2, i, xyz1, frag_aca[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, frag_bca[0][i])
        qcprot.SetDArray(1, i, xyz2, frag_bca[1][i])
        qcprot.SetDArray(2, i, xyz2, frag_bca[2][i])

    rot = qcprot.MakeDvector(9)

    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********

    if rmsd > cutoff:
        # exit without further computation
        # free memory
        qcprot.FreeDMatrix(xyz1)
        qcprot.FreeDMatrix(xyz2)
        qcprot.FreeDArray(rot)
        return rmsd, []

    rotmat = [None] * 9
    for i in range(0, 9):
        rotmat[i] = qcprot.GetDvector(i, rot)

    # translate the other SSE of the current smotif

    sse_2nd_coos = getcoo(native_fragb_2ndsse)
    cm_sse2nd = translateCM(sse_2nd_coos, b_cen)
    rot_sse_2nd = applyRot(cm_sse2nd, rotmat)
    trans_sse2nd = applyTranslation(rot_sse_2nd, a_cen)

    # append the translated coordinates
    temp_holder = (psmotif)[:]

    if direction == 'left':
        temp_holder.insert(0, trans_sse2nd)
    else:
        temp_holder.append(trans_sse2nd)

    # free memory
    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    return rmsd, temp_holder


def rmsdQCP4(pair, presse, alt_smotif_log,  csmotif, direction, cutoff):
    """
    TODO Delete this function
    :param pair:
    :param presse:
    :param alt_smotif_log:
    :param csmotif:
    :param direction:
    :param cutoff:
    :return:
    """

    from utility.alt_smotif_util import delete_last_sse

    if direction == 'left':
        frag_b = getcoo(csmotif[2])
        native_fragb_2ndsse = (csmotif[1])[:]
        frag_a = copy.deepcopy(presse[pair[1]])
        presse = delete_last_sse(presse, alt_smotif_log)
    else:
        frag_a = copy.deepcopy(presse[pair[0]])
        frag_b = getcoo(csmotif[1])
        native_fragb_2ndsse = (csmotif[2])[:]
        presse = delete_last_sse(presse, alt_smotif_log)

    frag_a, a_cen = centerCoo(frag_a)
    frag_b, b_cen = centerCoo(frag_b)

    frag_aca = getCAcoo(frag_a)
    frag_bca = getCAcoo(frag_b)

    fraglen = len(frag_aca[0])
    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, frag_aca[0][i])
        qcprot.SetDArray(1, i, xyz1, frag_aca[1][i])
        qcprot.SetDArray(2, i, xyz1, frag_aca[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, frag_bca[0][i])
        qcprot.SetDArray(1, i, xyz2, frag_bca[1][i])
        qcprot.SetDArray(2, i, xyz2, frag_bca[2][i])

    rot = qcprot.MakeDvector(9)

    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********

    if rmsd > cutoff:
        # exit without further computation
        # free memory
        qcprot.FreeDMatrix(xyz1)
        qcprot.FreeDMatrix(xyz2)
        qcprot.FreeDArray(rot)
        return rmsd, []

    rotmat = [None] * 9
    for i in range(0, 9):
        rotmat[i] = qcprot.GetDvector(i, rot)

    # translate the other SSE of the current smotif

    sse_2nd_coos = getcoo(native_fragb_2ndsse)
    cm_sse2nd = translateCM(sse_2nd_coos, b_cen)
    rot_sse_2nd = applyRot(cm_sse2nd, rotmat)
    trans_sse2nd = applyTranslation(rot_sse_2nd, a_cen)

    # append the translated coordinates
    temp_holder = presse[:]

    if direction == 'left':
        temp_holder.insert(0, trans_sse2nd)
    else:
        temp_holder.append(trans_sse2nd)

    # free memory
    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    return rmsd, temp_holder


def getKdist(sse_array, atom_type):
    """

    :param sse_array:
    :param atom_type:
    :return:
    """
    # atom_array = ['N', 'H', 'CA', 'C', 'O']
    smotif_type = sse_array[0][0][0] + sse_array[1][0][0]
    mean, sd = 0.0, 0.0
    dist_knowledge = {'ss': [['N'], [1.9, 0.16], [3.4, 0.07], [3.78, 0.1], ['O']],
                      'hh': [['N'], [2.68, 0.06], [3.8, 0.05], [3.0, 0.03], ['O']],
                      'hs': [['N'], [2.19, 0.06], [3.78, 0.00], [3.1, 0.06], ['O']],
                      'sh': [['N'], [3.16, 0.17], [3.73, 0.1], [2.88, 0.02], ['O']]}
    mean, sd = dist_knowledge[smotif_type][atom_type]
    return mean - (2*sd)


def kClashes(coo_arrays, sse_ordered, current_ss):

    """
    Known minimum distances between various atoms between the SSEs of the smotif
    the first entry is mean and second entry is SD, computed on a sample of 100 lowest distances observed over single
    digit smotifs that have the largest number of entries.
    TODO if doing both NH and CO overlaps alternate the index between the two for 2X speedup

    """

    # [['strand', 15, 10, 3, 59, 73], ['strand', 15, 3, 13, 77, 91], ['strand', 14, 11, 4, 103, 116]]
    # atom_array = ['N', 'H', 'CA', 'C', 'O']

    current_index = sse_ordered.index(current_ss)
    indices = range(0, len(sse_ordered))
    indices.remove(current_index)

    # Compute clashes for amide hydrogen pairs
    atom_type = 1
    sse1 = getXcoo(coo_arrays[current_index], atom_type)
    for p in indices:
        kdist = getKdist([sse_ordered[current_index], sse_ordered[p]], atom_type)
        sse2 = getXcoo(coo_arrays[p], atom_type)
        for j in range(0, len(sse1[0])):
            for k in range(0, len(sse2[0])):
                dist = get_dist([sse1[0][j], sse1[1][j], sse1[2][j]], [sse2[0][k], sse2[1][k], sse2[2][k]])
                dist = round(dist, 2)
                if dist < kdist:
                    return False

    atom_type = 3
    sse1 = getXcoo(coo_arrays[current_index], atom_type)
    for p in indices:
        kdist = getKdist([sse_ordered[current_index], sse_ordered[p]], atom_type)
        sse2 = getXcoo(coo_arrays[p], atom_type)
        for j in range(0, len(sse1[0])):
            for k in range(0, len(sse2[0])):
                dist = get_dist([sse1[0][j], sse1[1][j], sse1[2][j]], [sse2[0][k], sse2[1][k], sse2[2][k]])
                dist = round(dist, 2)
                if dist < kdist:
                    return False
    return True


def kClahsesAltSmotif(coo_arrays, sse_ordered, sse_pairs, alt_def):

    if alt_def[-1] == 'right':
        sse = sse_ordered[sse_pairs[1]]
    else:
        sse = sse_ordered[sse_pairs[0]]

    if kClashes(coo_arrays, sse_ordered, sse):
        return True
    else:
        return False

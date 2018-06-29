#!/usr/bin/env python

"""
Project_Name: main/utility, File_name: smotif_util.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:12:54 PM
"""

import utility.io_util as io


def getSmotif(s1, s2):
    """
    Return Smotif based on SS definition

    ('helix-helix', 6, 23, 4)
    new_smotif_file = hh_6_23.db
    ['helix', 6, 9, 4, 146, 151] ['helix', 23, 4, 1, 156, 178]
    """
    s1_type, s2_type = '', ''

    if s1[0] == 'helix':
        s1_type = 'h'
    if s1[0] == 'strand':
        s1_type = 's'
    if s2[0] == 'helix':
        s2_type = 'h'
    if s2[0] == 'strand':
        s2_type = 's'
    smotif_type = s1_type + s2_type
    s1_len = s1[1]
    s2_len = s2[1]
    smotif = [smotif_type, s1_len, s2_len]
    return smotif


def readSmotifDatabase(smotif, *database_cutoff):
    """

    :param smotif:
    :return:
    """

    import os

    cwd = (os.path.dirname(os.path.realpath(__file__))).split("/")
    root_dir = ''
    for entry in cwd[:-2]:
        if entry == '':
            pass
        else:
            root_dir = root_dir + '/' + entry
    smotif_db_path = root_dir + '/databases/database_cutoff_' + database_cutoff[0] + '/'
    file_name = smotif[0] + "_" + str(smotif[1]) + "_" + str(smotif[2]) + ".db.tar.gz"
    fin = smotif_db_path + file_name
    try:
        smotif_data = io.readTarPickle(fin)
        print "Reading in smotif database: ", fin

        return smotif_data
    except:
        print "Error in reading smotif database: ", fin
        return False


def orderSeq(previous_smotif, current_seq, direction):
    """

    :param previous_smotif:
    :param current_seq:
    :param direction:
    :return:
    """
    # This function is depreciated

    previous_seq = (previous_smotif[4][1])[:]

    if direction == 'left':
        concat_seq = current_seq + previous_seq
    else:
        concat_seq = previous_seq + current_seq

    return concat_seq


def array2string(array):
    string = ''
    for entry in array:
        for t in entry:
            string = string+str(t)
    return string


def orderCATH(previous_smotif, current_smotif, direction, alt_smotf_def):
    """

    :param previous_smotif:
    :param current_smotif:
    :param direction:
    :return:
    """

    previous_cath = (previous_smotif[3][1])[:]
    previous_sse = (previous_smotif[3][2])[:]

    previous_sse.append(alt_smotf_def)
    previous_cath.append(current_smotif)
    """
    if direction == 'left':
        previous_cath.insert(0, current_smotif)
    else:
        previous_cath.append(current_smotif)
    """
    return previous_cath, previous_sse


def getRMSDcutoff(smotif_def):
    rmsd_cutoff = 0

    s1_len, s2_len = smotif_def[1], smotif_def[2]

    if (s1_len < 7) or (s2_len < 7):
        rmsd_cutoff = 1.0
    if ( 6 < s1_len < 16) or (6 < s2_len < 16):
        rmsd_cutoff = 2.0
    if (s1_len > 15) or (s2_len > 15):
        rmsd_cutoff = 2.5
    if (s1_len > 20) or (s2_len > 20):
        rmsd_cutoff = 3.0
    if (s1_len > 30) or (s2_len > 30):
        rmsd_cutoff = 3.5

    return rmsd_cutoff
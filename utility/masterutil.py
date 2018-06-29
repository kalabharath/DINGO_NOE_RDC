import os

import io_util as io
import utility.smotif_util as sm

def enum(*sequential, **named):
    """

    :param sequential:
    :param named:
    :return:
    """
    # fake an enumerated type in Python

    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def getPairSSProfiles(s1, s2, ss_profile):
    """
    return corresponding +/- 2 profile
    definitons for the given given index
    order
    """
    s1_l = ss_profile[s1]
    s2_l = ss_profile[s2]
    return s1_l, s2_l


def getRunSeq():
    """
    generate run seq, a seq list of pairs of
    indexes of profiles for job scheduling
    """
    # TODO needed to generalize based on input data type
    map_route = []

    ss_profiles = io.getSSprofilesFile()
    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")

    s1, s2 = map_route[0][0], map_route[0][1]
    s1_list, s2_list = getPairSSProfiles(s1, s2, ss_profiles)
    run_seq = []
    for i in range(len(s1_list)):
        for j in range(len(s2_list)):
            run_seq.append([i, j])
    return run_seq


def getSSlist():
    map_route = []
    ss_profiles = io.readPickle("ss_profiles.pickle")
    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")

    if map_route:
        s1, s2 = map_route[0][0], map_route[0][1]
        s1_list, s2_list = getPairSSProfiles(s1, s2, ss_profiles)
        return s1_list, s2_list, [s1, s2]
    else:
        return False, False, False


def getSSdef(index_array):
    """

    :param index_array:
    :return:
    """
    s1_list, s2_list, sse_route = getSSlist()
    return s1_list[index_array[0]], s2_list[index_array[1]], sse_route


def getfromDB(previous_smotif, current_ss, direction, database_cutoff, stage, alt_smotif_def):

    """

    :param previous_smotif:
    :param current_ss:
    :param direction:
    :param database_cutoff:
    :param stage:
    :param alt_smotif_def:
    :return:
    """
    if stage == 2:

        previous_sse_index = previous_smotif[0][2]
        psmotif = previous_smotif[1][-1]

        if direction == 'left':
            previous_sse = alt_smotif_def[1]
            previous_ss_index = previous_sse_index.index(previous_sse)
            previous_ss = psmotif[previous_ss_index]
        else:
            previous_sse = alt_smotif_def[0]
            previous_ss_index = previous_sse_index.index(previous_sse)
            previous_ss = psmotif[previous_ss_index]
    else:
        searched_smotifs = (previous_smotif[1][1])[:]
        previous_sse_index = (previous_smotif[1][2])[:]

        if direction == 'left':
            previous_ss_def = previous_sse_index.index(alt_smotif_def[1])
        else:
            previous_ss_def = previous_sse_index.index(alt_smotif_def[0])
        previous_ss = searched_smotifs[previous_ss_def]
        # print "Get correct db:", searched_smotifs, previous_sse_index
        # print alt_smotif_def, previous_ss_def

    # current_ss, previous_ss
    if direction == 'left':  # double check this implementation
        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def, database_cutoff), smotif_def


def generate_refinement_order(sse_array):

    import itertools
    from utility.smotif_util import array2string
    computed_pairs = []
    indices = list(itertools.combinations(range(len(sse_array)), 2))
    refine_pairs = []
    for pair in indices:
        if abs(pair[0] - pair[1]) == 1:
            pass
        else:
            t_array = array2string([sse_array[pair[0]], sse_array[pair[1]]])
            computed_pairs.append(t_array)
            refine_pairs.append(pair)

    return refine_pairs, computed_pairs


def generate_refinement_order2(sse_array, computed_pairs):
    import itertools
    from utility.smotif_util import array2string

    indices = list(itertools.combinations(range(len(sse_array)), 2))
    refine_pairs = []
    for pair in indices:
        if abs(pair[0] - pair[1]) == 1:
            pass
        else:
            t_array = array2string([sse_array[pair[0]], sse_array[pair[1]]])

            if t_array in computed_pairs:
                pass
            else:
                computed_pairs.append(t_array)
                refine_pairs.append(pair)
    return refine_pairs, computed_pairs


def orderSSE(previous_smotif, current_sse, direction, stage, alt_smotif_def):

    """

    :param previous_smotif:
    :param current_sse:
    :param direction:
    :param stage:
    :param alt_smotif_def:
    :return:
    """
    if stage == 2:
        previous_sse_route = (previous_smotif[0][2])[:]
        order_sse_route = previous_sse_route[:]
        ordered_sse = (previous_smotif[1][1])[:]

        if direction == 'left':
            ordered_sse.insert(0, current_sse)
            order_sse_route.insert(0, alt_smotif_def[0])
            previous_sse_index = [direction, previous_sse_route.index(alt_smotif_def[1])]
        else:
            ordered_sse.append(current_sse)
            order_sse_route.append(alt_smotif_def[1])
            previous_sse_index = [direction, previous_sse_route.index(alt_smotif_def[0])]
        return ordered_sse, order_sse_route, previous_sse_route, previous_sse_index

    elif stage > 2:
        previous_sse_route = (previous_smotif[1][2])[:]
        order_sse_route = previous_sse_route[:]
        previous_sse = previous_smotif[2][2]

        if direction == 'left':
            previous_sse.insert(0, current_sse)
            order_sse_route.insert(0, alt_smotif_def[0])
            previous_sse_index = [direction, previous_sse_route.index(alt_smotif_def[1])]
        else:
            previous_sse.append(current_sse)
            order_sse_route.append(alt_smotif_def[1])
            previous_sse_index = [direction, previous_sse_route.index(alt_smotif_def[0])]

        return previous_sse, order_sse_route, previous_sse_route, previous_sse_index

    else:
        computed_pairs = previous_smotif[8][2]
        log_refine_smotif = previous_smotif[8][3]
        previous_sse = previous_smotif[2][2]

        if direction == 'left':
            previous_sse.insert(0, current_sse)
        else:
            previous_sse.append(current_sse)
        refine_pairs, computed_pairs = generate_refinement_order2(previous_sse, computed_pairs)

        return previous_sse, refine_pairs, computed_pairs, log_refine_smotif


def fetchNOEdata(previous_smotif):

    return previous_smotif[5][-3], previous_smotif[5][-2], previous_smotif[5][-1]

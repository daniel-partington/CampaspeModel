from operator import itemgetter
from itertools import groupby

first_entry = lambda x: x[0]
last_entry = lambda x: x[-1]
flatten = lambda l: [item for sublist in l for item in sublist]

import sys
if sys.version_info[0] == 2:
    range = xrange


def merging_reaches(river_seg2, river_seg_temp, rchlen_sum, merge_group, rchlen_weights):

    def weighted(col):
        return (col * rchlen_weights).sum()

    river_seg2.loc[merge_group[0], 'strtop'] = weighted(river_seg_temp['strtop'])
    river_seg2.loc[merge_group[0], 'rchlen'] = rchlen_sum
    river_seg2.set_value(merge_group[0], 'amalg_riv_points', first_entry(river_seg_temp['amalg_riv_points'].tolist()))
    river_seg2.loc[merge_group[0], 'Cumulative Length'] = last_entry(river_seg_temp['Cumulative Length'].tolist())
    river_seg2.loc[merge_group[0], 'strtop_raw'] = weighted(river_seg_temp['strtop_raw'])
    river_seg2.loc[merge_group[0], 'slope'] = weighted(river_seg_temp['slope'])
    river_seg2.loc[merge_group[0], 'k'] = first_entry(river_seg_temp['k'].tolist())
    river_seg2.loc[merge_group[0], 'i'] = first_entry(river_seg_temp['i'].tolist())
    river_seg2.loc[merge_group[0], 'j'] = first_entry(river_seg_temp['j'].tolist())
    river_seg2.set_value(merge_group[0], 'amalg_riv_points_collection', flatten(river_seg_temp['amalg_riv_points_collection']))
    river_seg2.loc[merge_group[0], 'strhc1'] = weighted(river_seg_temp['strhc1'])
    river_seg2.loc[merge_group[0], 'strthick'] = weighted(river_seg_temp['strthick'])
    river_seg2.set_value(merge_group[0], 'amalg_riv_points_tuple', first_entry(river_seg_temp['amalg_riv_points_tuple'].tolist()))

    river_seg2.drop(merge_group[1:], inplace=True)

    return river_seg2



def merge_collocated_stream_reaches(river_segment, max_length=3000.):
    """
    Sort out collocated stream reaches to avoid short circuiting:

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # THINGS THAT NEED TO BE DONE FOR EACH COLUMN OF river_seg WHEN MERGING ROW
    # BASED ON THE merge_group:
    #
    #  strtop = average weighted by rchlen
    #  rchlen = sum()
    #  amalg_riv_points = first entry
    #  Cumulative length = last entry
    #  strtop_raw = average weighted by rchlen
    #  slope = average weighted by rchlen
    #  k = first entry
    #  i = first entry
    #  j = first entry
    #  amalg_riv_points_collection = join lists of tuples into one
    #  strhc1 = average weighted by rchlen
    #  strthick = average weighted by rchlen
    #  amalg_riv_points_tuple = first entry
    #
    #  1. Create new row based on above rules
    #  2. Replace first row indexed in merge_group with new row
    #  3. Delete all other rows indexed in merge_group
    #
    """
    river_seg2 = river_segment.copy()

    collocation_tests = ['skipping a cell', 'adjacent']
    for test in collocation_tests:

        merge_row = []
        for ind in range(river_seg2.shape[0]):
            if ind == 0:
                continue
            elif ind == river_seg2.shape[0] - 1:
                prev = river_seg2.iloc[ind - 1]
                curr = river_seg2.iloc[ind]
            else:
                prev = river_seg2.iloc[ind - 1]
                curr = river_seg2.iloc[ind]
                nexx = river_seg2.iloc[ind + 1]
                #def loc_tup(row):
                #    return (row['i'], row['j'])
                if test == 'skipping a cell':
                    if prev['amalg_riv_points_tuple'] == nexx['amalg_riv_points_tuple']:
                        if curr['rchlen'] < max_length:
                            merge_row += [ind]
                elif test == 'adjacent':
                    if prev['amalg_riv_points_tuple'] == curr['amalg_riv_points_tuple']:
                        if curr['rchlen'] < max_length:
                            merge_row += [ind]

        merge_row_consec = []
        for k, g in groupby(enumerate(merge_row), lambda i_x:i_x[0] - i_x[1]):
            merge_row_consec.append(list(map(itemgetter(1), g)))

        for merge_group in merge_row_consec:
            index_list = river_seg2.index.tolist()
            index_dict = {x:index for index, x in enumerate(index_list)}
            if test == 'skipping a cell':
                merge_group = [index_list[index_dict[merge_group[0]] - 1]] + merge_group
                merge_group = merge_group + [index_list[index_dict[merge_group[-1]] + 1]]
            elif test == 'adjacent':
                merge_group = [index_list[index_dict[merge_group[0]] - 1]] + merge_group

            #merge_group = merge_group + [merge_group[-1] + 1]
            river_seg_temp = river_seg2.loc[merge_group]
            rchlen_temp = river_seg_temp['rchlen']
            rchlen_sum = rchlen_temp.sum()
            rchlen_weights = rchlen_temp / rchlen_sum

            # Is this is the same as `merging_reaches()` defined above?
            def weighted(col):
                return (col * rchlen_weights).sum()

            river_seg2.loc[merge_group[0], 'strtop'] = weighted(river_seg_temp['strtop'])
            river_seg2.loc[merge_group[0], 'rchlen'] = rchlen_sum
            river_seg2.set_value(merge_group[0], 'amalg_riv_points', first_entry(river_seg_temp['amalg_riv_points'].tolist()))
            river_seg2.loc[merge_group[0], 'Cumulative Length'] = last_entry(river_seg_temp['Cumulative Length'].tolist())
            river_seg2.loc[merge_group[0], 'strtop_raw'] = weighted(river_seg_temp['strtop_raw'])
            river_seg2.loc[merge_group[0], 'slope'] = weighted(river_seg_temp['slope'])
            river_seg2.loc[merge_group[0], 'k'] = first_entry(river_seg_temp['k'].tolist())
            river_seg2.loc[merge_group[0], 'i'] = first_entry(river_seg_temp['i'].tolist())
            river_seg2.loc[merge_group[0], 'j'] = first_entry(river_seg_temp['j'].tolist())
            river_seg2.set_value(merge_group[0], 'amalg_riv_points_collection', flatten(river_seg_temp['amalg_riv_points_collection']))
            river_seg2.loc[merge_group[0], 'strhc1'] = weighted(river_seg_temp['strhc1'])
            river_seg2.loc[merge_group[0], 'strthick'] = weighted(river_seg_temp['strthick'])
            river_seg2.set_value(merge_group[0], 'amalg_riv_points_tuple', first_entry(river_seg_temp['amalg_riv_points_tuple'].tolist()))

            river_seg2.drop(merge_group[1:], inplace=True)

        river_seg2.index = list(range(river_seg2.shape[0]))

    return river_seg2


def merge_very_short_stream_reaches(river_segment, min_length=500.):
    """

    """
    river_seg2 = river_segment.copy()

    min_length = min_length
    merge_row_too_short = []
    for ind in range(river_seg2.shape[0]):
        if ind == 0:
            continue
        elif ind == river_seg2.shape[0] - 1:
            prev = river_seg2.iloc[ind - 1]
            curr = river_seg2.iloc[ind]
        else:
            prev = river_seg2.iloc[ind - 1]
            curr = river_seg2.iloc[ind]
            nexx = river_seg2.iloc[ind + 1]
            if prev['amalg_riv_points_tuple'] == nexx['amalg_riv_points_tuple']:
                pass
            else:
                if curr['rchlen'] < min_length:
                    merge_row_too_short += [ind]

    merge_row_too_short_consec = []
    for k, g in groupby(enumerate(merge_row_too_short), lambda i_x1:i_x1[0]-i_x1[1]):
        merge_row_too_short_consec.append(list(map(itemgetter(1), g)))

    for merge_group in merge_row_too_short_consec:
        index_list = river_seg2.index.tolist()
        index_dict = {x:index for index, x in enumerate(index_list)}
        merge_group = [index_list[index_dict[merge_group[0]] - 1]] + merge_group
        #merge_group = merge_group + [merge_group[-1] + 1]
        river_seg_temp = river_seg2.loc[merge_group]
        rchlen_temp = river_seg_temp['rchlen']
        rchlen_sum = rchlen_temp.sum()
        rchlen_weights = rchlen_temp / rchlen_sum
        def weighted(col):
            return (col * rchlen_weights).sum()

        river_seg2.loc[merge_group[0], 'strtop'] = weighted(river_seg_temp['strtop'])
        river_seg2.loc[merge_group[0], 'rchlen'] = rchlen_sum
        river_seg2.set_value(merge_group[0], 'amalg_riv_points', first_entry(river_seg_temp['amalg_riv_points'].tolist()))
        river_seg2.loc[merge_group[0], 'Cumulative Length'] = last_entry(river_seg_temp['Cumulative Length'].tolist())
        river_seg2.loc[merge_group[0], 'strtop_raw'] = weighted(river_seg_temp['strtop_raw'])
        river_seg2.loc[merge_group[0], 'slope'] = weighted(river_seg_temp['slope'])
        river_seg2.loc[merge_group[0], 'k'] = first_entry(river_seg_temp['k'].tolist())
        river_seg2.loc[merge_group[0], 'i'] = first_entry(river_seg_temp['i'].tolist())
        river_seg2.loc[merge_group[0], 'j'] = first_entry(river_seg_temp['j'].tolist())
        river_seg2.set_value(merge_group[0], 'amalg_riv_points_collection', flatten(river_seg_temp['amalg_riv_points_collection']))
        river_seg2.loc[merge_group[0], 'strhc1'] = weighted(river_seg_temp['strhc1'])
        river_seg2.loc[merge_group[0], 'strthick'] = weighted(river_seg_temp['strthick'])
        river_seg2.set_value(merge_group[0], 'amalg_riv_points_tuple', first_entry(river_seg_temp['amalg_riv_points_tuple'].tolist()))

        river_seg2.drop(merge_group[1], inplace=True)

    return river_seg2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

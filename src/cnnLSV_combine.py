def parse_info(seq):
    info = {'SVTYPE': '', 'SVLEN': 0, 'END': 0}
    for i in seq.split(';'):  # SVTYPE=INV;END=139008588;SUPPORT=13;STD_SPAN=4.4;STD_POS=1.5
        if i.split('=')[0] in ['SVLEN', 'END']:
            try:
                info[i.split('=')[0]] = abs(int(i.split('=')[1]))
            except:
                pass
        if i.split('=')[0] == 'SVTYPE':
            info[i.split('=')[0]] = i.split('=')[1][0:3]
    return info

def evaluate(sv1, sv2, bias):
    start1, end1 = sv1[0], sv1[1]
    start2, end2 = sv2[0], sv2[1]
    overlap_len = min(end1, end2) - max(start1, start2)
    # if start1 - offect <= start2 <= end1 + offect or start1 - offect <= end2 <= end1 + offect or start2 - offect <= start1 <= end2 + offect:
    #     if min(end2 - start2, end1 - start1) * 1.0 / max(end2 - start2, end1 - start1) >= bias:
    if overlap_len / max(end1-start1, end2-start2) >= bias:
        return True
    return False

def combine(input, callers, bias):
    index_info = list()
    sv_temp = dict()  # {'del': {'chr' :[[bp1, bp2, index, caller]]}
    f = open(input).readlines()
    if callers:
        rank = dict()
        for i in range(len(callers)):
            rank[callers[i].lower()] = i
    for i, line in enumerate(f):
        if '#' in line:
            continue
        seq = line.strip().split('\t')
        chrname = seq[0]
        pos = int(seq[1])
        caller = seq[2].split('.')[0].lower()
        info = parse_info(seq[7])
        svtype, svlen, end = info['SVTYPE'], info['SVLEN'], info['END']
        if svlen == 0:
            svlen = end - pos
        if svtype not in sv_temp:
            sv_temp[svtype] = {}
        if chrname not in sv_temp[svtype]:
            sv_temp[svtype][chrname] = [[pos, pos + svlen, i, caller]]
        else:
            sv_temp[svtype][chrname].append([pos, pos + svlen, i, caller])

    for st, chr_list in sv_temp.items():
        for chrname, sv_list in sv_temp[st].items():  # 'chr': [pos, end, i, caller], [pos, end, i, caller]
            if len(sv_list) == 1:
                index_info.append([[sv_list[0][2], sv_list[0][3]]])
            else:
                last = sv_list[0]
                index_info.append([[last[2], last[3]]])
                for i in range(1, len(sv_list)):
                    # start, end = sv_list[i][0], sv_list[i][1]
                    # last_start, last_end = last[0], last[1]
                    if evaluate(last, sv_list[i], bias):
                        index_info[-1].append([sv_list[i][2], sv_list[i][3]])
                    else:
                        index_info.append([[sv_list[i][2], sv_list[i][3]]])
                    last = sv_list[i]
    index_list = list()
    for index in index_info:
        if callers:
            index = sorted(index, key=lambda x: rank[x[1]])
        index_list.append([line[0] for line in index])
        # index_list.append(index[0][0])
    return index_list

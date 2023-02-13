import cv2
import numpy as np
import cigar
import sys

def deal_sup(sup_info, min_mapq):
    # [[chr20,167320,+,2436S1758M55I,60,167], [chr20,166350,-,1810S970M40I1429S,60,106]]
    new_sup = []
    for info in sup_info:
        info = info.split(',')
        if int(info[4]) >= min_mapq:
            cigar_list = list(cigar.Cigar(info[3]).items())  # [(2436, 'S'), (1758, 'M'), (55, 'I')]
            ref_start = int(info[1]) - 1
            ref_end = 0
            ori = 0 if info[2] == '-' else 1  # 方向0左1右
            for i in cigar_list:
                if i[1] == 'M':
                    ref_end = ref_start + int(i[0])
            new_sup.append([info[0], ref_start, ref_end, ori])
    return new_sup

def is_diff_ori(info):
    ori = info[0][2]
    for i in info[1:]:
        if i[2] != ori:
            return True
    return False

def discard_del(ref_start, ref_end, cigar, value, sv_size=1):
    normal_list = list()
    start = end = ref_start
    for c in cigar:
        if c[0] in [0, 7, 8]:
            end += c[1]
        elif c[0] == 2:
            if c[1] >= sv_size:
                if end - 1 >= start:
                    normal_list.append([start, end - 1, value])
                start = end + c[1]
                end = start
            else:
                end += c[1]
    if ref_end >= start:
        normal_list.append([start, ref_end, value])
    return normal_list


def encode_indel(bam, chrname, bppair, min_mapq, svtype):
    bp1, bp2 = bppair[0], bppair[1]
    search_range = bp2 - bp1
    s1, e1 = bp1 - search_range, bp1
    s2, e2 = bp2, bp2 + search_range
    read_base = 0
    ref_base = 3 * search_range
    color_list = list()  # [start, end, value]
    read_dic = dict()
    # 1) analysis cigar
    for read in bam.fetch(chrname, s1, e2):
        if read.mapq <= min_mapq:
            continue
        name = read.query_name
        ref_start = read.reference_start
        ref_end = read.reference_end
        query_start = read.query_alignment_start
        query_end = read.query_alignment_end
        if name not in read_dic:
            read_dic[name] = [[query_start, query_end, ref_start, ref_end]]
        else:
            read_dic[name].append([query_start, query_end, ref_start, ref_end])
        cigar = read.cigartuples
        # statistic indel
        shift_ins = 0
        shift_del = 0
        start = ref_start
        ins_list = list()
        del_list = list()
        normal_list = list()
        for c in cigar:
            if c[0] in [0, 2, 7, 8]:
                shift_ins += c[1]
            if c[0] == 1:
                ins_list.append([ref_start + shift_ins, ref_start + shift_ins + c[1] - 1])
            if c[0] in [0, 7, 8]:
                shift_del += c[1]
            if c[0] == 2:
                normal_list.append([start, ref_start + shift_del - 1])
                del_list.append([ref_start + shift_del, ref_start + shift_del + c[1] - 1])
                shift_del += c[1]
                start = ref_start + shift_del
        if ref_end - start >= 0:
            normal_list.append([start, ref_end])
        sv_list = list()

    # 2) analysis split read
        for k, v in read_dic.items():
            if len(v) <= 1:
                continue
            v = sorted(v, key=lambda x: x[0])       # query_start, query_end, ref_start, ref_end
            for i in range(len(v)-1):
                if (v[i+1][0]-v[i][1]) - (v[i+1][2]-v[i][3]) >= 50:
                    ins_len = (v[i+1][0]-v[i][1]) - (v[i+1][2]-v[i][3])
                    ins_list.append([v[i][3], v[i][3]+ins_len])
                elif (v[i+1][2]-v[i][3]) - (v[i+1][0]-v[i][1]) >= 50:
                    del_list.append([v[i][3], v[i+1][2]])
        if svtype.upper() in ['INS']:
            sv_list = ins_list
        elif svtype.upper() in ['DEL']:
            sv_list = del_list
        for sv in sv_list:
            s, e = sv
            value = 0
            if e-s <= e2 - s1:
                value += 1
                if bp1 <= (s+e)//2 <= bp2:
                    value += 2
                    if e-s <= 2*(bp2 - bp1):
                        value += 4
                        read_base += e - s + 1
            color_list.append([s, e, value])
        for sv in normal_list:
            s, e = sv
            read_base += max(0, min(e, e2)-max(s, s1))
    depth = read_base//ref_base
    return color_list, depth

def encode_inv(bam, chrname, bppair, min_mapq):
    bp1, bp2 = bppair[0], bppair[1]
    search_range = bp2 - bp1  # bp1     bp2
    s1, e1 = bp1 - search_range, bp1  # s1      e1      s2      e2
    s2, e2 = bp2, bp2 + search_range
    read_base = 0
    ref_base = 2 * search_range
    color_list = list()  # [start, end, value]
    ori_state = dict()
    for read in bam.fetch(chrname, max(bp1 - search_range, 0), bp2 + search_range):
        read_name = read.query_name
        ref_start = read.reference_start
        ref_end = read.reference_end
        ref_ori = read.is_reverse
        cigar = read.cigartuples
        # value = 0
        # if cigar[0][0] == 4 or cigar[0][0] == 5 or cigar[-1][0] == 4 or cigar[-1][0] == 5:
        #     value += 1
        value = 1
        state = -1
        if read_name not in ori_state:
            tags = read.get_tags()
            is_sa = False
            sup_info = list()
            for tag in tags:
                if tag[0] == 'SA':
                    is_sa = True
                    sup_info = tag[1].split(';')[:-1]
                    break
            if is_sa:
                sup_info = deal_sup(sup_info, min_mapq)
                info = [[chrname, ref_start, ref_end, ref_ori]] + sup_info
                info_new = list()
                for i in info:
                    if i[0] == read.reference_name:
                        if min(e2, i[2]) - max(s1, i[1]) > 0:
                            info_new.append(i)
                diff_ori = is_diff_ori(info_new)
                state = 1 if diff_ori else 0
                ori_state[read_name] = state
        else:
            state = ori_state[read_name]
        if bp1 <= (ref_start + ref_end) // 2 <= bp2:
            value += 2
            if state == 1:
                value += 4
        color_list += discard_del(ref_start, ref_end, cigar, value)
        d1 = min(ref_end, e1) - max(ref_start, s1)
        d2 = min(ref_end, e2) - max(ref_start, s2)
        if d1 > 0:
            read_base += d1
        if d2 > 0:
            read_base += d2
    depth = read_base // ref_base
    return color_list, depth

def encode_dup(bam, chrname, bppair, min_mapq):
    bp1, bp2 = bppair[0], bppair[1]
    search_range = bp2 - bp1
    s1, e1 = bp1 - search_range, bp1
    s2, e2 = bp2, bp2 + search_range
    read_base = 0
    ref_base = 3 * search_range
    color_list = list()  # [start, end, value]
    read_dic = dict()
    # 1) analysis cigar
    for read in bam.fetch(chrname, s1, e2):
        if read.mapq <= min_mapq:
            continue
        name = read.query_name
        ref_start = read.reference_start
        ref_end = read.reference_end
        query_start = read.query_alignment_start
        query_end = read.query_alignment_end
        if name not in read_dic:
            read_dic[name] = [[query_start, query_end, ref_start, ref_end]]
        else:
            read_dic[name].append([query_start, query_end, ref_start, ref_end])
        cigar = read.cigartuples
        # statistic indel
        shift_ins = 0
        shift_del = 0
        start = ref_start
        ins_list = list()
        del_list = list()
        normal_list = list()
        for c in cigar:
            if c[0] in [0, 2, 7, 8]:
                shift_ins += c[1]
            if c[0] == 1:
                ins_list.append([ref_start + shift_ins, ref_start + shift_ins + c[1] - 1])
            if c[0] in [0, 7, 8]:
                shift_del += c[1]
            if c[0] == 2:
                normal_list.append([start, ref_start + shift_del - 1])
                del_list.append([ref_start + shift_del, ref_start + shift_del + c[1] - 1])
                shift_del += c[1]
                start = ref_start + shift_del
        if ref_end - start >= 0:
            normal_list.append([start, ref_end])
        dup_list = list()

        # analysis dup
        for k, v in read_dic.items():
            if len(v) <= 1:
                continue
            v = sorted(v, key=lambda x: x[0])  # query_start, query_end, ref_start, ref_end
            for i in range(len(v) - 1):
                if min(v[i+1][3], v[i][3]) - max(v[i+1][2], v[i][2]) >= 0:
                    dup_list += [[v[i][2], v[i][3]], [v[i+1][2], v[i+1][3]]]
                # if (v[i+1][0]-v[i][1]) - (v[i+1][2]-v[i][3]) >= 50:
                    # ins_len = (v[i+1][0]-v[i][1]) - (v[i+1][2]-v[i][3])
                    # ins_list.append([v[i][2], v[i][2]+ins_len])
        for sv in dup_list+ins_list:
            s, e = sv
            value = 0
            if e - s <= e2 - s1:
                value += 1
                if bp1 <= (s + e) // 2 <= bp2:
                    value += 2
                    if e - s <= 2 * (bp2 - bp1):
                        value += 4
                        read_base += e - s + 1
            color_list.append([s, e, value])
        for sv in normal_list:
            s, e = sv
            read_base += max(0, min(e, e2) - max(s, s1))
    depth = read_base // ref_base
    return color_list, depth

def shift_index_indel(color_list, bppair):
    s1, e1 = bppair[0], bppair[1]           #s2     s1      e1      e2
    sv_len = e1 - s1 + 1
    s2 = s1 - sv_len
    e2 = e1 + sv_len
    shift = list()
    for read in color_list:
        start, end, color = read
        if min(end, e2) - max(start, s2) <= 0:
            continue
        shift_start = max(start, s2) - s2
        shift_end = min(end, e2) - s2
        shift.append([shift_start, shift_end, color])
    return shift, sv_len*3

def shift_index_inv(color_list, bppair):
    threshold = 1000
    bp1, bp2 = bppair[0:2]
    shift = list()
    if bp2 - bp1 <= threshold:
        threshold = bp2 - bp1
        s, e = bp1 - threshold, bp2 + threshold
        for read in color_list:
            start, end, color = read
            if min(end, e) - max(start, s) <= 0:
                continue
            shift_start = max(start, s) - s
            shift_end = min(end, e) - s
            shift.append([shift_start, shift_end, color])
    else:
        s1, e1 = bp1 - threshold, bp1 + threshold // 2
        s2, e2 = bp2 - threshold // 2, bp2 + threshold
        for read in color_list:
            start, end, color = read
            if min(end, e1) - max(start, s1) > 0:
                shift_start1 = max(start, s1) - s1
                shift_end1 = min(end, e1) - s1
                shift.append([shift_start1, shift_end1, color])
            if min(end, e2) - max(start, s2) > 0:
                shift_start2 = max(start, s2) - s2 + threshold//2*3
                shift_end2 = min(end, e2) - s2 + threshold//2*3
                shift.append([shift_start2, shift_end2, color])
    return shift, threshold*3


def deal_draw_list(draw_list, draw_range):
    dl1 = list()
    dl2 = list()
    bound = (draw_range - 1) // 2
    for draw in draw_list:
        left, right, value = draw
        if 0 <= right <= bound:
            dl1.append(draw)
        elif bound + 1 <= left <= draw_range - 1:
            dl2.append([left - bound - 1, right - bound - 1, value])
        else:
            dl1.append([left, bound, value])
            dl2.append([0, right - bound - 1, value])
    return dl1, dl2

def generate_img(elelist, col, depth):
    imageback = np.zeros((len(elelist), col))
    index = 0
    for read in elelist:
        start, end, value = read
        imageback[index][start: end] = value
        index += 1
    imageback = np.sort(imageback, axis=0)[::-1]
    imageback = imageback[0:depth, :].tolist()  # depth = 是计算出来的
    for i in range(len(imageback)):
        for j in range(len(imageback[0])):
            x = imageback[i][j]
            if x == 0:
                imageback[i][j] = [0, 0, 0]
            elif x == 1:
                imageback[i][j] = [255, 0, 0]
            elif x == 2:
                imageback[i][j] = [0, 255, 0]
            elif x == 3:
                imageback[i][j] = [255, 255, 0]
            elif x == 4:
                imageback[i][j] = [0, 0, 255]
            elif x == 5:
                imageback[i][j] = [255, 0, 255]
            elif x == 6:
                imageback[i][j] = [0, 255, 255]
            else:
                imageback[i][j] = [255, 255, 255]
    imageback = np.array(imageback, dtype=np.uint8)
    return imageback

def draw_image(draw_list, draw_range, depth, length, width):
    bound = (draw_range-1)//2
    dl1, dl2 = deal_draw_list(draw_list, draw_range)
    img1 = generate_img(dl1, bound+1, depth)
    img2 = generate_img(dl2, draw_range-bound-1, depth)
    image = np.append(img1, img2, axis=1)
    image = cv2.resize(image, (length, width), interpolation=cv2.INTER_CUBIC)
    return image
import sys
import time
import logging
import argparse

def parse_info(seq):
    info = {'SVTYPE': '', 'SVLEN': 0, 'END': 0}  # 先考虑SVLEN，再考虑END
    for i in seq.split(';'):  # SVTYPE=INV;END=139008588;SUPPORT=13;STD_SPAN=4.4;STD_POS=1.5
        if i.split('=')[0] in ['SVLEN', 'END']:
            try:
                info[i.split('=')[0]] = abs(int(i.split('=')[1]))
            except:
                pass
        if i.split('=')[0] == 'SVTYPE':
            info[i.split('=')[0]] = i.split('=')[1][0:3]
    return info


def evaluate(sv1, sv2, offset, bias, svtype):
    if svtype == 'iidd':
        start1, end1 = sv1[0], sv1[1]
        start2, end2 = sv2[0], sv2[1]
        overlap = (min(end1, end2) - max(start1, start2))/ max(end2-start2, end1-start1)
        if overlap >= bias:
            return True
        return False
    elif svtype == 'bnd':  # [pos1, chr2, pos2, form, i, caller]
        if sv1[1] == sv2[1] and abs(sv1[0] - sv2[0]) < offset and abs(sv1[2] - sv2[2]) < offset:
            return True
        else:
            return False


def parse_bnd(info):
    if info[0] == ']':
        form = 1
        chr = info.split(':')[0][1:]
        pos = int(info.split(':')[1][:-2])
    elif info[0] == '[':
        form = 2
        chr = info.split(':')[0][1:]
        pos = int(info.split(':')[1][:-2])
    else:
        if info[1] == ']':
            form = 3
            chr = info.split(':')[0][2:]
            pos = int(info.split(':')[1][:-1])
        else:
            form = 4
            chr = info.split(':')[0][2:]
            pos = int(info.split(':')[1][:-1])
    return [chr, pos, form]


def parseArgs(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('input', type=str)
    parser.add_argument('output', type=str)
    parser.add_argument('-b', '--bias', type=float, default=0.7)
    parser.add_argument('-o', '--offset', type=int, default=1000)
    parser.add_argument('--callers', type=str, nargs='+', default=['cuteSV', 'Sniffles2', 'svim', 'pbsv'])  # 多参数输入
    parser.add_argument('--support', '-s', type=int, default=2)
    # parser.add_argument('--callers', type=str, nargs='+')  # 多参数输入
    # 多参数并且扩增
    # parser.add_argument('-i','--input',action='append',nargs=2, metavar=('url','name'))
    args = parser.parse_args(argv)
    return args


def main_ctrl(args):
    input, output = args.input, args.output
    callers = args.callers
    bias, offset, support = args.bias, args.offset, args.support

    index_info = list()  # 合并时用到的
    index_list = list()  # 最后用到的
    sv_temp = dict()  # {'del': {'chr' :[[bp1, bp2, index, caller]]}
    f = open(input).readlines()

    # 如果输入的有工具排序，则按照优先级排序
    if callers:
        rank = dict()
        for i in range(len(callers)):
            rank[callers[i]] = i

    for i, line in enumerate(f):
        if '#' in line:
            index_list.append(i)
            continue
        seq = line.strip().split('\t')
        chrname = seq[0]
        pos = int(seq[1])
        caller = seq[2].split('.')[0]
        if caller not in callers:
            continue
        info = parse_info(seq[7])
        svtype, svlen, end = info['SVTYPE'], info['SVLEN'], info['END']
        sv_info = list()
        if svtype in ['INS', 'DEL', 'INV', 'DUP']:
            if svlen == 0:
                svlen = end - pos
            sv_info = [pos, pos + svlen, i, caller]
        elif svtype in ['BND']:  # [pos1, chr2, pos2, form, i, caller]
            sv_info = [pos] + parse_bnd(seq[4]) + [i, caller]
        # 按照类型存到字典里面
        if svtype not in sv_temp:
            sv_temp[svtype] = dict()
        if chrname not in sv_temp[svtype]:
            sv_temp[svtype][chrname] = [sv_info]
        else:
            sv_temp[svtype][chrname].append(sv_info)

    # Begining combine
    for st, chr_list in sv_temp.items():
        if st in ['INS', 'DEL', 'INV', 'DUP']:
            for chrname, sv_list in sv_temp[st].items():  # 'chr': [pos, end, i, caller], [pos, end, i, caller]
                if len(sv_list) == 1:
                    index_info.append([[sv_list[0][2], sv_list[0][3]]])
                else:
                    last = sv_list[0]
                    index_info.append([[last[2], last[3]]])
                    for i in range(1, len(sv_list)):
                        # start, end = sv_list[i][0], sv_list[i][1]
                        # last_start, last_end = last[0], last[1]
                        if evaluate(last, sv_list[i], offset, bias, 'iidd'):
                            index_info[-1].append([sv_list[i][2], sv_list[i][3]])
                        else:
                            index_info.append([[sv_list[i][2], sv_list[i][3]]])
                        last = sv_list[i]
        elif st in ['BND']:  # [pos1, chr2, pos2, form, i, caller]
            for chrname, sv_list in sv_temp[st].items():
                if len(sv_list) == 1:
                    index_info.append([[sv_list[0][4], sv_list[0][5]]])
                else:                           # sort  form, chr2, pos1, pos2
                    sv_list = sorted(sv_list, key=lambda x: [x[1], x[0], x[2]])
                    last = sv_list[0]
                    index_info.append([[last[4], last[5]]])
                    for i in range(1, len(sv_list)):
                        if evaluate(last, sv_list[i], offset, bias, 'bnd'):
                            index_info[-1].append([sv_list[i][4], sv_list[i][5]])
                        else:
                            index_info.append([[sv_list[i][4], sv_list[i][5]]])
                        last = sv_list[i]

    for index in index_info:
        if len(index) < support:
            continue
        if callers:
            index = sorted(index, key=lambda x: rank[x[1]])
        index_list.append(index[0][0])
    index_list.sort()

    f = open(input, 'r').readlines()
    with open(output, 'w') as f1:
        for i in index_list:
            f1.write(f[i])
    f1.close()


def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    logging.info("Running %s" % " ".join(sys.argv))


def run(argv):
    setupLogging()
    args = parseArgs(argv)
    starttime = time.time()
    main_ctrl(args)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))


if __name__ == '__main__':
    run(sys.argv[1:])

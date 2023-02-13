import sys
import time
import os
import logging
import argparse
import pysam
import cv2
from multiprocessing import Pool
from cnnLSV_combine import combine
from cnnLSV_encode import *
from cnnLSV_net import pre_sv


def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', type=str, help='BAM file.')
    parser.add_argument('input', type=str, help='Input SVs dataset called by callers.')
    parser.add_argument('output', type=str, help='Output SVs dataset filtered by cnnLSV.')
    parser.add_argument('--dataset', type=str, default='sim')
    parser.add_argument('-t', '--threads', type=int, default=4)
    parser.add_argument('-T', '--tempdir', type=str, help='Save tempory files.', default='cnnLSV_temp')
    parser.add_argument('--svtype', type=str, nargs='+', default=['INS', 'DEL', 'INV', 'DUP'])
    # combine
    parser.add_argument('--minsup', '-s', type=int, default=5)
    parser.add_argument('--bias', type=float, default=0.5)
    parser.add_argument('--callers', type=str, nargs='+', default=['cuteSV', 'svim', 'Sniffles2', 'pbsv'])
    # others
    parser.add_argument('--model', type=str, default='simmodel.pt')
    parser.add_argument('--width', type=int, default=100)
    parser.add_argument('--length', type=int, default=200)
    parser.add_argument('--minmapq', type=int, default=20)
    parser.add_argument('--batchsize', type=int, default=50)
    parser.add_argument('--minsvlen', type=int, default=50)
    # cnn structural
    parser.add_argument('--dropout', type=float, default=0.5)
    parser.add_argument('--line', type=int, nargs='+', default=[900, 128, 128, 1])
    # SVTYPE: INDEL
    parser.add_argument('--maxsvlen', type=int, default=50000)
    # parser.add_argument('-i', '--input', action='append', nargs='+', default=[['INS, DEL']])
    # metavar = ('url', 'name')
    args = parser.parse_args(argv)
    return args

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

def filter_sv(bam_fp, tempdir, task_sv, length, width, minmapq):
    bam = pysam.AlignmentFile(bam_fp)
    task, sv_list = task_sv
    true_list = list()
    for sv in sv_list:
        chrname, start, end, svtype, info, index = sv  # [chr, pos, end, svtype]
        bp = [start, end]
        svlen = end - start
        minlen, maxlen, pt, dropout, line, mean, std = info['minsvlen'], info['maxsvlen'], info['model'], info['dropout'], info['line'], info['mean'], info['std']
        if svlen > maxlen or svlen < minlen:
            true_list.append(index)
        else:
            try:
                draw_list = list()
                draw_range, depth = 0, 0
                if svtype == 'DEL' or svtype == 'INS':
                    color_list, depth = encode_indel(bam, chrname, bp, minmapq, svtype)
                    draw_list, draw_range = shift_index_indel(color_list, bp)  # del没有width
                elif svtype == 'DUP':
                    color_list, depth = encode_dup(bam, chrname, bp, minmapq)
                    draw_list, draw_range = shift_index_indel(color_list, bp)
                elif svtype == 'INV':
                    color_list, depth = encode_inv(bam, chrname, bp, minmapq)
                    draw_list, draw_range = shift_index_inv(color_list, bp)
                img = draw_image(draw_list, draw_range, depth, length, width)
                pre = pre_sv(img, pt, mean, std, line, dropout, svtype)
                if pre:
                    true_list.append(index)
            except:
                true_list.append(index)
    with open(tempdir + '/' + task + '.txt', 'w') as f:
        f.write(str(true_list))
    f.close()
    bam.close()


def multi_filter_sv(args):
    return filter_sv(*args)

def main_ctrl(args):
    # 1. Preparing
    logging.info('Parse the parameters.')

    mean_std = eval(open('mean_std.txt').read())
    dataset = args.dataset

    sv_info = dict()
    for svtype in ['INS', 'DEL', 'DUP', 'INV']:
        sv_info[svtype] = {'minsvlen' : args.minsvlen,
                           'maxsvlen' : args.maxsvlen,
                           'model' : args.model,
                           'dropout' : args.dropout,
                           'line' : args.line,
                           'mean' : mean_std[dataset][svtype]['mean'],
                           'std' : mean_std[dataset][svtype]['std']}

    st_input = [i.upper() for i in args.svtype]         # Users' input
    st_filter = list(set(st_input) & set(sv_info.keys()))  # svtype should filter
    # for t in set(st_input) - set(sv_info.keys()):
    #     logging.info('cnnLSV does not support filter the {}!'.format(t.upper()))
    logging.info('cnnLSV will filter {} from the input file.'.format('&'.join(st_filter)))

    # 2. Analysis input_file
    # combine vcf
    logging.info('Analysis input file.')
    line_index = list()
    filter_index = list()
    comb_index = combine(args.input, args.callers, args.bias)
    for i in comb_index:
        if len(i) < args.minsup:
            filter_index += i
        else:
            line_index += i
    logging.info('{} SVs are supported by less than {} callers.'.format(len(filter_index), args.minsup))
    batchsize = args.batchsize
    task_dict = dict()  # [[[chr, pos, end, svtype],....]]
    count = 0
    for i, line in enumerate(open(args.input, 'r').readlines()):
        if '#' in line:
            line_index.append(i)
            continue
        elif i in filter_index:
            seq = line.split('\t')
            chrname, pos = seq[0], int(seq[1])
            info = parse_info(seq[7])
            svtype, svlen, end = info['SVTYPE'], info['SVLEN'], info['END']
            if svtype not in st_filter:
                line_index.append(i)
                continue
            if svlen == 0:
                svlen = end - pos
            if svtype in st_filter and svlen >= args.minsvlen:
                if count % batchsize == 0:
                    task_dict['task' + str(count // batchsize)] = []
                task_dict['task' + str(count // batchsize)].append([chrname, pos, pos + svlen, svtype, sv_info[svtype], i])
                count += 1
    logging.info('There are {} SVs waiting for filtering.'.format(count))

    # 3. Filter SV
    if not os.path.exists(args.tempdir):
        os.mkdir(args.tempdir)
    logging.info('Begining filter SVs.')
    pools = Pool(processes=int(args.threads))
    for key, value in task_dict.items():
        para = [(args.bam, args.tempdir, [key, value], args.length, args.width, args.minmapq)]
        pools.map_async(multi_filter_sv, para)
    pools.close()
    pools.join()

    # 4. Write VCF
    logging.info('Combine and write output file.')
    for f in os.listdir(args.tempdir):
        line_index += eval(open(args.tempdir + '/' + f).read())
        os.remove(args.tempdir + '/' + f)
    os.rmdir(args.tempdir)
    line_index.sort()
    input = open(args.input).readlines()
    with open(args.output, 'w') as f:
        for i in line_index:
            f.write(input[i])
    f.close()


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

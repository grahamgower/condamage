#!/usr/bin/env python
# Plot post-mortem deamination patterns, from `condamage' output.
#
# Copyright (c) 2018 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

from __future__ import print_function
import sys
import os.path
import collections

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

def parse_condamage(fn):
    data = collections.defaultdict(lambda: ([],[]))
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == "#":
                continue
            fields = line.split()
            ctx = fields[0]
            i, mm, n = map(int, fields[1:])

            if n == 0:
                continue
            data[ctx][0].append(i)
            data[ctx][1].append(float(mm)/n)
    return data

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot histogram for 1 column of input file")
    parser.add_argument("-u", action="store_false", default=True, help="plot unconditional mismatches [%(default)s]")
    parser.add_argument("-5", dest="c5", action="store_false", default=True, help="plot mismatches conditional on 5' C>T [%(default)s]")
    parser.add_argument("-3", dest="c3", action="store_false", default=True, help="plot mismatches conditional on 3' G>A [%(default)s]")
    parser.add_argument("-s", "--scale", type=float, default=1.5, help="scale the size of the plot [%(default)s]")
    parser.add_argument("-w", "--wide", action="store_false", default=True, help="plot widescreen ratio (16x9) [%(default)s]")
    parser.add_argument("-o", "--opdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("--title", type=str, help="plot title")
    parser.add_argument("infile", help="input file")
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = parse_args()
    data = parse_condamage(args.infile)

    pdf = PdfPages(args.opdf)
    if args.wide:
        fig_w, fig_h = plt.figaspect(9.0/16.0)
    else:
        fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(args.scale*fig_w, args.scale*fig_h))
    gs1 = gridspec.GridSpec(1, 2)
    ax1 = fig1.add_subplot(gs1[0])
    ax2 = fig1.add_subplot(gs1[1], sharey=ax1)
    ax2.invert_xaxis()

    if not args.title:
        args.title = os.path.basename(args.infile)

    ctxlist = []
    lslist = []
    mlist = []
    labels = []
    if args.u:
        ctxlist.extend(["C2T", "G2A"])
        lslist.extend(["-", "-"])
        mlist.extend(["o", "s"])
        labels.extend(["$C \\rightarrow T$",
                        "$G \\rightarrow A$"])
    if args.c5:
        ctxlist.extend(["COND5_C2T", "COND5_G2A"])
        lslist.extend([":", ":"])
        mlist.extend(["<", ">"])
        labels.extend(["$C \\rightarrow T\,\ \mid\ 5^{\prime}\ C \\rightarrow T$",
                        "$G \\rightarrow A\ \mid\ 5^{\prime}\ C \\rightarrow T$"])
    if args.c3:
        ctxlist.extend(["COND3_C2T", "COND3_G2A"])
        lslist.extend(["--", "--"])
        mlist.extend(["d", "*"])
        labels.extend(["$C \\rightarrow T\,\ \mid\ 3^{\prime}\ G \\rightarrow A$",
                        "$G \\rightarrow A\ \mid\ 3^{\prime}\ G \\rightarrow A$"])

    xmax = 0
    ymax = 0
    for ctx, ls, marker, label in zip(ctxlist,lslist,mlist,labels):
        ctx5 = ctx+"5"
        ctx3 = ctx+"3"
        assert ctx5 in data, "{} missing from {}".format(ctx5, args.infile)
        assert ctx3 in data, "{} missing from {}".format(ctx3, args.infile)
        ax1.plot(data[ctx5][0], data[ctx5][1], marker=marker, linestyle=ls, label=label)
        ax2.plot(data[ctx3][0], data[ctx3][1], marker=marker, linestyle=ls, label=label)

        if data[ctx5][0][-1] > xmax:
            xmax = data[ctx5][0][-1]
        if data[ctx3][0][-1] > xmax:
            xmax = data[ctx3][0][-1]

        if ctx != "COND5_C2T":
            if max(data[ctx5][1]) > ymax:
                ymax = max(data[ctx5][1])
        if ctx != "COND3_G2A":
            if max(data[ctx3][1]) > ymax:
                ymax = max(data[ctx3][1])

    ax2.set_ylim([0, ymax*1.1])

    xticks = list(range(1,xmax+1))
    xticklabels = [str(x) for x in xticks]
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticklabels, rotation='vertical')
    ax2.set_xticks(xticks)
    #ax2.set_xticklabels(xticklabels, rotation='vertical')
    xticklabels2 = [str(-1*x) for x in xticks]
    ax2.set_xticklabels(xticklabels2, rotation='vertical')

    fig1.suptitle(args.title)
    ax1.set_xlabel("Mismatches towards 5' end", labelpad=10)
    ax2.set_xlabel("Mismatches towards 3' end", labelpad=10)
    ax1.set_ylabel("Frequency")
    ax2.yaxis.tick_right()
    ax1.legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf.savefig()
    pdf.close()

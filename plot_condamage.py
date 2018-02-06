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
    parser = argparse.ArgumentParser(description="Plot conditional damage profile, from `condamage' output.")
    parser.add_argument("-a", "--all", action="store_true", default=False, help="plot all conditional traces [%(default)s]")
    parser.add_argument("-s", "--singlestranded", action="store_true", default=False, help="plot mismatches conditional on C>T at either end [%(default)s]")
    parser.add_argument("--scale", type=float, default=1.5, help="scale the size of the plot [%(default)s]")
    parser.add_argument("--wide", action="store_false", default=True, help="plot widescreen ratio (16x9) [%(default)s]")
    parser.add_argument("-o", "--opdf", metavar="out.pdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("-t", "--title", type=str, help="title for the plot")
    parser.add_argument("infile", metavar="condamage.txt", help="input file, produced by `condamage'")
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

    col = iter(plt.get_cmap("tab20").colors)
    alpha = 1.0

    plotmeta = collections.OrderedDict([
            # key: (edgecolour, facecolour, linestyle, marker, label)
            ("C2T?", (next(col), next(col), "-", "o", "$C \\rightarrow T$")),
            ("G2A?", (next(col), next(col), "--", "s", "$G \\rightarrow A$")),

            ("C2T?|5C2T", (next(col), next(col), ":", "<", "$C \\rightarrow T\,\ \mid\ 5^{\prime}\ C \\rightarrow T$")),
            ("G2A?|5C2T", (next(col), next(col), ":", ">", "$G \\rightarrow A\ \mid\ 5^{\prime}\ C \\rightarrow T$")),
            ("C2T?|3C2T", (next(col), next(col), ":", "x", "$C \\rightarrow T\,\ \mid\ 3^{\prime}\ C \\rightarrow T$")),
            ("G2A?|3C2T", (next(col), next(col), ":", "d", "$G \\rightarrow A\ \mid\ 3^{\prime}\ C \\rightarrow T$")),

            ("C2T?|5G2A", (next(col), next(col), ":", "^", "$C \\rightarrow T\,\ \mid\ 5^{\prime}\ G \\rightarrow A$")),
            ("G2A?|5G2A", (next(col), next(col), ":", "v", "$G \\rightarrow A\ \mid\ 5^{\prime}\ G \\rightarrow A$")),
            ("C2T?|3G2A", (next(col), next(col), ":", "*", "$C \\rightarrow T\,\ \mid\ 3^{\prime}\ G \\rightarrow A$")),
            ("G2A?|3G2A", (next(col), next(col), ":", "p", "$G \\rightarrow A\ \mid\ 3^{\prime}\ G \\rightarrow A$")),
        ])

    if args.all:
        k5list = plotmeta.keys()
        k3list = plotmeta.keys()
    elif args.singlestranded:
        k5list = ["C2T?", "G2A?", "C2T?|3C2T"]
        k3list = ["C2T?", "G2A?", "C2T?|5C2T"]
    else:
        k5list = ["C2T?", "G2A?", "C2T?|3G2A"]
        k3list = ["C2T?", "G2A?", "G2A?|5C2T"]

    xmax = 0
    ymax = 0
    for k in k5list:
        k5 = k.replace("?", "5", 1)
        assert k5 in data, "{} missing from {}".format(k5, args.infile)
        ecol, col, ls, m, lbl = plotmeta[k]
        ax1.plot(data[k5][0], data[k5][1], color=col, markeredgecolor=ecol, linestyle=ls, marker=m, label=lbl, alpha=alpha)

        if data[k5][0][-1] > xmax:
            xmax = data[k5][0][-1]

        kmax = max(data[k5][1])
        if kmax != 1.0 and kmax > ymax:
            ymax = kmax

    for k in k3list:
        k3 = k.replace("?", "3", 1)
        assert k3 in data, "{} missing from {}".format(k3, args.infile)
        ecol, col, ls, m, lbl = plotmeta[k]
        ax2.plot(data[k3][0], data[k3][1], color=col, markeredgecolor=ecol, linestyle=ls, marker=m, label=lbl, alpha=alpha)

        if data[k3][0][-1] > xmax:
            xmax = data[k3][0][-1]

        kmax = max(data[k3][1])
        if kmax != 1.0 and kmax > ymax:
            ymax = kmax

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
    ax1.legend(loc="upper right")
    ax2.legend(loc="upper left")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf.savefig()
    pdf.close()

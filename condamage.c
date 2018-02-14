/*
 * Score post-mortem damage patterns.  Also score patterns conditional on
 * deamination at the most 5' position, and patterns conditional on
 * deamination at the most 3' position.
 *
 * Copyright (c) 2018 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#define CONDAMAGE_VERSION "2"

typedef struct {
	char *bam_fn; // input filename
	char *bam_ofn; // output filename
	char *fasta_fn;

	int argc;
	char **argv;

	int c5, c3, g5, g3; // -C and -G params
	int lmax; // maximum length for the fraglen histograms

	// Size of the region we'll record mismatches for,
	// separately for the start and the end of reads.
	size_t window;

	int fwd_only;
	int rev_only;
} opt_t;

/*
 * Append a line to the bam header.
 *
 * From the SAMv1 spec:
 * hdr->l_text   Length of the header text, including any NUL padding.
 * hdr->text     Plain header text in SAM; not necessarily NUL-terminated.
 */
static int
bam_hdr_append(bam_hdr_t *hdr, const char *line)
{
	int oldlen;
	int newlen;
	int addnl;
	char *c, *tmp;

	for (oldlen=0, c=hdr->text;
		oldlen<hdr->l_text && *c!='\0';
		oldlen++, c++)
		;

	if (hdr->text[oldlen-1] != '\n')
		addnl = 1;
	else
		addnl = 0;

	newlen = oldlen + addnl + strlen(line) + 2;
	tmp = realloc(hdr->text, newlen);
	if (tmp == NULL) {
		perror("realloc:bam_hdr_append");
		return -1;
	}
	hdr->text = tmp;
	hdr->l_text = newlen;

	if (addnl) {
		sprintf(hdr->text+oldlen, "\n");
		oldlen++;
	}
	sprintf(hdr->text+oldlen, "%s\n", line);

	return 0;
}

/*
 * Load reference sequence, if required.
 */
static int
get_refseq(char **ref, faidx_t *fai, bam_hdr_t *bam_hdr, int tid)
{
	static int ref_tid = -1;
	static int ref_len = -1;

	if (ref_tid != tid) {
		if (*ref)
			free(*ref);
		*ref = NULL;
		ref_len = faidx_seq_len(fai, bam_hdr->target_name[tid]);
		if (ref_len == -1) {
			fprintf(stderr, "bam has region `%s', which is not in fasta file\n",
					bam_hdr->target_name[tid]);
			return -1;
		}
		*ref = faidx_fetch_seq(fai, bam_hdr->target_name[tid], 0, ref_len, &ref_len);
		if (*ref == NULL)
			return -1;
		ref_tid = tid;
	}

	return ref_len;
}

int
condamage(opt_t *opt)
{
	int i, j, k;
	int ret;

	samFile *bam_fp, *bam_ofp = NULL;
	bam_hdr_t *bam_hdr, *bam_ohdr = NULL;
	faidx_t *fai;
	char *ref = NULL;
	bam1_t *b;

	uint64_t *lhist, *lhist_cond; // fragment length counts

	enum {_5C2T=0, _3C2T, _5G2A, _3G2A};
#define COND_5C2T (1<<_5C2T)
#define COND_3C2T (1<<_3C2T)
#define COND_5G2A (1<<_5G2A)
#define COND_3G2A (1<<_3G2A)

	struct counts {
		// unconditional
		uint64_t c, c2t, g, g2a;

		// conditional
		struct {
			uint64_t c, c2t, g, g2a;
		} cond[4];

	} *counts5, // counts for the window towards the 5' end
	  *counts3; // counts for the window towards the 3' end

	counts5 = calloc(opt->window, sizeof(*counts5));
	if (counts5 == NULL) {
		perror("calloc:counts5");
		ret = -1;
		goto err0;
	}

	counts3 = calloc(opt->window, sizeof(*counts3));
	if (counts3 == NULL) {
		perror("calloc:counts3");
		ret = -2;
		goto err1;
	}

	lhist = calloc(opt->lmax, sizeof(*lhist));
	if (lhist == NULL) {
		perror("calloc:lhist");
		ret = -3;
		goto err2;
	}

	lhist_cond = calloc(opt->lmax, 4*sizeof(*lhist));
	if (lhist_cond == NULL) {
		perror("calloc:lhist_cond");
		ret = -4;
		goto err3;
	}

	bam_fp = sam_open(opt->bam_fn, "r");
	if (bam_fp == NULL) {
		fprintf(stderr, "bam_open: %s: %s\n", opt->bam_fn, strerror(errno));
		ret = -5;
		goto err4;
	}

	bam_hdr = sam_hdr_read(bam_fp);
	if (bam_hdr == NULL) {
		fprintf(stderr, "%s: couldn't read header\n", opt->bam_fn);
		ret = -6;
		goto err5;
	}

	fai = fai_load(opt->fasta_fn);
	if (fai == NULL) {
		ret = -7;
		goto err6;
	}

	b = bam_init1();
	if (b == NULL) {
		ret = -8;
		goto err7;
	}

	if (opt->bam_ofn) {
		bam_ofp = sam_open(opt->bam_ofn, "w");
		if (bam_ofp == NULL) {
			fprintf(stderr, "bam_open: %s: %s\n", opt->bam_ofn, strerror(errno));
			ret = -9;
			goto err8;
		}

		bam_ohdr = bam_hdr_dup(bam_hdr);
		if (bam_ohdr == NULL) {
			/*
			 * XXX: bam_hdr_dup doesn't properly check for
			 * allocation failures, so we'll crash before
			 * getting here.
			 */
			fprintf(stderr, "bam_hdr_dup: failed to allocate memory: %s\n", strerror(errno));
			ret = -10;
			goto err9;
		}

#define BUFLEN 8096
		char pgbuf[BUFLEN];
		int offs = snprintf(pgbuf, BUFLEN, "@PG\tID:condamage\tPN:condamage\tVN:%s\tCL:", CONDAMAGE_VERSION);
		for (i=0; i<opt->argc && offs < BUFLEN-1; i++)
			offs += snprintf(pgbuf+offs, BUFLEN-offs, "%s%s", i==0?"":" ", opt->argv[i]);

		if (bam_hdr_append(bam_ohdr, pgbuf) < 0) {
			ret = -11;
			goto err10;
		}

		if (sam_hdr_write(bam_ofp, bam_ohdr) < 0) {
			fprintf(stderr, "sam_hdr_write: %s: %s\n", opt->bam_ofn, strerror(errno));
			ret = -12;
			goto err10;
		}
	}

	while (1) {
		int r = sam_read1(bam_fp, bam_hdr, b);
		if (r < 0) {
			if (r == -1)
				break;
			fprintf(stderr, "sam_read1: %s: read failed\n", opt->bam_ofn);
			ret = -13;
			goto err10;
		}
		bam1_core_t *c = &b->core;

		if (c->flag & (BAM_FUNMAP|BAM_FQCFAIL|BAM_FDUP|
				BAM_FSECONDARY|BAM_FSUPPLEMENTARY))
			// skip these
			continue;

		if (c->flag & BAM_FPAIRED)
			continue;

		if (opt->fwd_only && bam_is_rev(b))
			continue;

		if (opt->rev_only && !bam_is_rev(b))
			continue;

		int op;
		int x, // offset in ref
		    y; // offset in query seq
		char c1, c2;

		int hclip = 0;
		int b_out = 0;
		int cond = 0;

		uint8_t *seq = bam_get_seq(b);
		uint32_t *cigar = bam_get_cigar(b);
		int ref_len = get_refseq(&ref, fai, bam_hdr, c->tid);

		if (ref_len == -1) {
			ret = -14;
			goto err10;
		}

		if (bam_endpos(b) > ref_len) {
			fprintf(stderr, "%s: read mapped outside the reference sequence: bam/ref mismatch?\n",
					bam_get_qname(b));
			continue;
		}

		// check for mismatch at left most position
		op = bam_cigar_op(cigar[0]);
		if (op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF) {
			c1 = seq_nt16_str[bam_seqi(seq, 0)];
			c2 = ref[c->pos];

			if (c2 == 'C' && c1 == 'T') {
				if (bam_is_rev(b))
					cond |= COND_3G2A;
				else
					cond |= COND_5C2T;
			}

			if (c2 == 'G' && c1 == 'A') {
				if (bam_is_rev(b))
					cond |= COND_3C2T;
				else
					cond |= COND_5G2A;
			}
		}

		// check for mismatch at right most position
		op = bam_cigar_op(cigar[c->n_cigar-1]);
		if (op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF) {
			c1 = seq_nt16_str[bam_seqi(seq, c->l_qseq-1)];
			c2 = ref[bam_endpos(b)-1];

			if (c2 == 'G' && c1 == 'A') {
				if (bam_is_rev(b))
					cond |= COND_5C2T;
				else
					cond |= COND_3G2A;
			}

			if (c2 == 'C' && c1 == 'T') {
				if (bam_is_rev(b))
					cond |= COND_5G2A;
				else
					cond |= COND_3C2T;
			}
		}

		for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
			int op = bam_cigar_op(cigar[i]);
			int l = bam_cigar_oplen(cigar[i]);

			if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
				for (j = 0; j < l; ++j) {
					int z1 = y + j;
					int z2 = c->l_qseq - (z1 + 1);


					if (z1>=opt->window && z2>=opt->window)
						continue;

					c1 = seq_nt16_str[bam_seqi(seq, z1)];
					c2 = ref[x+j];

// update counts
#define c_update(cnt, z, var) \
	do { \
		unsigned k; \
		cnt[z].var++; \
		for (k=0; k<4; k++) { \
			if (cond & (1<<k)) \
				cnt[z].cond[k].var++; \
		} \
	} while (0)

					// ref has C
					if (c2 == 'C') {

						if (bam_is_rev(b)) {
							if (z1 < opt->window) {
								c_update(counts3, z1, g);
							} else if (z2 < opt->window) {
								c_update(counts5, z2, g);
							}
						} else {
							if (z1 < opt->window) {
								c_update(counts5, z1, c);
							} else if (z2 < opt->window) {
								c_update(counts3, z2, c);
							}
						}

						// read has T: C->T for fwd reads, G->A for rev reads
						if (c1 == 'T') {
							if (bam_is_rev(b)) {
								if (z1 < opt->window) {
									c_update(counts3, z1, g2a);
								} else if (z2 < opt->window) {
									c_update(counts5, z2, g2a);
								}
								if (z1 < opt->g3)
									b_out = 1;
								if (z2 < opt->g5)
									b_out = 1;
							} else {
								if (z1 < opt->window) {
									c_update(counts5, z1, c2t);
								} else if (z2 < opt->window) {
									c_update(counts3, z2, c2t);
								}
								if (z1 < opt->c5)
									b_out = 1;
								if (z2 < opt->c3)
									b_out = 1;
							}
						}
					}

					// ref has G
					if (c2 == 'G') {

						if (bam_is_rev(b)) {
							if (z1 < opt->window) {
								c_update(counts3, z1, c);
							} else if (z2 < opt->window) {
								c_update(counts5, z2, c);
							}
						} else {
							if (z1 < opt->window) {
								c_update(counts5, z1, g);
							} else if (z2 < opt->window) {
								c_update(counts3, z2, g);
							}
						}

						// read has A: G->A for fwd reads, C->T for rev reads
						if (c1 == 'A') {
							if (bam_is_rev(b)) {
								if (z1 < opt->window) {
									c_update(counts3, z1, c2t);
								} else if (z2 < opt->window) {
									c_update(counts5, z2, c2t);
								}
								if (z1 < opt->c3)
									b_out = 1;
								if (z2 < opt->c5)
									b_out = 1;
							} else {
								if (z1 < opt->window) {
									c_update(counts5, z1, g2a);
								} else if (z2 < opt->window) {
									c_update(counts3, z2, g2a);
								}
								if (z1 < opt->g5)
									b_out = 1;
								if (z2 < opt->g3)
									b_out = 1;
							}
						}
					}

				}
				x += l;
				y += l;

			} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {
				y += l;
			} else if (op == BAM_CREF_SKIP || op == BAM_CDEL) {
				x += l;
			} else if (op == BAM_CHARD_CLIP) {
				hclip += l;
			}

		}

		if (bam_ofp && b_out) {
			if (sam_write1(bam_ofp, bam_ohdr, b) < 0) {
				fprintf(stderr, "sam_write1: %s: write failed\n", opt->bam_ofn);
				ret = -15;
				goto err10;
			}
		}

		y += hclip;

		if (y < opt->lmax) {
			lhist[y]++;
			for (k=0; k<4; k++) {
				if (cond & (1<<k))
					lhist_cond[y<<2 | k]++;
			}
		}

	}

	printf("#condamage version %s\n", CONDAMAGE_VERSION);
	printf("#cmdline:");
	for (i=0; i<opt->argc; i++)
		printf(" %s", opt->argv[i]);
	printf("\n\n");

	// unconditional stats
	printf("#C2T5\ti\tmm\tn\n");
	printf("# C2T5  C to T mismatches towards the 5' end\n");
	printf("# i     distance from 5' end\n");
	printf("# mm    number of mismatches\n");
	printf("# n     matches+mismatches (ref has C)\n\n");
	for (i=0; i<opt->window; i++)
		printf("C2T5\t%d\t%jd\t%jd\n", i+1, (uintmax_t)counts5[i].c2t, (uintmax_t)counts5[i].c);
	printf("\n");

	printf("#C2T3\ti\tmm\tn\n");
	printf("# C2T3  C to T mismatches towards the 3' end\n");
	printf("# i     distance from 3' end\n");
	printf("# mm    number of mismatches\n");
	printf("# n     matches+mismatches (ref has C)\n\n");
	for (i=0; i<opt->window; i++)
		printf("C2T3\t%d\t%jd\t%jd\n", i+1, (uintmax_t)counts3[i].c2t, (uintmax_t)counts3[i].c);
	printf("\n");

	printf("#G2A5\ti\tmm\tn\n");
	printf("# G2A5  G to A mismatches towards the 5' end\n");
	printf("# i     distance from 5' end\n");
	printf("# mm    number of mismatches\n");
	printf("# n     matches+mismatches (ref has G)\n\n");
	for (i=0; i<opt->window; i++)
		printf("G2A5\t%d\t%jd\t%jd\n", i+1, (uintmax_t)counts5[i].g2a, (uintmax_t)counts5[i].g);
	printf("\n");

	printf("#G2A3\ti\tmm\tn\n");
	printf("# G2A3  G to A mismatches towards the 3' end\n");
	printf("# i     distance from 3' end\n");
	printf("# mm    number of mismatches\n");
	printf("# n     matches+mismatches (ref has G)\n\n");
	for (i=0; i<opt->window; i++)
		printf("G2A3\t%d\t%jd\t%jd\n", i+1, (uintmax_t)counts3[i].g2a, (uintmax_t)counts3[i].g);
	printf("\n");


	// conditional stats
	int win;
	for (win=0; win<2; win++) {
		char ch_win = "53"[win];
		struct counts *cnts;

		if (win == 0)
			cnts = counts5;
		else
			cnts = counts3;

		for (k=0; k<4; k++) {
			char *str_cond = ((char *[]){"5C2T", "3C2T", "5G2A", "3G2A"})[k];

			printf("#C2T%c|%s\ti\tmm\tn\n", ch_win, str_cond);
			printf("# C2T%c|%s  C to T mismatches towards the %c' end,\n", ch_win, str_cond, ch_win);
			printf("#            conditional on a %c to %c mismatch at the most %c' position\n", str_cond[1], str_cond[3], str_cond[0]);
			for (i=0; i<opt->window; i++)
				printf("C2T%c|%s\t%d\t%jd\t%jd\n", ch_win, str_cond, i+1, (uintmax_t)cnts[i].cond[k].c2t, (uintmax_t)cnts[i].cond[k].c);
			printf("\n");

			printf("#G2A%c|%s\ti\tmm\tn\n", ch_win, str_cond);
			printf("# G2A%c|%s  G to A mismatches towards the %c' end,\n", ch_win, str_cond, ch_win);
			printf("#            conditional on a %c to %c mismatch at the most %c' position\n", str_cond[1], str_cond[3], str_cond[0]);
			for (i=0; i<opt->window; i++)
				printf("G2A%c|%s\t%d\t%jd\t%jd\n", ch_win, str_cond, i+1, (uintmax_t)cnts[i].cond[k].g2a, (uintmax_t)cnts[i].cond[k].g);
			printf("\n");
		}
	}

	// fragment length histograms
	int lmax;
	for (lmax=opt->lmax; lmax>0 && lhist[lmax-1]==0; lmax--)
		;
	if (lmax > 0) {
		printf("#FL\tj\tk\tx1\tx2\tx3\tx4\n");
		printf("# FL  count of fragments with a given length\n");
		printf("# j   fragment length\n");
		printf("# k   number of fragments of length j\n");
		printf("# x1   number of fragments of length j with a 5' C->T \n");
		printf("# x2   number of fragments of length j with a 3' G->A \n");
		printf("# x3   number of fragments of length j with a 5' C->T \n");
		printf("# x4   number of fragments of length j with a 3' G->A \n");
		for (i=1; i<lmax; i++)
			printf("FL\t%d\t%zd\t%zd\t%zd\t%zd\t%zd\n", i, (uintmax_t)lhist[i],
					(uintmax_t)lhist_cond[i<<2 | _5C2T],
					(uintmax_t)lhist_cond[i<<2 | _3C2T],
					(uintmax_t)lhist_cond[i<<2 | _5G2A],
					(uintmax_t)lhist_cond[i<<2 | _3G2A]);
	}


	ret = 0;
err10:
	if (bam_ohdr)
		bam_hdr_destroy(bam_ohdr);
err9:
	if (bam_ofp)
		sam_close(bam_ofp);
err8:
	bam_destroy1(b);
//
	if (ref)
		free(ref);
err7:
	if (opt->fasta_fn)
		fai_destroy(fai);
err6:
	bam_hdr_destroy(bam_hdr);
err5:
	sam_close(bam_fp);
err4:
	free(lhist_cond);
err3:
	free(lhist);
err2:
	free(counts3);
err1:
	free(counts5);
err0:
	return ret;
}

static void
usage(const opt_t *opt)
{
	fprintf(stderr, "condamage v%s\n", CONDAMAGE_VERSION);
	fprintf(stderr, "usage: %s [...] in.bam ref.fasta\n", opt->argv[0]);
	fprintf(stderr, "\n");
	fprintf(stderr, "  -w INT       Size of the region for which (mis)matches are recorded [%zd]\n", opt->window);
	fprintf(stderr, "  -o FILE      BAM output filename [%s]\n", opt->bam_ofn?opt->bam_ofn:"");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -C INT,INT   Output reads with a C->T mismatch within INT bases of terminal\n");
	fprintf(stderr, "                position (5',3') [%d,%d]\n", opt->c5, opt->c3);
	fprintf(stderr, "                E.g. -C 3,3 outputs reads with C->T mismatches in any of the 3 \n");
	fprintf(stderr, "                positions at the start or end of a read.\n");
	fprintf(stderr, "                -C 1,1 is appropriate for single stranded libaries.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -G INT,INT   Output reads with a G->A mismatch within INT bases of terminal\n");
	fprintf(stderr, "                position (5',3') [%d,%d]\n", opt->g5, opt->g3);
	fprintf(stderr, "                E.g. -C 3,0 -G 0,3 outputs reads with C->T mismatches in any of\n");
	fprintf(stderr, "                the 3 positions at the start of the read, or a G->A mismatch in\n");
	fprintf(stderr, "                any of the 3 positions at the end of the read.\n");
	fprintf(stderr, "                -C 1,0 -G 0,1 is appropriate for double stranded libaries.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -l INT       Maximum length for fragment length histograms [%d]\n", opt->lmax);

	//fprintf(stderr, "  -f           Only consider reads aligned to the forward (ref) strand\n");
	//fprintf(stderr, "  -r           Only consider reads aligned to the reverse (non ref) strand\n");
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;

	memset(&opt, 0, sizeof(opt_t));
	opt.window = 30;
	opt.lmax = 1024;
	opt.argc = argc;
	opt.argv = argv;

	while ((c = getopt(argc, argv, "w:o:C:G:fr")) != -1) {
		switch (c) {
			case 'w':
				{
					unsigned long w = strtoul(optarg, NULL, 0);
					if (w < 0 || w > 100) {
						fprintf(stderr, "-w `%s' is invalid\n", optarg);
						usage(&opt);
					}
					opt.window = w;
				}
				break;
			case 'o':
				opt.bam_ofn = optarg;
				break;
			case 'l':
				{
					unsigned long l = strtoul(optarg, NULL, 0);
					if (l < 100 || l > 1024*1024) {
						fprintf(stderr, "-l `%s' is invalid\n", optarg);
						usage(&opt);
					}
					opt.lmax = l;
				}
				break;
			case 'f':
				opt.fwd_only = 1;
				break;
			case 'r':
				opt.rev_only = 1;
				break;
			case 'C':
				{
					char *tmp;
					unsigned long x1, x2;
					x1 = strtoul(optarg, &tmp, 0);
					if (x1 < 0 || x1 > 100 || tmp[0] != ',') {
						fprintf(stderr, "-C `%s' is invalid\n", optarg);
						usage(&opt);
					}

					x2 = strtoul(tmp+1, NULL, 0);
					if (x2 < 0 || x2 > 100) {
						fprintf(stderr, "-C `%s' is invalid\n", optarg);
						usage(&opt);
					}

					opt.c5 = x1;
					opt.c3 = x2;
				}
				break;
			case 'G':
				{
					char *tmp;
					unsigned long x1, x2;
					x1 = strtoul(optarg, &tmp, 0);
					if (x1 < 0 || x1 > 100 || tmp[0] != ',') {
						fprintf(stderr, "-G `%s' is invalid\n", optarg);
						usage(&opt);
					}

					x2 = strtoul(tmp+1, NULL, 0);
					if (x2 < 0 || x2 > 100) {
						fprintf(stderr, "-G `%s' is invalid\n", optarg);
						usage(&opt);
					}

					opt.g5 = x1;
					opt.g3 = x2;
				}
				break;
			default:
				usage(&opt);
		}
	}

	int cgsum = opt.c5+opt.c3+opt.g5+opt.g3;
	if (cgsum && opt.bam_ofn == NULL) {
		fprintf(stderr, "-C/-G specified, but no -o FILE given\n");
		usage(&opt);
	}
	if (opt.bam_ofn && cgsum == 0) {
		fprintf(stderr, "-o FILE specified, but no -C/-G\n");
		usage(&opt);
	}

	if (opt.fwd_only && opt.rev_only) {
		fprintf(stderr, "-f and -r flags are mutually incompatible\n");
		usage(&opt);
	}

	if (argc-optind != 2) {
		usage(&opt);
	}

	opt.bam_fn = argv[optind];
	opt.fasta_fn = argv[optind+1];

	return (condamage(&opt) < 0);
}

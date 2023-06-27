import pysam
import glob
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("inp", help = "input file created from x.vcf by grep -v '#' x.vcf | cut -f1,2,4,5 ")
parser.add_argument("bam_pattern", help = "path and pattern of bam files to analyze")
parser.add_argument("out", help = "output table name")
args = parser.parse_args()

out = open(args.out, "w", buffering = 1)
p = sorted(glob.glob(args.bam_pattern))
pp = [x.split("/")[-1].split("_")[0] for x in p]
out.write("Chrom" + " " + "Pos" + " " + "Ref" + " " + "Alt" + " " + " ".join([x + "_Supp" for x in pp]) + " " + " ".join([x + "_Cov" for x in pp]) + "\n")

with open(args.inp, "r") as f:
	for line in f:
#		a,c,mref,malt = [],[],[],[]
		a = [0] * len(p)
		c = [0] * len(p)
		chrom,pos,ref,alt = line.strip().split("\t")
		if "," in alt:
			if "*" in alt:
				alt = alt[:-2]
			else:
				alt = alt.split(",")[0]
		if len(ref) == len(alt):
			if len(alt) > 1:
				alt = alt[0]
			for b in p:
				x = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
#				m = []
				nn = b.split("/")[-1].split("_")[0]
				bamfile = pysam.AlignmentFile(b, "rb")
				for pu in bamfile.pileup(chrom, int(pos), int(pos)+1, stepper = "nofilter"):
					if pu.pos == int(pos) - 1:
						for read in pu.pileups:
							if not read.is_del and not read.is_refskip:
								x[read.alignment.query_sequence[read.query_position]] += 1
#								m.append(read.alignment.mapping_quality)
				a[p.index(b)] = str(x[alt])
				c[p.index(b)] = str(sum(x.values()))
#				mref.append(str(np.mean([y for y,z in zip(m, x) if z == ref])))
#				if len([y for y in x if y == alt]) > 0:
#					malt.append(str(round(np.mean([y for y,z in zip(m, x) if z == alt]))))
#				else:
#					malt.append("NA")
				bamfile.close()
#			print(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + " " + " ".join(mref) + " " + " ".join(malt))
#			out.write(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + " " + " ".join(mref) + " " + " ".join(malt) + "\n")
#			print(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + "\n")
			out.write(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + "\n")
		elif len(ref) > len(alt):
			for b in p:
#				x, m = [], []
				x = {"a": 0, "b": 0}
				bamfile = pysam.AlignmentFile(b, "rb")
				for pu in bamfile.pileup(chrom, int(pos), int(pos)+1, stepper = "nofilter"):
					if pu.pos == int(pos):
						for read in pu.pileups:
							if not read.is_refskip:
								if read.is_del:
									x["a"] += 1
								else:
									x["b"] += 1
#								m.append(read.alignment.mapping_quality)
				a[p.index(b)] = str(x["a"])
				c[p.index(b)] = str(x["a"] + x["b"])

#				mref.append(str(np.mean([y for y,z in zip(m, x) if z == "+"])))
#				if len([y for y in x if y == "-"]) > 0:
#					malt.append(str(round(np.mean([y for y,z in zip(m, x) if z == "-"]))))
#				else:
#					malt.append("NA")
				bamfile.close()
#			print(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + " " + " ".join(mref) + " " + " ".join(malt))
#			out.write(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + " " + " ".join(mref) + " " + " ".join(malt) + "\n")
#			print(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + "\n")
			out.write(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + "\n")
		else:
			for b in p:
#				x, m = [],[]
				x = {alt[1:]: 0}
				nn = b.split("/")[-1].split("_")[0]
				bamfile = pysam.AlignmentFile(b, "rb")
				for pu in bamfile.pileup(chrom, int(pos), int(pos)+1, stepper = "nofilter"):
					if pu.pos == int(pos):
						for read in pu.pileups:
							if not read.is_refskip:
								if read.query_position is not None:
									try:
										x[read.alignment.query_sequence[read.query_position:(read.query_position + len(alt) - 1)]] += 1
									except KeyError:
										x[read.alignment.query_sequence[read.query_position:(read.query_position + len(alt) - 1)]] = 1
#									m.append(read.alignment.mapping_quality)
				a[p.index(b)] = str(x[alt[1:]])
				c[p.index(b)] = str(sum(x.values()))
#				mref.append(str(np.mean([y for y,z in zip(m, x) if z != alt[1:]])))
#				if len([y for y in x if y == alt[1:]]) > 0:
#					malt.append(str(round(np.mean([y for y,z in zip(m, x) if z == alt[1:]]))))
#				else:
#					malt.append("NA")
				bamfile.close()
#			print(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + " " + " ".join(mref) + " " + " ".join(malt))
#			out.write(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + " " + " ".join(mref) + " " + " ".join(malt) + "\n")
#			print(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + "\n")
			out.write(chrom + " " + pos + " " + ref + " " + alt + " "  + " ".join(a) + " " + " ".join(c) + "\n")
out.close

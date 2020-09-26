# SiteStripper

SiteStripper is a Perl program that removes fast-evolving sites from alignments. In doing this, it takes data partitioning and secondary structure into account so that you can still run phylogenetic analyses with partitioned models or a doublet model on the output alignments.

SiteStripper does not calculate site rates. You need to use other software such as HyPhy or IQ-Tree to do this.

### User guide

SiteStripper is a command-line program. Just run the program and it will show you a list of command-line flags you can use to control the program (see also below). The alignment file has to be in fasta format. The rates file has to be in the HyPhy or IQ-tree format. See the example files for more details.

SiteStripper is designed to keep track of which partition (character set) any give site is in. After removing the fast-evolving sites, it will return the character set definitions that allow you to run phylogenetic analyses with the same partitioning strategy as before but with a subset of the characters. If you want to use this option, you need to make a file in which the partitions are defined. See the parts.txt file of the example for formatting. To run the example with partitions, type 

```perl sitestripper.pl -a alignment.fas -r siterates.txt -f .8 -o out.nex -pt parts.txt```

If you have paired bases (e.g. stem regions in rRNA genes) you can specify this and the paired bases will either be retained together or discarded together. To determine whether a pair is discarded or retained, the mean of its base rates is calculated and compared to the rates of other bases (and base pairs). If you want to use this option, you need a file in which you specify base pairings as follows:

```
1:12
2:11
3:10
4:9
18:28
19:27
20:26
```

If you want to see whether removing fast sites gives you better results than removing random sites, you can have the program generate alignments from which randomly picked sites were removed. See https://doi.org/10.1093/molbev/msq091 for an example.

Here is a complete list of the command-line flags you can use.

```
mandatory parameters
   -a    alignment_file  (must be in fasta format)
   -r    rates_file (not mandatory if you use -ra flag)
   -f    fraction of characters to keep
   -o    output alignment  (in nexus format)

optional parameters
   -rf   rates file format (hyphy or iqtree, default hyphy)
   -pa   file with list of paired bases
   -pt   file with definitions of data partitions  (one per line)
   -of   output format: fasta|nexus  (default: nexus)
   -ra   number of alignments with random sites stripped to be generated
```		

### Citation
If you find these tools useful, please cite them in your work. I recommend citing them as follows:
Verbruggen H. (2018) SiteStripper version 1.02. https://github.com/hverbruggen/SiteStripper

### Notes and disclaimer
SiteStripper is in development and has not been tested extensively. It is quite plausible that incorrectly formatted input could lead to nonsensical output. In such cases, you should double-check your input, compare it to the example files and try again. If this still doesn't work, please feel free to write me an email (heroen.verbruggen@gmail.com).
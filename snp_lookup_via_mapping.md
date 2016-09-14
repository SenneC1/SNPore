

    # Init
    #
    
    snpFile   = '/home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta'  # REMOVE FIRST LINE IN ORIGINAL FILE (BREAKS FASTA FORMAT)
    readFile  = '/home/senne/nanopore/SNP/Nanopore_data/2_potential_snp_amplicons_3mism.fasta'
    resultDir = '/home/senne/nanopore/SNP/results_yannick'
    fastq_file_name= '/home/senne/nanopore/SNP/Nanopore_data/ligatedSNPs.fastq'
    bwa       = '/opt/tools/bwa-0.7.15'            # v0.7.5
    samtools  = '/opt/tools/samtools-1.3.1' # v1.3.1
    bcftools  = '/opt/tools/bcftools-1.3.1' # v1.3.1
    
    # Check
    !ls {snpFile}
    print('Number of SNPs in reference file:')
    !grep -c ">" {snpFile}


    /home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta
    Number of SNPs in reference file:
    52


## Info on the fastq file


    #from collections import Counter
    #fastq_file_name = 'ligatedSTRs1.fastq'
    #histogram_data = Counter()
    hist_array = []
    hist_arrayG = []
    hist_arrayA = []
    hist_arrayC = []
    hist_arrayT = []
    with open(fastq_file_name, "rb") as infile:
        for line in infile:
            if line.startswith(b'G'):
                read_length = len(line[:-1]) #Last char is \n
                #histogram_data[read_length] += 1
                hist_arrayG.append(read_length)
    with open(fastq_file_name, "rb") as infile:
        for line in infile:
            if line.startswith(b'A'):
                read_length = len(line[:-1]) #Last char is \n
                #histogram_data[read_length] += 1
                hist_arrayA.append(read_length) 
    with open(fastq_file_name, "rb") as infile:
        for line in infile:
            if line.startswith(b'C'):
                read_length = len(line[:-1]) #Last char is \n
                #histogram_data[read_length] += 1
                hist_arrayC.append(read_length) 
    with open(fastq_file_name, "rb") as infile:
        for line in infile:
            if line.startswith(b'T'):
                read_length = len(line[:-1]) #Last char is \n
                #histogram_data[read_length] += 1
                hist_arrayT.append(read_length)
    
    hist_array.extend(hist_arrayA)
    hist_array.extend(hist_arrayC)
    hist_array.extend(hist_arrayG)
    hist_array.extend(hist_arrayT)
    
    import numpy as np
    
    len (hist_array), np.mean(hist_array), np.median(hist_array) , max(hist_array), np.std(hist_array)




    (14324, 507.74413571628037, 411.0, 5850, 324.5996162457601)




    # Map reads to reference sequences
    #
    
    # Build index of the references
    !{bwa} index {snpFile}
    
    # Map reads
    !{bwa} mem -t 10 -k 14 -W 20 -r 10 -O 6,6 -E 1 {snpFile} {readFile} > {resultDir}/ONTGAP4ttt.sam


    [bwa_index] Pack FASTA... 0.00 sec
    [bwa_index] Construct BWT for the packed sequence...
    [bwa_index] 0.00 seconds elapse.
    [bwa_index] Update BWT... 0.00 sec
    [bwa_index] Pack forward-only FASTA... 0.00 sec
    [bwa_index] Construct SA from BWT and Occ... 0.00 sec
    [main] Version: 0.7.15-r1140
    [main] CMD: /opt/tools/bwa-0.7.15 index /home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta
    [main] Real time: 0.016 sec; CPU: 0.006 sec
    [M::bwa_idx_load_from_disk] read 0 ALT contigs
    [M::process] read 23057 sequences (1988758 bp)...
    [M::mem_process_seqs] Processed 23057 reads in 1.265 CPU sec, 0.141 real sec
    [main] Version: 0.7.15-r1140
    [main] CMD: /opt/tools/bwa-0.7.15 mem -t 10 -k 14 -W 20 -r 10 -O 6,6 -E 1 /home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta /home/senne/nanopore/SNP/Nanopore_data/2_potential_snp_amplicons_3mism.fasta
    [main] Real time: 0.195 sec; CPU: 1.311 sec



    # Make sorted bam and index
    #
    
    !{samtools} view -Sbu {resultDir}/ONTGAP4ttt.sam | {samtools} sort -o {resultDir}/ONTGAP4ttt_sorted.bam -
    !{samtools} index {resultDir}/ONTGAP4ttt_sorted.bam {resultDir}/ONTGAP4ttt_sorted.bam.bai
    print('Done')


    Done



    # Display some stats
    #
    
    !{samtools} flagstat {resultDir}/ONTGAP4ttt_sorted.bam


    23062 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    5 + 0 supplementary
    0 + 0 duplicates
    8034 + 0 mapped (34.84% : N/A)
    0 + 0 paired in sequencing
    0 + 0 read1
    0 + 0 read2
    0 + 0 properly paired (N/A : N/A)
    0 + 0 with itself and mate mapped
    0 + 0 singletons (N/A : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)



    # Generate vcf file from bam file. Needs the reference and its index file 
    #
    # Note: the commands below are for samtools and bcftools v1.3.1 (will not work on v0.1.19!)
    
    # Reporting all positions
    !{samtools} mpileup -d 100000 -uf {snpFile} {resultDir}/ONTGAP4ttt_sorted.bam | {bcftools} call -V indels -m - > {resultDir}/ONTGAP4ttt_sorted.bam.vcf
    
    # Reporting variants only (excludes SNPs homozygous for reference allele)
    !{samtools} mpileup -d 100000 -uf {snpFile} {resultDir}/ONTGAP4ttt_sorted.bam | {bcftools} call -V indels -mv - > {resultDir}/ONTGAP4yyy_sorted.bam.vcf
    
    print('Done')

    Note: Neither --ploidy nor --ploidy-file given, assuming all sites are diploid
    [mpileup] 1 samples in 1 input files
    Note: Neither --ploidy nor --ploidy-file given, assuming all sites are diploid
    [mpileup] 1 samples in 1 input files
    Done



    # Check the vcf file
    #
    
    !head -n 100 {resultDir}/ONTGAP2ttt_sorted.bam.vcf


    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##samtoolsVersion=1.3.1+htslib-1.3.1
    ##samtoolsCommand=samtools mpileup -d 100000 -uf /home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta /home/senne/nanopore/SNP/results_yannick/ONTGAP2ttt_sorted.bam
    ##reference=file:///home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta
    ##contig=<ID=rs1490413,length=51>
    ##contig=<ID=rs876724,length=51>
    ##contig=<ID=rs1357617,length=51>
    ##contig=<ID=rs2046361,length=51>
    ##contig=<ID=rs717302,length=51>
    ##contig=<ID=rs1029047,length=51>
    ##contig=<ID=rs917118,length=51>
    ##contig=<ID=rs763869,length=51>
    ##contig=<ID=rs1015250,length=51>
    ##contig=<ID=rs735155,length=51>
    ##contig=<ID=rs901398,length=51>
    ##contig=<ID=rs2107612,length=51>
    ##contig=<ID=rs1886510,length=51>
    ##contig=<ID=rs1454361,length=51>
    ##contig=<ID=rs2016276,length=51>
    ##contig=<ID=rs729172,length=51>
    ##contig=<ID=rs740910,length=51>
    ##contig=<ID=rs1493232,length=51>
    ##contig=<ID=rs719366,length=51>
    ##contig=<ID=rs1031825,length=51>
    ##contig=<ID=rs722098,length=51>
    ##contig=<ID=rs733164,length=51>
    ##contig=<ID=rs826472,length=51>
    ##contig=<ID=rs2831700,length=51>
    ##contig=<ID=rs873196,length=51>
    ##contig=<ID=rs1382387,length=51>
    ##contig=<ID=rs2111980,length=51>
    ##contig=<ID=rs2056277,length=51>
    ##contig=<ID=rs1024116,length=51>
    ##contig=<ID=rs727811,length=51>
    ##contig=<ID=rs1413212,length=51>
    ##contig=<ID=rs938283,length=51>
    ##contig=<ID=rs1979255,length=51>
    ##contig=<ID=rs1463729,length=51>
    ##contig=<ID=rs2076848,length=51>
    ##contig=<ID=rs1355366,length=51>
    ##contig=<ID=rs907100,length=51>
    ##contig=<ID=rs354439,length=51>
    ##contig=<ID=rs2040411,length=51>
    ##contig=<ID=rs737681,length=51>
    ##contig=<ID=rs2830795,length=51>
    ##contig=<ID=rs251934,length=51>
    ##contig=<ID=rs914165,length=51>
    ##contig=<ID=rs10495407,length=51>
    ##contig=<ID=rs1360288,length=51>
    ##contig=<ID=rs964681,length=51>
    ##contig=<ID=rs1005533,length=51>
    ##contig=<ID=rs8037429,length=51>
    ##contig=<ID=rs891700,length=51>
    ##contig=<ID=rs1335873,length=51>
    ##contig=<ID=rs1028528,length=51>
    ##contig=<ID=rs1528460,length=51>
    ##ALT=<ID=*,Description="Represents allele(s) other than observed.">
    ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
    ##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
    ##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
    ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
    ##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
    ##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
    ##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
    ##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
    ##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
    ##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
    ##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
    ##bcftools_callVersion=1.3.1+htslib-1.3.1
    ##bcftools_callCommand=call -V indels -m -
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	/home/senne/nanopore/SNP/results_yannick/ONTGAP2ttt_sorted.bam
    rs1490413	1	.	A	.	0	.	DP=121;MQSB=0.80569;MQ0F=0;AN=2;DP4=44,77,0,0;MQ=55	GT	0/0
    rs1490413	2	.	C	.	0	.	DP=124;MQSB=0.813884;MQ0F=0;AN=2;DP4=45,79,0,0;MQ=55	GT	0/0
    rs1490413	3	.	T	.	0	.	DP=127;MQSB=0.861099;MQ0F=0;AN=2;DP4=45,82,0,0;MQ=55	GT	0/0
    rs1490413	4	.	G	.	0	.	DP=128;MQSB=0.875116;MQ0F=0;AN=2;DP4=45,83,0,0;MQ=55	GT	0/0
    rs1490413	5	.	G	.	0	.	DP=135;MQSB=0.954599;MQ0F=0;AN=2;DP4=47,88,0,0;MQ=56	GT	0/0
    rs1490413	6	.	G	.	0	.	DP=136;MQSB=0.939625;MQ0F=0;AN=2;DP4=48,88,0,0;MQ=56	GT	0/0
    rs1490413	7	.	C	.	0	.	DP=138;SGB=-0.379885;RPB=1;MQB=1;MQSB=0.956775;BQB=1;MQ0F=0;AN=2;DP4=48,89,0,1;MQ=56	GT	0/0
    rs1490413	8	.	T	.	0	.	DP=142;MQSB=0.981754;MQ0F=0;AN=2;DP4=48,94,0,0;MQ=56	GT	0/0
    rs1490413	9	.	G	.	30.0896	.	DP=143;VDB=4.54206e-05;SGB=-0.651104;RPB=0.845307;MQB=0.909799;MQSB=0.986198;BQB=1;MQ0F=0;AN=2;DP4=48,87,0,8;MQ=56	GT	0/0
    rs1490413	10	.	A	.	0	.	DP=146;SGB=-0.379885;RPB=1;MQB=1;MQSB=0.995617;BQB=1;MQ0F=0;AN=2;DP4=48,97,0,1;MQ=56	GT	0/0
    rs1490413	11	.	T	.	0	.	DP=151;MQSB=0.993448;MQ0F=0;AN=2;DP4=50,101,0,0;MQ=56	GT	0/0
    rs1490413	12	.	G	.	0	.	DP=152;VDB=0.02;SGB=-0.453602;RPB=0.846667;MQB=0.773333;MQSB=0.98759;BQB=1;MQ0F=0;AN=2;DP4=51,99,0,2;MQ=56	GT	0/0
    rs1490413	13	.	T	.	0	.	DP=154;MQSB=0.984343;MQ0F=0;AN=2;DP4=52,102,0,0;MQ=56	GT	0/0
    rs1490413	14	.	G	.	0	.	DP=155;VDB=0.0433485;SGB=-0.511536;RPB=0.315182;MQB=0.803014;MQSB=0.988001;BQB=1;MQ0F=0;AN=2;DP4=50,102,2,1;MQ=56	GT	0/0
    rs1490413	15	.	G	.	0	.	DP=160;SGB=-0.379885;RPB=1;MQB=1;MQSB=0.985418;BQB=1;MQ0F=0;AN=2;DP4=53,106,1,0;MQ=56	GT	0/0
    rs1490413	16	.	G	.	0	.	DP=162;SGB=-0.379885;RPB=1;MQB=1;MQSB=0.982256;BQB=1;MQ0F=0;AN=2;DP4=54,107,1,0;MQ=56	GT	0/0
    rs1490413	17	.	T	.	0	.	DP=167;SGB=-0.379885;RPB=1;MQB=1;MQSB=0.989433;BQB=1;MQ0F=0;AN=2;DP4=55,111,1,0;MQ=56	GT	0/0
    rs1490413	18	.	T	.	0	.	DP=168;MQSB=0.992071;MQ0F=0;AN=2;DP4=56,112,0,0;MQ=56	GT	0/0
    rs1490413	19	.	C	.	0	.	DP=168;MQSB=0.992071;MQ0F=0;AN=2;DP4=56,112,0,0;MQ=56	GT	0/0
    rs1490413	20	.	T	.	0	.	DP=167;MQSB=0.999173;MQ0F=0;AN=2;DP4=56,111,0,0;MQ=57	GT	0/0



    # Get SNP profile
    
    
    snpData = {}
    
    with open(resultDir+'/ONTGAP4ttt_sorted.bam.vcf') as f:
        for l in f:
            if l.startswith('#'):
                continue
                
            snp, pos, id, ref, alt, qual, filter, info, d, dd = l.split()
            
            # Our SNP of interest is always at position 26 of the reference
            if int(pos) != 26:
                continue
    
            par = {}
            for p in info.split(';'):
                pv = p.split('=')
                par[pv[0]] = pv[1]
            
            snpData[snp] = {'pos': pos, 'ref': ref, 'alt': alt, 'qual': qual, 'filter': filter, 'info': par}
    
    # DEBUG
    print('Got data for {} SNPs:'.format(len(snpData)))
    
    # Save/print results
    with open(resultDir + '/ttt_profile.csv', 'w') as f:
        # Table header
        f.write('snp, coverage, ref_allele, ref_percent, alt_allele, alt_percent, genotype\n')
        
        # Table data
        for s in sorted(snpData.keys()):
            totalDepth = int(snpData[s]['info']['DP'])
            depthList  = [int(d) for d in snpData[s]['info']['DP4'].split(',')]
            refDepth   = sum(depthList[0:2])
            altDepth   = sum(depthList[2:4])
            
            # Estimate the diploid genotype: when the minor allele is more than 10 times weaker than the major allele,
            # we should ignore it for a pure sample?
            if refDepth > altDepth and altDepth/refDepth < 0.1:
                genotype = snpData[s]['ref'] + snpData[s]['ref']
            elif altDepth > refDepth and refDepth/altDepth < 0.1:
                genotype = snpData[s]['alt'] + snpData[s]['alt']
            else:
                genotype = snpData[s]['ref'] + snpData[s]['alt']
            
            if snpData[s]['alt'] == '.':
                # Only 1 allele was observed
                f.write(','.join([s, str(totalDepth), snpData[s]['ref'], '{:.1f}'.format(100*refDepth/totalDepth), '', '', snpData[s]['ref']+snpData[s]['ref']]) + '\n')
                # DEBUG
                print('  {} ({})  {} ({:.1f} %)'.format(s, totalDepth, snpData[s]['ref'], 100*refDepth/totalDepth))
            else:
                # Two alleles were observed
                f.write(','.join([s, str(totalDepth), snpData[s]['ref'], '{:.1f}'.format(100*refDepth/totalDepth), snpData[s]['alt'], '{:.1f}'.format(100*altDepth/totalDepth), genotype]) + '\n')
                # DEBUG
                print('  {} ({})  {} ({:.1f} %)  {} ({:.1f} %)'.format(s, totalDepth, snpData[s]['ref'], 100*refDepth/totalDepth, snpData[s]['alt'], 100*altDepth/totalDepth))


    Got data for 52 SNPs:
      rs1005533 (99)  A (61.6 %)  G (38.4 %)
      rs1015250 (67)  C (3.0 %)  G (97.0 %)
      rs1024116 (111)  A (38.7 %)  G (61.3 %)
      rs1028528 (136)  A (2.9 %)  G (97.1 %)
      rs1029047 (24)  A (100.0 %)
      rs1031825 (82)  A (23.2 %)  C (76.8 %)
      rs10495407 (188)  A (60.6 %)  G (39.4 %)
      rs1335873 (120)  A (1.7 %)  T (98.3 %)
      rs1355366 (166)  A (5.4 %)  G (94.6 %)
      rs1357617 (220)  A (100.0 %)
      rs1360288 (184)  C (59.8 %)  T (40.2 %)
      rs1382387 (207)  G (2.4 %)  T (97.6 %)
      rs1413212 (168)  A (64.9 %)  G (35.1 %)
      rs1454361 (200)  A (100.0 %)
      rs1463729 (157)  A (97.5 %)
      rs1490413 (235)  A (98.3 %)
      rs1493232 (78)  A (67.9 %)  C (32.1 %)
      rs1528460 (85)  C (2.4 %)  T (97.6 %)
      rs1886510 (144)  C (4.9 %)  T (95.1 %)
      rs1979255 (301)  C (98.0 %)
      rs2016276 (70)  A (95.7 %)
      rs2040411 (72)  A (0.0 %)  G (100.0 %)
      rs2046361 (111)  A (3.6 %)  T (96.4 %)
      rs2056277 (199)  C (47.7 %)  T (52.3 %)
      rs2076848 (127)  A (72.4 %)  T (27.6 %)
      rs2107612 (49)  A (0.0 %)  G (100.0 %)
      rs2111980 (212)  A (61.3 %)  G (38.7 %)
      rs251934 (136)  C (2.2 %)  T (97.8 %)
      rs2830795 (74)  A (74.3 %)  G (25.7 %)
      rs2831700 (299)  A (99.0 %)
      rs354439 (128)  A (73.4 %)  T (26.6 %)
      rs717302 (113)  A (9.7 %)  G (90.3 %)
      rs719366 (49)  C (53.1 %)  T (46.9 %)
      rs722098 (232)  A (99.1 %)
      rs727811 (113)  A (6.2 %)  C (93.8 %)
      rs729172 (231)  A (65.4 %)  C (34.6 %)
      rs733164 (177)  A (21.5 %)  G (78.5 %)
      rs735155 (159)  A (6.3 %)  G (93.7 %)
      rs737681 (237)  C (100.0 %)
      rs740910 (22)  A (0.0 %)  G (100.0 %)
      rs763869 (163)  C (61.3 %)  T (38.7 %)
      rs8037429 (69)  C (58.0 %)  T (42.0 %)
      rs826472 (144)  C (100.0 %)
      rs873196 (149)  C (17.4 %)  T (82.6 %)
      rs876724 (118)  C (99.2 %)
      rs891700 (39)  A (64.1 %)  G (35.9 %)
      rs901398 (204)  C (56.9 %)  T (43.1 %)
      rs907100 (150)  C (95.3 %)
      rs914165 (165)  A (6.7 %)  G (93.3 %)
      rs917118 (407)  C (99.8 %)
      rs938283 (101)  C (71.3 %)  T (28.7 %)
      rs964681 (81)  C (3.7 %)  T,G (96.3 %)


## Locus count in fasta


    %%bash
    
    cd /home/senne/nanopore/SNP/Nanopore_data/
    
    grep -c 'rs2056277' 2_potential_snp_amplicons_3mism.fasta
    grep -c 'rs1413212' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs2107612' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs2111980' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs251934' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1028528' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs2831700' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs901398' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs722098' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs2076848' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs1493232' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs735155' 2_potential_snp_amplicons_3mism.fasta 
     grep -c 'rs1528460' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs1005533' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs733164' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1029047' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs727811' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1024116' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs1015250' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs907100' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs737681' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs717302' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs740910' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs1979255' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs719366' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs2040411' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs2016276' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs135761' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs8037429' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1335873' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs876724' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs354439' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1031825' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs873196' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1463729' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1886510' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs763869' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs2830795' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1454361' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1355366' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs938283' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1490413' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs964681' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs826472' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs729172' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1382387' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs10495407' 2_potential_snp_amplicons_3mism.fasta 
      grep -c 'rs891700' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs917118' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs914165' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs2046361' 2_potential_snp_amplicons_3mism.fasta
      grep -c 'rs1360288' 2_potential_snp_amplicons_3mism.fasta


    399
    592
    246
    783
    407
    603
    521
    478
    424
    523
    264
    512
    337
    396
    676
    606
    246
    383
    490
    443
    485
    269
    228
    585
    252
    456
    372
    422
    196
    504
    306
    463
    301
    557
    463
    559
    372
    293
    419
    687
    389
    488
    392
    304
    598
    589
    518
    197
    606
    623
    356
    479




# Plot mapping


    import numpy as np
    import matplotlib.pyplot as plt
    
    with open("/home/senne/nanopore/SNP/test") as f:
        data = f.read()
    
    data = data.split('\n')
    print (data)
    x = [row.split(' ')[0] for row in data]
    y = [row.split(' ')[1] for row in data]
    
    
    fig = plt.figure()

    ['1 55,1', '2 10,7', '3 38,2', '4 4,2', '5 78,1', '6 23,9', '7 59,4', '8 5,1', '9 7,4', '10 95,6', '11 53,3', '12 3', '13 63,7', '14 99,1', '15 94,9', '16 95,7', '17 69,5', '18 8,4', '19 8,7', '20 93,8', '21 85,5', '22 1,6', '23 6,2', '24 48', '25 66,4', '26 2,3', '27 55,5', '28 4,6', '29 66,3', '30 96,6', '31 56,8', '32 11,7', '33 54,7', '34 97', '35 9,7', '36 61,3', '37 28,4', '38 7,8', '39 96,2', '40 2,3', '41 58,1', '42 55,5', '43 97,7', '44 23,4', '45 96,5', '46 61,6', '47 56,4', '48 90,5', '49 9,1', '50 98,7', '51 65', '52 14,2', '']



    ---------------------------------------------------------------------------

    IndexError                                Traceback (most recent call last)

    <ipython-input-52-d21dcd508237> in <module>()
          8 print (data)
          9 x = [row.split(' ')[0] for row in data]
    ---> 10 y = [row.split(' ')[1] for row in data]
         11 
         12 


    <ipython-input-52-d21dcd508237> in <listcomp>(.0)
          8 print (data)
          9 x = [row.split(' ')[0] for row in data]
    ---> 10 y = [row.split(' ')[1] for row in data]
         11 
         12 


    IndexError: list index out of range



    import numpy as np
    import matplotlib.pyplot as plt
    
    
    
    N = 5
    
    Allele1 = (4,5,3,2,1)
    
    Allele2 = (1,2,3,5,4)
    
    #Other= (2,25,25,25,2)
    
    width = 0.25
    ind = np.arange(N)  
    
    rects1 = plt.bar(ind, Allele1, width, color='r')
    
    rects2 = plt.bar(ind+width, Allele2, width, color='b')
    
    #rects3 = plt.bar(ind+width+width, Other, width, color='y')
    
    
    
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Relative frequency')
    #ax.set_title('Scores by group and gender')
    ax.set_xticks(ind)
    ax.set_xticklabels(('a,b,c,d,g'))
    ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))
    
    plt.show()


    

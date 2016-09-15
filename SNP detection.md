#SNP Detection

## Setting the dependencies
```python
    snpFile   = '/home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta'
    readFile  = '/home/senne/nanopore/SNP/Nanopore_data/2_potential_snp_amplicons_3mism.fasta'
    resultDir = '/home/senne/nanopore/SNP/results_yannick'
    fastqFile = '/home/senne/nanopore/SNP/Nanopore_data/ligatedSNPs.fastq'
    bwa       = '/opt/tools/bwa-0.7.15'     # v0.7.5
    samtools  = '/opt/tools/samtools-1.3.1' # v1.3.1
    bcftools  = '/opt/tools/bcftools-1.3.1' # v1.3.1
```
```python
    # Check
    !ls {snpFile}
    print('Number of SNPs in reference file:')
    !grep -c ">" {snpFile}
```

    /home/senne/nanopore/SNP/known_SNP_sequence/SNP_sequence.fasta
    Number of SNPs in reference file:
    52


## Info on the fastq file
```python
    hist_array = []
    hist_arrayG = []
    hist_arrayA = []
    hist_arrayC = []
    hist_arrayT = []
    with open(fastqFile, "rb") as infile:
        for line in infile:
            if line.startswith(b'G'):
                read_length = len(line[:-1]) #Last char is \n
                #histogram_data[read_length] += 1
                hist_arrayG.append(read_length)
    with open(fastqFile, "rb") as infile:
        for line in infile:
            if line.startswith(b'A'):
                read_length = len(line[:-1]) #Last char is \n
                #histogram_data[read_length] += 1
                hist_arrayA.append(read_length) 
    with open(fastqFile, "rb") as infile:
        for line in infile:
            if line.startswith(b'C'):
                read_length = len(line[:-1]) #Last char is \n
                #histogram_data[read_length] += 1
                hist_arrayC.append(read_length) 
    with open(fastqFile, "rb") as infile:
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
```

## Map reads to reference sequences
```    
    # Build index of the references
    !{bwa} index {snpFile}
    
    # Map reads
    !{bwa} mem -t 10 -k 14 -W 20 -r 10 -O 6,6 -E 1 {snpFile} {readFile} > {resultDir}/ONTGAP4ttt.sam

    # Make sorted bam and index
    
    !{samtools} view -Sbu {resultDir}/ONTGAP4ttt.sam | {samtools} sort -o {resultDir}/ONTGAP4ttt_sorted.bam -
    !{samtools} index {resultDir}/ONTGAP4ttt_sorted.bam {resultDir}/ONTGAP4ttt_sorted.bam.bai
    print('Done')

    # Display some stats
    
    !{samtools} flagstat {resultDir}/ONTGAP4ttt_sorted.bam
```

 ## Generate vcf file from bam file. Needs the reference and its index file 
```    
    # Note: the commands below are for samtools and bcftools v1.3.1
    
    # Reporting all positions
    !{samtools} mpileup -d 100000 -uf {snpFile} {resultDir}/ONTGAP4ttt_sorted.bam | {bcftools} call -V indels -m - > {resultDir}/ONTGAP4ttt_sorted.bam.vcf
    
    print('Done')
```

## Check the content of the vcf file
```
    !head -n 100 {resultDir}/ONTGAP2ttt_sorted.bam.vcf
```

## Get SNP profile
```python
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
```
```
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
```

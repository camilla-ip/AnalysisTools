#!/usr/bin/env python

# ------------------------------------------------------------------------------
# splitpops
#
# Split a file of Illumina paired-end reads for a Hepatitis C Virus (HCV) sample
# into one BAM file per sub-genotype population. Only sub-genotypes that attract
# at least M% of the read pairs (default 1%) are reported as a population.
#
# Input
# - Sample reads (BAM, paired-end)
# - HCV target references (FASTA, each reference named SUBGENOTYPE_REFID)
# - The minimum percentage of sample reads that match a target reference
#   sub-genotype (range 0 to 100, default 1.0 for 1%)
#
# Method
# - Map each of the sample reads to the target references of known sub-genotype,
#   discarding any non-primary or supplementary mappings.
# - If  at least one read of the pair were able to be mapped to a target
#   reference genome, the pair was allocated the sub-genotype of the highest
#   scoring match. If both mapped to different sub-genotypes, the allocated
#   sub-genotype read pair was marked as "ambiguous". If neither read mapped,
#   the sub-genotype was designated "unclassified".
#
# Output
# - Summary stats on each population (DATAID.splitpops.stats.txt file, containing
#   tab-separated fields "dataid, popN, readclass, readcnt, readpct, bampath"
#   where N is a number starting from 1, readclass is the sub-genotype (e.g., 1a)
# - A set of BAM files, one per sub-genotype (e.g., DATAID.splitpops.pop1.bam)
# 
# Camilla Ip
# The Wellcome Centre for Human Genetics
# University of Oxford
# ------------------------------------------------------------------------------

import os
import sys
import time

import snorkversion

_processname = 'splitpops'
_outsuffix = 'splitpops'
_outfilesuffixD = {
    'readclassbam'   : { 'suffix':'.'+_outsuffix+'.readclass.bam', 'istmp':False, 'isfinal':True,  'msg':'Reads tagged with RG:Z:readclass (BAM)' },
    'statspath'      : { 'suffix':'.'+_outsuffix+'.stats.txt',     'istmp':False, 'isfinal':True,  'msg':'Output statistics for target pops (TSV with header)' },
    'sentinelpath'   : { 'suffix':'.'+_outsuffix+'.sentinel.txt',  'istmp':False, 'isfinal':True,  'msg':'Date and time of last computation (TXT)' },
    'fq1'            : { 'suffix':'_1.fastq',                      'istmp':True,  'isfinal':False, 'msg':'Forward reads in name sort order (FASTQ)' },
    'fq2'            : { 'suffix':'_2.fastq',                      'istmp':True,  'isfinal':False, 'msg':'Reverse reads in name sort order (FASTQ)' },
    'tmpbamstempath' : { 'suffix':'_namesorted_tmp',               'istmp':True,  'isfinal':False, 'msg':'Reads in name sort order (BAM)' },
    'mapbam'         : { 'suffix':'.'+_outsuffix+'.mapped.bam',    'istmp':True,  'isfinal':False, 'msg':'Reads mapped to refs in name sort order (BAM)' },
    'mapbamstem'     : { 'suffix':'.'+_outsuffix+'.mapped',        'istmp':True,  'isfinal':False, 'msg':'Reads mapped to refs in name sort order (BAM stem)' },
    'mapsam'         : { 'suffix':'.'+_outsuffix+'.mapped.sam',    'istmp':True,  'isfinal':False, 'msg':'Reads mapped to refs in name sort order (SAM)' },
    'readclasspath'  : { 'suffix':'.'+_outsuffix+'.readclass.txt', 'istmp':True,  'isfinal':False, 'msg':'Readclass of each read pair: readid readclass (TSV with header)' },
    'readclasssam'   : { 'suffix':'.'+_outsuffix+'.readclass.sam', 'istmp':True,  'isfinal':False, 'msg':'Reads tagged with RG:Z:readclass (SAM)' },
    'popNsam'        : { 'suffix':'.'+_outsuffix+'.popN.sam',      'istmp':True,  'isfinal':False, 'msg':'Reads from target population N (SAM)' },
    'popNbam'        : { 'suffix':'.'+_outsuffix+'.popN.bam',      'istmp':False, 'isfinal':True , 'msg':'Reads from target population N (BAM)' },
    'samsubsetpath'  : { 'suffix':'.'+_outsuffix+'.mapbam.subset', 'istmp':True,  'isfinal':False, 'msg':'Subset of mapbam columns (TSV)' }
}

_didsomecomputation = False	# Used to determine when to generate the sentinel.txt file

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Check pre-conditions.'
    mylogger.debug('Checking pre-requisites')
    mylogger.debug('Checking pre-requisites finished')
    return 0

def File_Path(outdir, dataid, filekey, N=0):
    'Return path for file type specified by filekey.'
    filepath = os.path.join(outdir, dataid+_outfilesuffixD[filekey]['suffix'])
    if filekey in ['popNsam', 'popNbam'] and N:
        filepath = filepath.replace('popN', 'pop{N}'.format(N=N))
    return filepath

def Convert_Input_Reads_BAM_To_FASTQ(args, P, mylogger, myhandler, _processname):
    'Name-sort the BAM file, then convert to a pair of FASTQ files. Exit program on error.'
  # Output files
    fq1 = File_Path(args.outdir, args.dataid, 'fq1')
    fq2 = File_Path(args.outdir, args.dataid, 'fq2')

  # Exit if outfiles exist and no overwrite
    if os.path.exists(fq1) and os.path.getsize(fq1) > 0 and \
       os.path.exists(fq2) and os.path.getsize(fq2) > 0 and \
       not args.overwrite:
        mylogger.info('Skipping BAM to FASTQ conversion')
        return fq1, fq2
    global _didsomecomputation
    _didsomecomputation = True
    mylogger.info('Converting input BAM to FASTQ')
  # Variables
    samtools = P.option['prog']['samtools']
    inbam = args.bampath
    tmpbamstempath = File_Path(args.outdir, args.dataid, 'tmpbamstempath')
    java = P.option['prog']['java']
    samtofastq = P.option['prog']['samtofastq']
    intermediateL = [tmpbamstempath+'.bam']
  # Sort BAM reads by name
    cmd = '{samtools} sort -n {inbam} {tmpbamstempath}'.format(
        samtools=samtools, inbam=inbam, tmpbamstempath=tmpbamstempath)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to name-sort BAM ({0})'.format(inbam))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
  # Convert BAM to FASTQ, losing any readgroups associated with each BAM record
    cmd = '{java} -jar {samtofastq} INPUT={bam} FASTQ={fq1} SECOND_END_FASTQ={fq2}'.format(
        java=java, samtofastq=samtofastq, bam=inbam, fq1=fq1, fq2=fq2)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to convert BAM to FASTQ ({0}, {1})'.format(
            inbam, fq1.replace('_1.fastq', '[1|2].fastq')))
        mylogger.error('Command: {0}'.format(cmd))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
  # Remove intermediate files
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Return FASTQ paths
    return fq1, fq2

def Map_Reads_To_Target_References(args, P, mylogger, myhandler, _processname, fq1, fq2):
    'Map the FASTQ reads to the single FASTA file of references. Exit program on error.'
  # Output files
    mapbam = File_Path(args.outdir, args.dataid, 'mapbam')
    mapbamstem = File_Path(args.outdir, args.dataid, 'mapbamstem')
  # Exit if outfiles exist and no overwrite
    if os.path.exists(mapbam) and os.path.getsize(mapbam) > 0 and \
       not args.overwrite:
        mylogger.info('Skipping read mapping to reference')
        return mapbam
    global _didsomecomputation
    _didsomecomputation = True
    mylogger.info('Mapping reads to references of known genotype classification')
  # Variables
    samtools = P.option['prog']['samtools']
    mapsam = File_Path(args.outdir, args.dataid, 'mapsam')
    bwa = P.option['prog']['bwa']
    refpath = args.targetrefpath.replace('.fasta', '.fa')
    intermediateL = [mapsam]
  # Run "bwa mem" keeping shorter, non-primary split hits
    cmd = '{bwa} mem -M {refpath} {fq1} {fq2}'.format(
        bwa=bwa, refpath=refpath, fq1=fq1, fq2=fq2)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to map BAM to target refs ({0}, {1})'.format(
            fq1.replace('_1.fastq', '[1|2].fastq'), refpath))
        mylogger.error('{0}'.format(re))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
    with open(mapsam, 'w') as sam_fp:
        sam_fp.write('{0}\n'.format(ro.rstrip('\n')))
  # Convert SAM output to BAM, excluding (-F) non-primary (256) and supplementary (2048) alignments
    cmd = '{samtools} view -Sb -F 2304 {insam} | {samtools} sort -n - {outbamstem}'.format(
        samtools=samtools, insam=mapsam, outbamstem=mapbamstem)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to convert SAM to BAM ({0}, {1})'.format(mapsam, mapbam))
        mylogger.error('Command: {0}'.format(cmd))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
  # Remove intermediate files
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Return path of mapped BAM
    return mapbam

def Classify_Read_Pairs_To_Target_Genotypes(args, P, mylogger, myhandler, _processname, mapbam):
    '''
    Read the pairs of name-sorted reads.
    Genotype classification is from the mapped mate with the highest MAPQ.
    If neither mate is mapped the read pair is "unclassified".
    '''
  # Output files
    readclasspath = File_Path(args.outdir, args.dataid, 'readclasspath')
    readclassbam = File_Path(args.outdir, args.dataid, 'readclassbam')
    statspath = File_Path(args.outdir, args.dataid, 'statspath')
  # Exit if outfiles exist and no overwrite
    if os.path.exists(readclasspath) and os.path.getsize(readclasspath) > 0 and \
       os.path.exists(readclassbam) and os.path.getsize(readclassbam) > 0 and \
       os.path.exists(statspath) and os.path.getsize(statspath) > 0 and \
       not args.overwrite:
        mylogger.info('Skipping read pair classification')
        return readclasspath, readclassbam, statspath
    global _didsomecomputation
    _didsomecomputation = True
    mylogger.info('Classifying each read pair')
  # Variables
    samtools = P.option['prog']['samtools']
    samsubsetpath = File_Path(args.outdir, args.dataid, 'samsubsetpath')
    readclasssam = File_Path(args.outdir, args.dataid, 'readclasssam')
    intermediateL = [samsubsetpath, readclasssam]
  # Create intermediate file with just readid, refname, mapq
    cmd = '{samtools} view {inbam} | cut -f1,2,3,5'.format(
        samtools=samtools, inbam=mapbam)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to extract QNAME, FLAG, RNAME and MAPQ from BAM ({0})'.format(mapbam))
        mylogger.error('Command: {0}'.format(cmd))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
    with open(samsubsetpath, 'w') as out_fp:
        out_fp.write('{0}\n'.format(ro.rstrip('\n')))
  # Read each pair of records (in name-sorted order) and output allocation.
  # During this pass, also collate the number of read pairs in each readclass.
    classlistcnt = {}
    with open(readclasspath, 'w') as out_fp:
        R = ['readid', 'readclass']
        out_fp.write('{0}\n'.format('\t'.join(R)))
        try:
            with open(samsubsetpath) as in_fp:
                for line1 in in_fp:
                    line2 = in_fp.next()
                    r1_readid, r1_flag, r1_rname, r1_mapq = line1.split('\t')
                    r2_readid, r2_flag, r2_rname, r2_mapq = line2.split('\t')
                    r1_ismapped = not int(r1_flag)&4
                    r2_ismapped = not int(r2_flag)&4
                    if not r1_ismapped and not r2_ismapped:
                        readclass = "unclassified"
                    elif not r1_ismapped and r2_ismapped:
                        readclass = r2_rname.split('_')[0]
                    elif r1_ismapped and not r2_ismapped:
                        readclass = r1_rname.split('_')[0]
                    else: # if both mapped, return gt of highest scoring mapping,
                          # if a tie and reads map to different genotypes return "ambiguous"
                        r1_mapq = int(r1_mapq)
                        r2_mapq = int(r2_mapq)
                        if r1_mapq > r2_mapq:
                            readclass = r1_rname.split('_')[0]
                        elif r1_mapq < r2_mapq:
                            readclass = r2_rname.split('_')[0]
                        else:
                            r1_readclass = r1_rname.split('_')[0]
                            r2_readclass = r2_rname.split('_')[0]
                            if r1_readclass == r2_readclass:
                                readclass = r1_readclass
                            else:
                                readclass = "ambiguous"
                    R = [r1_readid, readclass]
                    out_fp.write('{0}\n'.format('\t'.join(R)))
                    if readclass not in classlistcnt:
                        classlistcnt[readclass] = 0
                    classlistcnt[readclass] += 1
        except StopIteration:
            pass
  # Create .OUTSUFFIX.stats.txt file
  # Print allocations with >args.mintargetpoppct% (typically 1.0) of reads in decr freq order
    rpscnt = sum([classlistcnt[key] for key in classlistcnt.keys()])
    minrpscnt = int(args.mintargetpoppct / 100.0 * rpscnt)
    keyL_decr = sorted(classlistcnt, key=classlistcnt.get, reverse=True)
    with open(statspath, 'w') as out_fp:
        R = ['dataid', 'popN', 'readclass', 'readcnt', 'readpct', 'bampath']
        out_fp.write('{0}\n'.format('\t'.join(R)))
        otherclasscnt = 0
        N = 0
      # Readclasses with enough reads to be considered a class of its own
        for readclass in keyL_decr:
            if readclass not in ['ambiguous', 'unclassified']:
                if classlistcnt[readclass] >= minrpscnt:
                    N += 1
                    readcnt = classlistcnt[readclass]*2
                    readpct = round(classlistcnt[readclass] / float(rpscnt) * 100.0 if rpscnt else 0.0, 1)
                    bampath = File_Path(args.outdir, args.dataid, 'popNbam', N)
                    R = [args.dataid, N, readclass, readcnt, readpct, bampath]
                    out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
                else:
                    otherclasscnt += classlistcnt[readclass]
      # readclass "othertargetclasses"
        N += 1
        readcnt = otherclasscnt*2
        readpct = round(otherclasscnt / float(rpscnt) * 100.0, 1) \
            if otherclasscnt and rpscnt else 0.0
        bampath = File_Path(args.outdir, args.dataid, 'popNbam', N)
        R = [args.dataid, N, 'othertargetclasses', readcnt, readpct, bampath]
        out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
      # readclass "ambiguous"
        N += 1
        readcnt = classlistcnt['ambiguous']*2 if 'ambiguous' in classlistcnt else 0
        readpct = round(classlistcnt['ambiguous'] / float(rpscnt) * 100.0, 1) \
            if 'ambiguous' in classlistcnt and rpscnt else 0.0
        bampath = File_Path(args.outdir, args.dataid, 'popNbam', N)
        R = [args.dataid, N, 'ambiguous', readcnt, readpct, bampath]
        out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
      # readclass "unclassified"
        N += 1
        readcnt = classlistcnt['unclassified']*2 if 'unclassified' in classlistcnt else 0
        readpct = round(classlistcnt['unclassified'] / float(rpscnt) * 100.0, 1) \
            if 'unclassified' in classlistcnt and rpscnt else 0.0
        bampath = File_Path(args.outdir, args.dataid, 'popNbam', N)
        R = [args.dataid, N, 'unclassified', readcnt, readpct, bampath]
        out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
  # Create a new file .OUTSUFFIX.readclass.bam with all the mapped reads from
  # .OUTSUFFIX.mapped.bam annotated with the inferred readclass in the readgroup.
  # Retrieve existing bam header
    cmd = '{samtools} view -H {inbam}'.format(samtools=samtools, inbam=mapbam) 
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to retrieve header from BAM ({0})'.format(mapbam))
        mylogger.error('Command: {0}'.format(cmd))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
    mapbamheader = ro.rstrip('\n')
  # Retrieve list of readclasses
    cmd = 'tail -n +2 {readclasspath} | cut -f2 | sort -u'.format(readclasspath=readclasspath)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to retrieve unique list of readclasses ({0})'.format(readclasspath))
        mylogger.error('Command: {0}'.format(cmd))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
    readclassL = ro.strip('\n').split('\n') if len(ro.strip('\n')) else []
    readclassL.sort()
  # Retrieve entire file of readclass allocations for each readid
    with open(readclasspath, 'r') as in_fp:
        readclassD = dict(x.rstrip().split(None, 1) for x in in_fp)
  # Print all this info to a SAM file
    with open(readclasssam, 'w') as out_fp:
      # Append existing header from mapped bam containing @HD and @SQ records
        out_fp.write('{0}\n'.format(mapbamheader))
      # Add one @RG line for each readclass in the .OUTSUFFIX.readclass.txt file
        for readclass in readclassL:
            out_fp.write('@RG\tID:{rgid}\n'.format(rgid=readclass))
      # Iterate through records of the mapped BAM file and output record appended with RG:Z:rgid
        cmd = '{samtools} view {inbam}'.format(samtools=samtools, inbam=mapbam)
        mylogger.debug('Command: {0}'.format(cmd))
        rv, ro, re = P.sys_exec(cmd)
        if rv != 0:
            mylogger.error('Failed to retrieve contents of BAM ({0})'.format(mapbam))
            mylogger.error('Command: {0}'.format(cmd))
            Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
            sys.exit(P.err_code('ErrorSysCall'))
        recordL = ro.rstrip('\n').split('\n') if len(ro.rstrip('\n')) else []
        for record in recordL:
            R = record.split('\t')
            readid = R[0]
            R.append('RG:Z:{rgid}'.format(rgid=readclassD[readid]))
            out_fp.write('{0}\n'.format('\t'.join(R)))
  # Convert the SAM file to a (name-sorted) BAM file
    cmd = '{samtools} view -Sb -o {outbam} {insam} '.format(
        samtools=samtools, outbam=readclassbam, insam=readclasssam)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to convert SAM to BAM ({0}, {1})'.format(readclasssam, readclassbam))
        mylogger.error('Command: {0}'.format(cmd))
        Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
        sys.exit(P.err_code('ErrorSysCall'))
  # Remove intermediate files
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Return path of file containing classification of each read pair id
    return readclasspath, readclassbam, statspath

def Statspath_Info(statspath):
    'Return contents of stats file plus sampath and placeholder for samfp (default None).'
    readclassL = []
    readclassD = {}
    with open(statspath, 'r') as in_fp:
        linecnt = 0
        for line in in_fp:
            linecnt += 1
            if linecnt == 1 or not len(line.strip()):
                continue
            dataid, N, readclass, readcnt, readpct, bampath = line.rstrip('\n').split('\t')
            readclassL.append(readclass)
            readclassD[readclass] = {
                'N': int(N),
                'readcnt': int(readcnt),
                'readpct': float(readpct),
                'bampath': bampath,
                'sampath': bampath.replace('.bam', '.sam'),
                'samfp': None,
            }
    return readclassL, readclassD

def Split_Classified_Bam_Into_Population_Bams(args, P, mylogger, myhandler, _processname, readclassbam, statspath):
    'Create one BAM file for each population documented in the stats file.'
  # Read in the stats file
    readclassL, readclassD = Statspath_Info(statspath)
  # Exit if all outfiles exist and no overwrite
    allexist = True
    for i in range(0, len(readclassL)):
        readclass = readclassL[i]
        allexist = allexist & \
           (os.path.exists(readclassD[readclass]['bampath']) and \
            os.path.getsize(readclassD[readclass]['bampath']) > 0)
        pass
    if allexist and not args.overwrite:
        mylogger.info('Skipping splitting BAM into population BAMs')
        return True
    global _didsomecomputation
    _didsomecomputation = True
    mylogger.info('Splitting single BAM into population BAMs')
  # Variables
    samtools = P.option['prog']['samtools']
    intermediateL = []

  # Open the file pointers to the SAM file for each
    for readclass in readclassL:
        try:
            readclassD[readclass]['samfp'] = open(readclassD[readclass]['sampath'], 'w')
        except:
            mylogger.error('Failed to open SAM file for writing ({0})'.format(readclassD[readclass]['sampath']))
            Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
            sys.exit(P.err_code('ErrorFileOpen'))

  # Write header and records to each .OUTSUFFIX.popN.sam file
    for i, readclass in enumerate(readclassL):

      # Print header to SAM file
        cmd = '{samtools} view -H {inbam} | grep -v "^@RG"'.format(samtools=samtools, inbam=readclassbam)
        mylogger.debug('Command: {0}'.format(cmd))
        rv, ro, re = P.sys_exec(cmd)
        if rv != 0:
            mylogger.error('Failed to retrieve @HD, @SQ and @PG records of header from BAM ({0})'.format(readclassbam))
            mylogger.error('Command: {0}'.format(cmd))
            Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
            sys.exit(P.err_code('ErrorSysCall'))
        readclassD[readclass]['samfp'].write('{0}\n'.format(ro.rstrip('\n')))
        if readclass != 'othertargetclasses':
            cmd = '{samtools} view -H {inbam} | grep "^@RG" | grep "{readclass}"'.format(
                samtools=samtools, inbam=readclassbam, readclass=readclass)
        else:
            cmd = '{samtools} view -H {inbam} | grep "^@RG" | grep -v "ambiguous" | grep -v "unclassified"'.format(
                samtools=samtools, inbam=readclassbam)
            for rc in readclassL:
                if rc not in ['ambiguous', 'unclassified', 'othertargetclasses']:
                    cmd += ' | grep -v "{readclass}"'.format(readclass=rc)
        mylogger.debug('Command: {0}'.format(cmd))
        rv, ro, re = P.sys_exec(cmd)
        if len(ro):
            readclassD[readclass]['samfp'].write('{0}\n'.format(ro.rstrip('\n')))

      # Print reads to SAM file
        if readclass != 'othertargetclasses':
            cmd = '{samtools} view {inbam} | grep "RG:Z:{readclass}"'.format(
                samtools=samtools, inbam=readclassbam, readclass=readclass)
        else:
            cmd = '{samtools} view {inbam} | grep -v "RG:Z:ambiguous" | grep -v "RG:Z:unclassified"'.format(
                samtools=samtools, inbam=readclassbam)
            for rc in readclassL:
                if rc not in ['ambiguous', 'unclassified', 'othertargetclasses']:
                    cmd += ' | grep -v "RG:Z:{readclass}"'.format(readclass=rc)
        mylogger.debug('Command: {0}'.format(cmd))
        rv, ro, re = P.sys_exec(cmd)
        if len(ro):
            readclassD[readclass]['samfp'].write('{0}\n'.format(ro.rstrip('\n')))

  # Close the file pointers to the SAM file for each
    for readclass in readclassL:
        try:
            readclassD[readclass]['samfp'].close()
        except:
            mylogger.error('Failed to close SAM file ({0})'.format(readclassD[readclass]['sampath']))
            Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
            sys.exit(P.err_code('ErrorFileClose'))

  # Convert each .OUTSUFFIX.popN.sam file to BAM
    for readclass in readclassL:
        outbamstem=readclassD[readclass]['bampath'].replace('.bam', '')
        cmd = '{samtools} view -Sb {insam} | {samtools} sort -n - {outbamstem}'.format(
            samtools=samtools, insam=readclassD[readclass]['sampath'], outbamstem=outbamstem)
        mylogger.debug('Command: {0}'.format(cmd))
        rv, ro, re = P.sys_exec(cmd)
        if rv != 0:
            mylogger.error('Failed to convert SAM to BAM ({0}, {1})'.format(mapsam, mapbam))
            mylogger.error('Command: {0}'.format(cmd))
            Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, False)
            sys.exit(P.err_code('ErrorSysCall'))

  # Remove intermediate files
    if args.deleteints:
        for readclass in readclassL:
            intermediateL.append(readclassD[readclass]['sampath'])
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Return path of file containing key to the .OUTSUFFIX.popN.bam files
    return True

def Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, finishedok=True):
    'Delete all intermediate (i.e., temporary) files if deleteints is True.'
    global _didsomecomputation
  # Get list of all possible intermediate files
    readclassL, readclassD = Statspath_Info(File_Path(args.outdir, args.dataid, 'statspath'))
    intermediateL = [File_Path(args.outdir, args.dataid, filekey) \
        for filekey in _outfilesuffixD.keys() \
        if _outfilesuffixD[filekey]['istmp']]
    for readclass in readclassL:
        intermediateL.append(readclassD[readclass]['sampath'])
  # Print start message
    mylogger.info('Tidying up')
  # Delete any remaining intermediate files
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
                    _didsomecomputation = True
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Create the sentinel path if it doesn't already exist XXXX Need to think about when might need to create this
    sentinel_path = File_Path(args.outdir, args.dataid, 'sentinelpath')
    if finishedok:
        if not os.path.exists(sentinel_path) or args.overwrite or _didsomecomputation:
            now = time.localtime()
            with open(sentinel_path, 'w') as out_fp:
                out_fp.write('{0}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', now)))
    else:
        try:
            if os.path.exists(sentinel_path):
                os.remove(sentinel_path)
        except:
            mylogger.info('Failed to remove sentinel file during tidying ({0})'.format(sentinel_path))
    return 0

def Process(args, P, mylogger, myhandler, processname, suffixL):
    'Execute main processing steps of this subtool.'
    mylogger.info('Processing')
  # Exit if overwrite is False and all output files already exist
    if not args.overwrite:
        alreadydone = True
        for suffix in suffixL:
            outpath = os.path.join(args.outdir, args.dataid+suffix)
            if not os.path.exists(outpath) or os.path.getsize(outpath) == 0:
                alreadydone = False
        statspath = File_Path(args.outdir, args.dataid, 'statspath')
        if alreadydone and os.path.exists(statspath):
            readclassL, readclassD = Statspath_Info(statspath)
            for readclass in readclassL:
                bampath = readclassD[readclass]['bampath']
                if not os.path.exists(bampath) or os.path.getsize(bampath) == 0:
                    alreadydone = False
        if alreadydone:
            mylogger.info('Skipping processing as all output files already exist and are non-empty')
            if args.deleteints:
                fq1 = File_Path(args.outdir, args.dataid, 'fq1')
                fq2 = File_Path(args.outdir, args.dataid, 'fq2')
                mapbam = File_Path(args.outdir, args.dataid, 'mapbam')
                readclasspath = File_Path(args.outdir, args.dataid, 'readclasspath')
                Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, True)
            return 0
  # Do the processing
    fq1, fq2 = Convert_Input_Reads_BAM_To_FASTQ(args, P, mylogger, myhandler, _processname)
    mapbam = Map_Reads_To_Target_References(args, P, mylogger, myhandler, _processname, fq1, fq2)
    readclasspath, readclassbam, statspath = Classify_Read_Pairs_To_Target_Genotypes(
        args, P, mylogger, myhandler, _processname, mapbam)
    splitok = Split_Classified_Bam_Into_Population_Bams(args, P, mylogger, myhandler, _processname, readclassbam, statspath)
    Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, True)
    mylogger.debug('Processing finished')
    return 0

def main(parser, args, P, mylogger, myhandler, argv):
    'Execute the code for this subtool.'
    mylogger.info('Started')
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    finalsuffixL = [_outfilesuffixD[filekey]['suffix'] for filekey in _outfilesuffixD.keys() if _outfilesuffixD[filekey]['isfinal']]
    Process(args, P, mylogger, myhandler, _processname, finalsuffixL)
    mylogger.info('Finished')
    if args.deleteints:
        filekeyL = [filekey for filekey in _outfilesuffixD.keys() if not _outfilesuffixD[filekey]['istmp']]
    else:
        filekeyL = _outfilesuffixD.keys()
    filekeyL.sort()
    for filekey in filekeyL:
        msg = 'Output: {filepath} : {description}'.format(
            filepath=File_Path('', args.dataid, filekey),
            description=_outfilesuffixD[filekey]['msg'])
        if _outfilesuffixD[filekey]['istmp']:
            msg += ' [tmp]'
        mylogger.info(msg)
    return 0

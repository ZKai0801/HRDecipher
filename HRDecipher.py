__doc__ = """ calculate HRD genomic scar score, as well as visualise scar regions """
__author__ = "Zhentian Kai"
__date__ = "2020/10/26"
__version__ = "v1.0"

import pandas as pd
import os
import argparse


CENTRO = {'chr1': [121500000, 128900000], 'chr2': [90500000, 96800000], 
          'chr3': [87900000, 93900000], 'chr4': [48200000, 52700000], 
          'chr5': [46100000, 50700000], 'chr6': [58700000, 63300000], 
          'chr7': [58000000, 61700000], 'chr8': [43100000, 48100000], 
          'chr9': [47300000, 50700000], 'chr10': [38000000, 42300000], 
          'chr11': [51600000, 55700000], 'chr12': [33300000, 38200000], 
          'chr13': [16300000, 19500000], 'chr14': [16100000, 19100000], 
          'chr15': [15800000, 20700000], 'chr16': [34600000, 38600000], 
          'chr17': [22200000, 25800000], 'chr18': [15400000, 19000000], 
          'chr19': [24400000, 28600000], 'chr20': [25600000, 29400000], 
          'chr21': [10900000, 14300000], 'chr22': [12200000, 17900000], 
          'chrX': [58100000, 63000000], 'chrY': [11600000, 13400000]}


def main(fname, methods, ofname):
    df = pd.read_csv(fname, sep="\t")
    segments = preprocess(df)
    seg = segments.copy()

    loh = calc_loh(seg)
    loh['HRD_tag'] = 'LOH'

    seg = segments.copy()
    tai = calc_tai(seg)

    seg = segments.copy()
    lst = calc_lst(seg)
    
    hrd = pd.DataFrame()
    hrd = hrd.append(loh, sort=False)
    hrd = hrd.append(tai, sort=False)
    hrd = hrd.append(lst, sort=False)
    hrd.reset_index(drop=True, inplace=True)
    hrd.to_csv(ofname, sep="\t", index=None)


def preprocess(df):
    # check integrity of input file
    for index, row in df.iterrows():
        if row['total_cn'] != row['A_cn'] + row['B_cn']:
            raise IOError("The total_cn is not equal to sum of A_cn and B_cn. Please check your input file")
        if row['Start_position'] > row['End_position']:
            raise IOError("Start_position is larger than End_position. Please check your input file")
    
    # check if the input panel covers all chromosomes
    no_chroms = [i for i in ['chr'+str(i) for i in range(1,23)] if i not in df['Chromosome'].tolist()]
    if no_chroms:
        no_chroms = str(no_chroms)
        raise IOError(f"Input dataframe does not contain {no_chroms}. Please make sure your assay cover all human chromsome")

    # remove NA
    df = df.dropna()

    # remove chromosome X/Y
    df = df[(df['Chromosome'] != 'chrX') & (df['Chromosome'] != 'chrY')]

    # keep A_cn column >= B_cn column
    for index, row in df.iterrows():
        if row['A_cn'] < row['B_cn']:
            df.loc[index, 'A_cn'], df.loc[index, 'B_cn'] = int(row['B_cn']), int(row['A_cn'])
    
    df = combine_segments(df)
    return df


def combine_segments(segments):
    """
    combine nearby segments if possible
    i.e. If A_cn & B_cn is identical between nearby segments,
    then merge those segments into one, with End_position changed to 
    the End_position value from the second segment. 
    total_cn is also changed to the sum of those merged segments
    """
    new_segments = pd.DataFrame(columns=segments.columns)
    for index, row in segments.iterrows():
        if index == 0:
            old_row = row
            continue
        if ((old_row["Chromosome"] == row["Chromosome"]) and 
            (old_row["A_cn"] == row["A_cn"]) and 
            (old_row["B_cn"] == row["B_cn"])):
            old_row['End_position'] = row['End_position']
        else:
            # if dataframe reach the last row, add the last
            # row if it cannot be combined
            if index != segments.shape[0] - 1:
                new_segments = new_segments.append(old_row)
            else:
                new_segments = new_segments.append(old_row)
                new_segments = new_segments.append(row)
            old_row = row
    
    new_segments.reset_index(drop=True, inplace=True)
    return new_segments


def calc_loh(segments, size_limit=15000000):
    """
    HRD-LOH was defined as LOH regions that exceed 15Mb and 
    does not cover the whole chromosome;  
    """
    # combine all non-deletion segments
    segments['A_cn'].values[segments['A_cn'].values >= 1] = 1
    segments = combine_segments(segments)

    # find LOH that cover the whole chromosome
    # Note: This logic is only true when the assay does cover whole chromosomes
    # May need further modification 
    chrom_loh = []
    all_chroms = set(segments['Chromosome'])
    for chrom in all_chroms:
        chrom_seg = segments[segments['Chromosome'] == chrom]
        if (sum(chrom_seg['B_cn']) == 0) and (sum(chrom_seg['A_cn']) != 0):
            chrom_loh.append(chrom)
    segments = segments[~segments['Chromosome'].isin(chrom_loh)]

    # get HRD-LOH
    all_loh = segments[(segments['B_cn'] == 0) & (segments['A_cn'] != 0)]
    hrd_loh = all_loh[all_loh['End_position'] - all_loh['Start_position'] > size_limit]
    return hrd_loh


def determine_lst(segments, size_limit=3000000, segment_cutoff=10000000):
    """
    Large scale transition (LST) was defined as chromosomal break between adjacent regions of at least 10 Mb, 
    with a distance between them not larger than 3Mb.
    """
    lst = pd.DataFrame(columns=['Chromosome', 'Start_position', 'End_position', 'HRD_tag'])
    i = 0
    segments.reset_index(drop=True, inplace=True)
    for index, row in segments.iterrows():
        if index == 0:
            old_row = row
            continue
        if ((old_row['End_position'] - old_row['Start_position'] >= segment_cutoff) and 
            (row['End_position'] - row['Start_position'] >= segment_cutoff) and 
            (row['Start_position'] - old_row['End_position'] <= size_limit)):
            lst.loc[i, 'Chromosome'] = row['Chromosome']
            lst.loc[i, 'Start_position'] = old_row['End_position']
            lst.loc[i, 'End_position'] = row['Start_position']
            lst.loc[i, 'HRD_tag'] = 'LST'
            i += 1
        old_row = row
    return lst


def calc_lst(segments):
    """
    Chromosomes need to split into arms;
    and remove short segments in advance
    """
    hrd_lst = pd.DataFrame(columns=['Chromosome', 'Start_position', 'End_position', 'HRD_tag'])
    
    for chrom in ['chr'+str(i) for i in range(1,23)]:
        chrom_seg = segments[segments['Chromosome'] == chrom]
        
        # skip chromosome with only one segments
        if chrom_seg.shape[0] < 2:
            continue

        # split into chromosome arms
        p_arm = chrom_seg[chrom_seg['Start_position'] <= CENTRO[chrom][0]]
        q_arm = chrom_seg[chrom_seg['End_position'] >= CENTRO[chrom][1]]
        p_arm.reset_index(drop=True, inplace=True)
        q_arm.reset_index(drop=True, inplace=True)

        # reset segments' coordinates
        if p_arm.shape[0] > 0:
            p_arm = combine_segments(p_arm)
            p_arm.loc[p_arm.shape[0]-1, "End_position"] = CENTRO[chrom][0]
        
        if q_arm.shape[0] > 0:
            q_arm = combine_segments(q_arm)
            q_arm.loc[0, "Start_position"] = CENTRO[chrom][1]

        # remove short segments <= 3,000,000bp
        no_3mb = p_arm[p_arm['End_position'] - p_arm['Start_position'] < 3000000].shape[0]
        while no_3mb > 0:
            p_arm = p_arm[p_arm['End_position'] - p_arm['Start_position'] >= 3000000]
            p_arm = combine_segments(p_arm)
            no_3mb = p_arm[p_arm['End_position'] - p_arm['Start_position'] < 3000000].shape[0]

        no_3mb = q_arm[q_arm['End_position'] - q_arm['Start_position'] < 3000000].shape[0]
        while no_3mb > 0:
            q_arm = q_arm[q_arm['End_position'] - q_arm['Start_position'] >= 3000000]
            q_arm = combine_segments(q_arm)
            no_3mb = q_arm[q_arm['End_position'] - q_arm['Start_position'] < 3000000].shape[0]
        
        # determine LST
        if p_arm.shape[0] >= 2 :
            lst = determine_lst(p_arm)
            hrd_lst = hrd_lst.append(lst, sort=False)
        if q_arm.shape[0] >= 2:
            lst = determine_lst(q_arm)
            hrd_lst = hrd_lst.append(lst, sort=False)
          
    return hrd_lst


def calc_tai(segments, size_limit=1000000, ploidy_by_chrom=True):
    """
    Telomeric Allelic Imbalance (TAI) was defined as regions of unequal distribution of
    parental allele sequence with or without changes in the overall copy number that 
    extend to the telomeric end of a chromosome.
    However, the size of limit of TAI is not defined in the original article:
    DOI: 10.1158/2159-8290.CD-11-0206
    Here, I applied the same limit used in scarHRD (Birkbak et al., 2012)
    
    :param ploidy_by_chrom -- option to determine ploidy per chromosome
    """
    segments = segments[segments['End_position'] - segments['Start_position'] > size_limit]
    segments.reset_index(drop=True, inplace=True)
    segments = combine_segments(segments)
    
    hrd_tai = pd.DataFrame()
    for chrom in ['chr'+str(i) for i in range(1,23)]:
        chrom_seg = segments[segments['Chromosome'] == chrom]

        # skip chromosome with only one segments
        if chrom_seg.shape[0] < 2:
            continue
        
        if ploidy_by_chrom:
            # find A_cn for the longest segments
            tmp = (chrom_seg['End_position'] - chrom_seg['Start_position']).astype(int)
            ploidy = chrom_seg.loc[tmp.idxmax(), 'A_cn']

        chrom_seg.loc[:,'ploidy'] = ploidy

        # determine allelic imbalance
        if (ploidy == 1) or (ploidy % 2 == 0):
            chrom_seg.loc[chrom_seg['A_cn'] != chrom_seg['B_cn'], 'HRD_tag'] = 'AI'
        if (ploidy != 1) and (ploidy %2 == 1):
            chrom_seg.loc[chrom_seg['A_cn'] + chrom_seg['B_cn'] != ploidy, 'HRD_tag'] = 'AI'
            chrom_seg.loc[chrom_seg['B_cn'] == 0, 'HRD_tag'] = 'AI'

        # check if short arms have TAI
        chrom_seg.reset_index(drop=True, inplace=True)
        if ((chrom_seg.loc[0, 'HRD_tag'] == 'AI') and 
            (chrom_seg.shape[0] != 1 ) and 
            (chrom_seg.loc[0, 'End_position'] < CENTRO[chrom][0])):
           chrom_seg.loc[0, 'HRD_tag'] = 'TAI'

        # check if long arms have TAI
        last_index = chrom_seg.shape[0] - 1
        if ((chrom_seg.loc[last_index, 'HRD_tag'] == 'AI') and 
            (chrom_seg.shape[0] != 1 ) and 
            (chrom_seg.loc[last_index, 'Start_position'] > CENTRO[chrom][1])):
           chrom_seg.loc[last_index, 'HRD_tag'] = 'TAI'
        
        hrd_tai = hrd_tai.append(chrom_seg[chrom_seg['HRD_tag'] == 'TAI'])
    
    return hrd_tai
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", 
                         help="""sampleID.pre_hrd.tsv, must contain following columns: 
                                 Chromosome, Start_position, End_position, total_cn, 
                                 A_cn, B_cn, ploidy""")
    parser.add_argument("-m", "--methods", type=str, default='sum', choices=['sum', 'myriad'], 
                        help="""Two methods are commonly used in HRD score calculation,
                        sum is the most widely used one, which simply add three score together;
                        myriad is the method described in Timms et al., 2014. 
                        DOI: 10.1186/s13058-014-0475-x""")
    parser.add_argument("-o", "--output", help="Output filename")
    args = parser.parse_args()

    wrkdir = os.path.dirname(os.path.abspath(args.input))
    sampleID = os.path.basename(args.input).split(".")[0]
    if not args.output:
        ofname = os.path.join(wrkdir, sampleID+".hrd.tsv")
    else:
        ofname = args.output    

    main(args.input, args.methods, ofname)


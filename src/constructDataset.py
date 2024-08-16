import sys
from typing import List, Dict, Tuple
import numpy as np
from dataclasses import dataclass
import csv
from tqdm import tqdm

@dataclass
class segment:
    indiv1: str
    hap1: int
    indiv2: str
    hap2: int
    start: int
    end: int


# split segment into chunks

def getIntervals(s: int, e: int, interval_size: int = 10) -> List[Tuple[int]]:
    """_summary_

    Args:
        s (int): start point for segment
        e (int): end point for segment
        interval_size (int): desired chunk size (will cause more features in final dataset)

    Returns:
        List[Tuple[int]]: _description_
    """
    return [(s + i * (e - s) // interval_size, s + (i + 1) * (e - s) // interval_size) for i in range(interval_size)]

# calculate mismatch rate for chunk
def loadVCFData(vcf_file_path: str):

    f = open(vcf_file_path, 'r')
    vcf_data = f.readlines()
    f.close()
    i = 0 # going to use this to determine where to cutoff header, can be different based on what has been done via bcf tools
    while vcf_data[i][0:6] != "#CHROM":
        i += 1
    vcf_data = vcf_data[i:]
    vcf_data = [i.strip().split() for i in vcf_data]
    ids_dict = {vcf_data[0][i]: i for i in range(9, len(vcf_data[0]))} # used to determine which column index to access based on sample id
    sites_array = [int(vcf_data[i][1]) for i in range(1, len(vcf_data))]
    return vcf_data, ids_dict, sites_array


def searchForStart(bp_start: int, sites_array: List[int]) -> int:
    """
    binary search for starting point when iterating through vcf
    """
    l, r = 0, len(sites_array) - 1
    while l < r:
        mid = (l + r) // 2
        if sites_array[mid] == bp_start:
            return mid
        elif sites_array[mid] > bp_start:
            r = mid - 1
        else:
            l = mid + 1
    return l

def processSegment(input_string: str) -> segment:
    """
    create segment object from string
    """
    lst = input_string.strip().split()
    hap_id1 = int(lst[0])
    hap_id2 = int(lst[1])
    start = int(float(lst[2]))
    end = int(float(lst[3]))
    haplotype1 = hap_id1 % 2
    haplotype2 = hap_id2 % 2
    tsk_id1 = "tsk_" + str(hap_id1 // 2)
    tsk_id2 = "tsk_" + str(hap_id2 // 2)
    return segment(tsk_id1, haplotype1, tsk_id2, haplotype2, start, end)

def findMismatches(segment, sites_array: List[int], vcf_data: List[List], ids_dict: Dict, start = None, end = None):
    if start == None:
        start = segment.start
    if end == None:
        end = segment.end
    mismatches = 0 # counter for number of mismatches
    m = searchForStart(start, sites_array) + 1 # starting index for iterating
    m_init = m # save starting index
    idx1 = ids_dict[segment.indiv1]
    idx2 = ids_dict[segment.indiv2]
    while m < len(vcf_data) - 1 and int(vcf_data[m][1]) <= end:
        allele1 = vcf_data[m][idx1].split("|")[segment.hap1]
        allele2 = vcf_data[m][idx2].split("|")[segment.hap2]
        """
        looks complicated, but we're taking the ith entry in the vcf list (one row of vcf file/one site),
        then taking the (idx)th entry in that list which will be a pair of alleles
        those alleles are a string in the format "x|x", so we split on that, and the corresponding allele from the haplotype
        given from hap-ibd is the allele we're looking for
        """
        if allele1 != allele2:
            mismatches += 1 #increment number of mismatches
        m += 1
        # print(i)
    return float(mismatches) / float(m - m_init) # divide number of mismatches found by the number of sites visited

def findCorrections(segment, vcf:List[List], vcf_smooth: List[List], ids_dict: dict, sites_array: List[int], start=None, end = None):
    if start == None:
        start = segment.start
    if end == None:
        end = segment.end
    vcf_idx1 = ids_dict[segment.indiv1]
    vcf_idx2 = ids_dict[segment.indiv2]
    m = searchForStart(start, sites_array) + 1
    m_init = m
    n_corrections = 0
    while m < len(vcf) - 1 and int(vcf[m][1]) <= end:
        vcf_val1 = vcf[m][vcf_idx1].split("|")[segment.hap1]
        vcf_val2 = vcf[m][vcf_idx2].split("|")[segment.hap2]
        vcf_val_s1 = vcf_smooth[m][vcf_idx1].split("|")[segment.hap1]
        vcf_val_s2 = vcf_smooth[m][vcf_idx2].split("|")[segment.hap2]
        if vcf_val1 != vcf_val_s1:
            n_corrections += 1
        if vcf_val2 != vcf_val_s2:
            n_corrections += 1
        m += 1
    return n_corrections / (m - m_init)
    
    


# write in format: hap1 hap2 start end mmr1 ... mmr10


def main():
    vcf_fp = sys.argv[1] # unsmoothed vcf
    vcf_fp_smooth = sys.argv[2] # smoothed vcf
    gt_segments_fp = sys.argv[3] # file path to gt segments
    rp_segments_fp = sys.argv[4] # file path to "true"  reported segments
    fi_segments_fp = sys.argv[5] # file path to "false" reported segments
    chunk_size = int(sys.argv[6]) # desired chunk size i.e. how many chunks you want to split the semgents into for analysis
    output_file = sys.argv[7]
    vcf, ids_dict, sites = loadVCFData(vcf_fp)
    vcf_s, _, _ = loadVCFData(vcf_fp_smooth)
    fieldnames = ["id1", "hap1", "id2", "hap2", "start", "end"]
    fieldnames += ["mmr_" + str(i) for i in range(1, chunk_size + 1)]
    fieldnames += ["cr_" + str(i) for i in range(1, chunk_size + 1)]
    fieldnames.append("classification")
    with open(output_file, "w+", newline="\n") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(fieldnames)
        # create gt_segments data
        f = open(gt_segments_fp, 'r')
        for line in tqdm(f.readlines()):
            seg = processSegment(line)
            intervals = getIntervals(seg.start, seg.end, chunk_size)
            segment_rates = []
            correction_rates = []
            for interval in intervals:
                r = findMismatches(seg, sites, vcf_s, ids_dict, interval[0], interval[1])
                segment_rates.append(r)
            for interval in intervals:
                c = findCorrections(seg, vcf, vcf_s, ids_dict, sites, interval[0], interval[1])
                correction_rates.append(c)
            out = [seg.indiv1, seg.hap1, seg.indiv2, seg.hap2, seg.start, seg.end]
            out += segment_rates
            out += correction_rates
            out.append(1)
            writer.writerow(out)
        f.close()
        f = open(rp_segments_fp, 'r')
        for line in tqdm(f.readlines()):
            seg = processSegment(line)
            intervals = getIntervals(seg.start, seg.end, chunk_size)
            segment_rates = []
            correction_rates = []
            for interval in intervals:
                r = findMismatches(seg, sites, vcf_s, ids_dict, interval[0], interval[1])
                segment_rates.append(r)
            for interval in intervals:
                c = findCorrections(seg, vcf, vcf_s, ids_dict, sites, interval[0], interval[1])
                correction_rates.append(c)
            out = [seg.indiv1, seg.hap1, seg.indiv2, seg.hap2, seg.start, seg.end]
            out += segment_rates
            out += correction_rates
            out.append(2)
            writer.writerow(out)
        f.close()
        f = open(fi_segments_fp, 'r')
        for line in tqdm(f.readlines()):
            seg = processSegment(line)
            intervals = getIntervals(seg.start, seg.end, chunk_size)
            segment_rates = []
            correction_rates = []
            for interval in intervals:
                r = findMismatches(seg, sites, vcf_s, ids_dict, interval[0], interval[1])
                segment_rates.append(r)
            for interval in intervals:
                c = findCorrections(seg, vcf, vcf_s, ids_dict, sites, interval[0], interval[1])
                correction_rates.append(c)
            out = [seg.indiv1, seg.hap1, seg.indiv2, seg.hap2, seg.start, seg.end]
            out += segment_rates
            out += correction_rates
            out.append(3)
            writer.writerow(out)
        f.close()
            


if __name__ == "__main__":
    main()

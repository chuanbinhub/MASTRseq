'''
MASTR-seq STR counting script <run_str_counter.py>
Written by Han-Seul Ryu (Last modified Sept 16th 2024)
Modified 8/4/25 by Tmalach to condense STR calling scripts
Modified by Chuanbin Su to integrate STR calling scripts into snakemake container (Aug 2025)

INPUT: 
aligned, demultiplexed, and sorted bam file processed from Dorado pipeline, 

OUTPUT: 
(1) filters.txt file that stores your filtering parameters, 
(2) counts.txt files with read names and STR count
'''
#!/usr/bin/env python
import pysam
import os
import sys
import subprocess

def main():

    strType = sys.argv[1]
    input_bam = sys.argv[2]
    output_path = sys.argv[3].rstrip('/')

    #parameters below were used for C9orf72 G2C4 counting. 
    if strType.upper() == 'C9ORF72':
        chrom = 'chr9'
        fetch_start = 27560000
        fetch_end = 27574000
        str_seq = 'GGCCCC' #Please enter in forward sequence (e.g. 'GGCCCC' for C9orf72 G4C2 repeat)
        min_fragment_len = 4000 #enter minimum length of expected Cas9 cleaved fragment
        ## any read shorter than min_fragment_len will be discarded

        ## Would you like to filter for reads that include a highly accurate sequence around the STR? 
        filter_accurate_sequence = 'N'         #Enter 'Y' or 'N'
        accurate_sequence = 'None'             #Enter sequence if 'Y'; Enter 'None' if 'N' above
        accurate_seq_where = 'None'            #Enter 'U' if the sequence is upstream of the STR and 'D' if sequence is downstream of STR
                                               #Enter 'None' if 'N' above
        
        ## Would you like to filter for reads that have a particular sequence directly upstream of the STR ? 
        filter_direct_upstr = 'Y'              #Enter 'Y' or 'N'
        direct_upstr_seq = 'GCCCC'             #Enter sequence if 'Y'; Enter 'None' if 'N' above
        
        ## Would you like to filter for reads that have a particular sequence directly downstream of the STR ?
        filter_direct_downstr = 'Y'            #Enter 'Y' or 'N'
        direct_downstr_seq = 'TAGCG'           #Enter sequence if 'Y'; Enter 'None' if 'N' above
        
        ## Would you like to filter out reads that have lots of non-STR bases within the STR?
        filter_interrupters = 'Y'              #Enter 'Y' or 'N'
        num_interrupting_bases_allowed = 100   #Enter max number of bases allowed between two [STR]x4

        ## Would you like to check if there are normal-length reads that pass all filters above and have <4 repeats?
        ## if you enter 'N', you will not retain any reads that have <4 repeats 
        check_for_less_than_four_repeats = 'Y'     #Enter 'Y' or 'N'
    
    #[HTT_CAG] parameters used for CAG counting.    
    elif strType.upper() == 'HTT':
        chrom = 'chr4'
        fetch_start = 3073000
        fetch_end = 3077000
        str_seq = 'CAG'          #Please enter in forward sequence (e.g. 'GGCCCC' for C9orf72 G4C2 repeat)
        min_fragment_len = 6000  #enter minimum length of expected Cas9 cleaved fragment
        ## any read shorter than min_fragment_len will be discarded

        ## Would you like to filter for reads that include a highly accurate sequence around the STR? 
        filter_accurate_sequence = 'N'         #Enter 'Y' or 'N'
        accurate_sequence = 'None'             #Enter sequence if 'Y'; Enter 'None' if 'N' above
        accurate_seq_where = 'None'            #Enter 'U' if the sequence is upstream of the STR and 'D' if sequence is downstream of STR
                                               #Enter 'None' if 'N' above

        ## Would you like to filter for reads that have a particular sequence directly upstream of the STR ? 
        filter_direct_upstr = 'Y'              #Enter 'Y' or 'N'
        direct_upstr_seq = 'CTTC'             #Enter sequence if 'Y'; Enter 'None' if 'N' above
    
        ## Would you like to filter for reads that have a particular sequence directly downstream of the STR ?
        filter_direct_downstr = 'N'            #Enter 'Y' or 'N'
        direct_downstr_seq = 'None'           #Enter sequence if 'Y'; Enter 'None' if 'N' above

        ## Would you like to filter out reads that have lots of non-STR bases within the STR?
        filter_interrupters = 'Y'              #Enter 'Y' or 'N'
        num_interrupting_bases_allowed = 100   #Enter max number of bases allowed between two [STR]x4

        ## Would you like to check if there are normal-length reads that pass all filters above and have <4 repeats?
        ## if you enter 'N', you will not retain any reads that have <4 repeats 
        check_for_less_than_four_repeats = 'N'     #Enter 'Y' or 'N' 

    #[FMR1_CGG] for parameters used for CGG counting.    
    elif strType.upper() == 'FMR1':
        chrom = 'chrX'
        fetch_start = 147911000
        fetch_end = 147918000
        str_seq = 'CGG' #Please enter in forward sequence (e.g. 'GGCCCC' for C9orf72 G4C2 repeat)
        min_fragment_len = 5900 #enter minimum length of expected Cas9 cleaved fragment
        ## any read shorter than min_fragment_len will be discarded

        ## Would you like to filter for reads that include a highly accurate sequence around the STR? 
        filter_accurate_sequence = 'Y'         #Enter 'Y' or 'N'
        accurate_sequence = 'ACCAAACCAA'             #Enter sequence if 'Y'; Enter 'None' if 'N' above
        accurate_seq_where = 'U'            #Enter 'U' if the sequence is upstream of the STR and 'D' if sequence is downstream of STR
                                               #Enter 'None' if 'N' above
        
        ## Would you like to filter for reads that have a particular sequence directly upstream of the STR ? 
        filter_direct_upstr = 'N'              #Enter 'Y' or 'N'
        direct_upstr_seq = 'None'             #Enter sequence if 'Y'; Enter 'None' if 'N' above
    
        ## Would you like to filter for reads that have a particular sequence directly downstream of the STR ?
        filter_direct_downstr = 'N'            #Enter 'Y' or 'N'
        direct_downstr_seq = 'None'           #Enter sequence if 'Y'; Enter 'None' if 'N' above
        
        ## Would you like to filter out reads that have lots of non-STR bases within the STR?
        filter_interrupters = 'Y'              #Enter 'Y' or 'N'
        num_interrupting_bases_allowed = 120   #Enter max number of bases allowed between two [STR]x4

        ## Would you like to check if there are normal-length reads that pass all filters above and have <4 repeats?
        ## if you enter 'N', you will not retain any reads that have <4 repeats 
        check_for_less_than_four_repeats = 'N'     #Enter 'Y' or 'N' 
        
    else:
        raise Exception('Please enter FMR1, HTT, or C9orf72!')

    ####################### End of filter parameter input. Do not modify anything after this line #########################

    if filter_accurate_sequence == 'Y':
        if accurate_seq_where not in ['U','D']:
            raise Exception("Error: Please enter either 'U' or 'D'")
    
    # make output files to store all counts and names of filtered reads 
    base_file_name = os.path.basename(input_bam).replace('.sorted.bam', '').replace('.bam', '')
    
    if not os.path.exists(output_path): 
        os.makedirs(output_path)
        
    all_revcounts = open(output_path + '/' + base_file_name+'_revcounts.txt','w')
    all_fwdcounts = open(output_path + '/' + base_file_name+'_fwdcounts.txt','w')
    all_counts = open(output_path + '/' + base_file_name+'_allcounts.txt','w')

    filter_txt = open(output_path + '/' + base_file_name+'_filters_used.txt','w')

    print('\nNow running STR counter!')
    bamfile = pysam.AlignmentFile(input_bam, "rb",check_sq=False)
    region = bamfile.fetch(chrom,fetch_start,fetch_end)

    total_reads = 0
    pass_filter = 0
    no_pass = 0

    str_seq_four = str_seq*4
    for x in region:
        total_reads += 1
        x = str(x)
        x_split = x.split('\t')
        name = x_split[0]
        fwd_rev = x_split[1]
        seq = x_split[9]

        if len(seq) < min_fragment_len:
            no_pass += 1
            continue 
        
        already_counted = False
        first = seq.find(str_seq_four) #find first instance of four consecutive repeats

        #if there's no instance of four consecutive repeats, check if the user wanted to rescue the reads with 2-3 repeats
        if first == -1: 
            if check_for_less_than_four_repeats == 'Y': 
                str_seq_two = str_seq * 2
                first_two_repeats = seq.find(str_seq_two)

                if first_two_repeats == -1:     #if there are no instances of 2 repeats, discard the read and move on to the next read
                    no_pass += 1
                    continue
                else:       # if there is an instance of 2 repeats, determine if read passes all other filters specified by user
                    if filter_accurate_sequence == 'Y':
                        test_accurate_seq = pass_accurate_seq_two(seq,str_seq,accurate_sequence,accurate_seq_where)
                        if test_accurate_seq == False:
                            no_pass += 1
                            continue
                    if filter_direct_upstr == 'Y':
                        test_direct_upstr = pass_direct_upstr_two(seq,str_seq,direct_upstr_seq)
                        if test_direct_upstr == False:
                            no_pass += 1
                            continue
                    if filter_direct_downstr == 'Y':
                        test_direct_downstr = pass_direct_downstr_two(seq,str_seq,direct_downstr_seq)
                        if test_direct_downstr == False:
                            no_pass += 1
                            continue
                    last_two_repeats = seq.rfind(str_seq_two)
                    cut_seq = seq[first_two_repeats:last_two_repeats+len(str_seq_two)]
                    already_counted = True
            
            else:  # if there's no instance of four repeats and user did not want to rescue reads with 2-3 repeats, discard read and move on to the next read        
                no_pass += 1
                continue

        # for reads with at least 4 consecutive repeats in tandem, determine if read passes all user-defined filters
        if already_counted == False: 
            last = seq.rfind(str_seq_four)
            if filter_accurate_sequence == 'Y':
                test_accurate_seq = pass_accurate_seq(seq,str_seq,accurate_sequence,accurate_seq_where)
                if test_accurate_seq == False:
                    no_pass += 1
                    continue
        
            if filter_direct_upstr == 'Y' and already_counted == False:
                test_direct_upstr = pass_direct_upstr(seq,str_seq,direct_upstr_seq)
                if test_direct_upstr == False:
                    no_pass += 1
                    continue
            
            if filter_direct_downstr == 'Y' and already_counted == False:
                test_direct_downstr = pass_direct_downstr(seq,str_seq,direct_downstr_seq)
                if test_direct_downstr == False:
                    no_pass += 1
                    continue
    
            #trim the sequence: remove everything upstream of the first 4 STRs in tandem and downstream of the last 4 STRs in tandem
            cut_seq = seq[first:last+len(str_seq_four)] 
    
            if filter_interrupters == 'Y':
                test_interrupters = pass_interrupters(cut_seq,str_seq,num_interrupting_bases_allowed)
                if test_interrupters == False:
                    no_pass += 1
                    continue

        # for reads that passed all user-defined filters, count STR
        count = count_str(cut_seq,str_seq)
            
        str_2_write = name + '\t' + str(count) + '\n'
        all_counts.write(str_2_write)
        pass_filter += 1
        if int(fwd_rev) == 16:
            all_revcounts.write(str_2_write)
        else:
            all_fwdcounts.write(str_2_write)

    all_revcounts.close()
    all_fwdcounts.close()
    all_counts.close()

    print('\n-------------------Done! Here are your filtering results-------------------')
    print('Total number of on-target reads from bamfile: '+str(total_reads))
    print('Number of reads passing your filter requests: '+str(pass_filter))
    print('Number of reads that did not pass filter requests: '+str(no_pass))
    print('\n-------------------Your counts are in this folder: ' + output_path + '-------------------\n')


    filter_txt.write('------------All filters used for counts.txt files in this directory-----------\n')
    filter_txt.write(f'input_file = {input_bam}\n')
    filter_txt.write(f'chrom = {chrom}\n')
    filter_txt.write(f'fetch_start = {fetch_start}\n')
    filter_txt.write(f'fetch_end = {fetch_end}\n')
    filter_txt.write(f'STR sequence = {str_seq}\n')
    filter_txt.write(f'Minimum fragment length = {min_fragment_len}\n\n')

    filter_txt.write('## Would you like to filter for reads that include a highly accurate sequence around the STR?\n')
    filter_txt.write(f'filter_accurate_sequence = {filter_accurate_sequence}\n')
    filter_txt.write(f'accurate_sequence = {accurate_sequence}\n')
    filter_txt.write(f'accurate_seq_where = {accurate_seq_where}\n\n')
    
    filter_txt.write('## Would you like to filter for reads that have a particular sequence directly upstream of the STR?\n')
    filter_txt.write(f'filter_direct_upstr = {filter_direct_upstr}\n')
    filter_txt.write(f'direct_upstr_seq = {direct_upstr_seq}\n\n')
    
    filter_txt.write('## Would you like to filter for reads that have a particular sequence directly downstream of the STR?\n')
    filter_txt.write(f'filter_direct_downstr = {filter_direct_downstr}\n')
    filter_txt.write(f'direct_downstr_seq = {direct_downstr_seq}\n\n')

    filter_txt.write('## Would you like to filter out reads that have lots of non-STR bases within the STR?\n')
    filter_txt.write(f'filter_interrupters = {filter_interrupters}\n')
    filter_txt.write(f'num_interrupting_bases_allowed = {num_interrupting_bases_allowed}\n')
    if filter_interrupters == 'Y':
         filter_txt.write(f'Any reads that have more than {num_interrupting_bases_allowed} between two instances of 4 STRs have been removed\n')

    filter_txt.write('\n## Would you like to check if there are normal-length reads that pass all filters above and have <4 repeats?\n')
    filter_txt.write(f'check_for_less_than_four_repeats = {check_for_less_than_four_repeats}\n')
    if check_for_less_than_four_repeats == 'Y':
         filter_txt.write(f'Any reads with 2 or 3 repeats in tandem have been included as a normal-length read \n')

    filter_txt.write('\n-------------------Here are your filtering results-------------------\n')
    filter_txt.write(f'Total number of on-target reads from bamfile: {total_reads}\n')
    filter_txt.write(f'Number of reads passing your filter requests: {pass_filter}\n')
    filter_txt.write(f'Number of reads that did not pass filter requests: {no_pass}')

    filter_txt.close()


def pass_accurate_seq(seq,str_seq,accurate_seq,accurate_seq_where):
    str_seq_four = str_seq*4
    res = False
    find = seq.find(accurate_seq) #find the accurate sequence provided by user
    if find != -1: # if accurate_sequence exists in the read, determine location compared to STR 
        if accurate_seq_where == 'U':
            first = seq.find(str_seq_four) #find first instance of four repeats in tandem
            if find + len(accurate_seq) <= first:
                res = True
        else:
            last = seq.rfind(str_seq_four) #find last instance of four repeats in tandem
            if last + len(str_seq_four) <= find:
                res = True
    return res

def pass_accurate_seq_two(seq,str_seq,accurate_seq,accurate_seq_where):
    str_seq_two = str_seq*2
    res = False
    find = seq.find(accurate_seq) #find the accurate sequence provided by user
    if find != -1: # if accurate_sequence exists in the read, determine location compared to STR 
        if accurate_seq_where == 'U':
            first = seq.find(str_seq_two) #find first instance of two repeats in tandem
            if find + len(accurate_seq) <= first:
                res = True
        else:
            last = seq.rfind(str_seq_two) #find last instance of two repeats in tandem
            if last + len(str_seq_two) <= find:
                res = True
    return res

def pass_direct_upstr(seq,str_seq,direct_upstr_seq):
    str_seq_four = str_seq*4
    first = seq.find(str_seq_four) #find first instance of four consecutive repeats
    res = False
    if first != -1 and first >= len(direct_upstr_seq): 
        res = seq[first-len(direct_upstr_seq):first] == direct_upstr_seq
    return res

def pass_direct_upstr_two(seq,str_seq,direct_upstr_seq):
    str_seq_two = str_seq*2
    first = seq.find(str_seq_two) #find first instance of two consecutive repeats
    res = False
    if first != -1 and first >= len(direct_upstr_seq): 
        res = seq[first-len(direct_upstr_seq):first] == direct_upstr_seq
    return res

def pass_direct_downstr(seq,str_seq,direct_downstr_seq):
    str_seq_four = str_seq*4
    last = seq.rfind(str_seq_four) #find last instance of four consecutive repeats
    res = False
    if last != -1 and last <= len(seq)-len(direct_downstr_seq)-len(str_seq_four): 
        res = seq[last+len(str_seq_four):last+len(str_seq_four)+len(direct_downstr_seq)] == direct_downstr_seq
    return res

def pass_direct_downstr_two(seq,str_seq,direct_downstr_seq):
    str_seq_two = str_seq*2
    last = seq.rfind(str_seq_two) #find last instance of four consecutive repeats
    res = False
    if last != -1 and last <= len(seq)-len(direct_downstr_seq)-len(str_seq_two): 
        res = seq[last+len(str_seq_two):last+len(str_seq_two)+len(direct_downstr_seq)] == direct_downstr_seq
    return res

def pass_interrupters(cut_seq,str_seq,num_interrupting_bases_allowed):
    str_seq_four = str_seq*4
    current_index = 0
    store = []
    write = True
    while current_index != -1:
        store.append(current_index)
        new_index = cut_seq.find(str_seq_four,current_index+1) #keep looking for the next four repeats
        if new_index > current_index+len(str_seq_four)+num_interrupting_bases_allowed: 
            write = False      #change boolean to False if the next four repeats are farther away than the number of bases allowed by user
            break              #also break the while loop 
        current_index = new_index
    
    return write


def count_str(cut_seq,str_seq):
    current_index = 0
    store = []
    while current_index != -1:
        store.append(current_index)
        current_index = cut_seq.find(str_seq,current_index+1) #keep looking for the next repeat
    count = len(store) #calculate how many instances of STR we found

    return count

if __name__ == '__main__':
    try:
        main()
        sys.exit(0)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)




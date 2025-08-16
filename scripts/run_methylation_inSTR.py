'''
MASTR-seq computational pipeline
PART5_run_methylation_inSTR.py
Written by Hanseul Ryu (Aug 2024)

Analyzes CpG methylation within STRs. Outputs .tsv files that store CpG methylation likelihood within STR and STR length for each on-target read. For each on-target read, it also generates a plot showing methylation likelihood and interrupting bases within the STR, base by base, so that users can visually do a quality check.

INPUT: (1) original modified bam file prior to alignment (output of PART1a), (2) fwd counts.txt file from PART2, (3) rev counts.txt file from PART2, (4) STR sequence in forward strand, (5) desired output folder

OUTPUT: .tsv file and plots within desired output folder

IMPORTANT: 
(1) Please make sure you have activated the correct conda environment as specified in README
(2) For input (1), please remember to provide original modified bam file (right after basecalling, before alignment). This way, we read in the CpG methylation likelihood (ML) for each read directly from Dorado basecalling and avoid the risk of reading in bam files where the ML information is changed or missing while aligning or demultiplexing. 
'''

import pysam
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess

def main():

    """
    example command 
    ----------------
    python PART5_run_methylation_inSTR.py [input_bam_prior_to_alignment] [_fwdcounts.txt] [_revcounts.txt] [STR_seq_fwd] [desired_output_dir]

    example command for G2C4 in C9orf72 
    ----------------
    python PART5_run_methylation_inSTR.py example_c9orf72_ALS_basecalled_sup_wCpG_unaligned.bam ./example_c9orf72_ALS_demultiplexed_barcode19_filter1/example_c9orf72_ALS_demultiplexed_barcode19_fwdcounts.txt ./example_c9orf72_ALS_demultiplexed_barcode19_filter1/example_c9orf72_ALS_demultiplexed_barcode19_revcounts.txt GGCCCC ./example_c9orf72_ALS_demultiplexed_barcode19_filter1/methylation_inSTR 
    """

    input_bam = sys.argv[1]
    input_fwd_counts = sys.argv[2]
    input_rev_counts = sys.argv[3]
    STR_entered = sys.argv[4]
    output_folder = sys.argv[5].rstrip('/')

    sorted_input_bam = input_bam.rstrip('.bam') + '_sorted.bam'
    if not os.path.isfile(sorted_input_bam):
        print('First sorting and indexing your original, unaligned bam file...')
    
        cmd = f'samtools sort {input_bam} > {sorted_input_bam}'
        process = subprocess.run(cmd, shell = True, executable="/bin/bash")

        cmd = f'samtools index {sorted_input_bam}'
        process = subprocess.run(cmd,shell = True, executable="/bin/bash")


    if not os.path.exists(output_folder): 
        os.makedirs(output_folder) 
        os.makedirs(f'{output_folder}/plots_for_each_read') 

    read_type_arr = ['fwd','rev']
    count_file_arr = [input_fwd_counts,input_rev_counts]
    base_bam_name = input_bam[input_bam.rfind('/')+1:].rstrip('.bam')

    df = pd.DataFrame(columns =['read_type','name','STR_length','methylation_likelihood'])

    for i in [0,1]:
        input_counts_df = pd.read_csv(count_file_arr[i], sep='\t', names = ['name','count'])

        filtered_bam_name = f'{base_bam_name}_{read_type_arr[i]}_filtered_unaligned.bam'

        print(f'-------------- Now filtering your original, unaligned bam file for {read_type_arr[i]} --------------')
        cmd = f'samtools view -o {output_folder}/{filtered_bam_name} -N {count_file_arr[i]} {input_bam}'
        process = subprocess.run(cmd,shell = True, executable="/bin/bash")

        print(f'-------------- Now analyzing methylation calls within STR for {read_type_arr[i]} --------------')

        bamfile = pysam.AlignmentFile(f'{output_folder}/{filtered_bam_name}', "rb", check_sq=False)

        if read_type_arr[i] == 'rev':
            str_seq = convert_str_seq_fwd_to_rev(STR_entered)
        else:
            str_seq = STR_entered
            
        all_names = [] #store all read names in this bamfile
        for s in bamfile:
            x = str(s)
            x_split = x.split('\t')
            name = x_split[0]
            seq = x_split[9]

            #save mc and ml information from bamfile
            MM_all = s.get_tag('MM')
            start_mc = MM_all.index('C+m')
            mc_str = MM_all[start_mc:].lstrip('C+m?,').rstrip(';')
            mc_arr = mc_str.split(',')
            mc_arr_np = np.array(mc_arr).astype('int') #mc_arr_np stores data about which C's have methylation calls
            ML_all = s.get_tag('ML')
            ml_arr = ML_all[int(len(ML_all)/2):] #ml_arr stores methylation likelihood (from 0 to 255) for each methylation call

            if name in all_names: #there have been instances where the same read is represented twice in some bam files
                continue #skip to next iteration of loop if we already processed the read
                    
            all_names.append(name)

            locate = input_counts_df.loc[input_counts_df['name'] == name]
            count_recorded = int(locate['count'].iloc[0])
            
            if count_recorded >= 4:  #if recorded count is >=4, detect STR using units of 4 repeats in tandem
                #left trim this sequence (i.e. remove everything upstream of STR)
                (seq_ltrim,mc_ltrim,ml_ltrim) = left_trim_read(seq,str_seq,mc_arr_np,ml_arr,4)
                
                #right trim this sequence (i.e. remove everything downstream of STR)
                (seq_rtrim,mc_rtrim,ml_rtrim) = right_trim_read(seq_ltrim,str_seq,mc_ltrim,ml_ltrim,4)

            else: #if recorded count is <4, look for CpG methylation in sequences with 2 repeats in tandem
                #left trim this sequence (i.e. remove everything upstream of STR) by looking for 2 consecutive repeats
                (seq_ltrim,mc_ltrim,ml_ltrim) = left_trim_read(seq,str_seq,mc_arr_np,ml_arr,2)
                
                #right trim this sequence (i.e. remove everything downstream of STR) by looking for 2 consecutive repeats
                (seq_rtrim,mc_rtrim,ml_rtrim) = right_trim_read(seq_ltrim,str_seq,mc_ltrim,ml_ltrim,2)


            #make methylation plot for this read and also store methylation likelihood values within the STR (mC)
            plot_name = f'{output_folder}/plots_for_each_read/{read_type_arr[i]}_{name}.png'
            mC = plot_mC_and_interrupters(seq_rtrim,str_seq,mc_rtrim,ml_rtrim,plot_name)
                
            #add STR methylation info for this read to dataframe
            df = df._append({'read_type':read_type_arr[i], 'name':name, 'STR_length':count_recorded, 'methylation_likelihood': mC},ignore_index = True)
            

    #print final dataframe to tsv file in the same directory as the plots
    tsv_name = output_folder + '/STRlength_5mC-likelihood.tsv'
    df.to_csv(tsv_name, sep="\t") 
    print(f'Done! Please find plots and tsv file for methylation calls within STR at: {output_folder}\n')

    return


def convert_str_seq_fwd_to_rev(str_seq_fwd):
    str_seq_rev = "x" * len(str_seq_fwd)
    str_seq_rev = [None] * len(str_seq_fwd)

    for i in range(0,len(str_seq_fwd)):
        if str_seq_fwd[i] == 'A':
            str_seq_rev[len(str_seq_rev)-i-1] = 'T'
        elif str_seq_fwd[i] == 'C':
            str_seq_rev[len(str_seq_rev)-i-1] = 'G'
        elif str_seq_fwd[i] == 'G':
            str_seq_rev[len(str_seq_rev)-i-1] = 'C'
        elif str_seq_fwd[i] == 'T':
            str_seq_rev[len(str_seq_rev)-i-1] = 'A'
        else:
            raise Exception("Error: STR sequence must only have a combination of A/C/G/T.")
    str_seq_rev = ''.join(str_seq_rev)
    return (str_seq_rev)


def left_trim_read(seq,str_seq,mc_arr,ml_arr,howmany):
    """ discard all info upstream of the repeat
    howmany variable determines whether we detect the STR in units of 4 or 2 repeats in tandem.
    4 is default unless the count is recorded as 2 or 3 repeats """

    find = seq.find(str_seq*howmany) #find first instance of 4 or 2 consecutive repeats
    count_Cs_here = seq[:find] 
    count = 0
    found = 0
    found = count_Cs_here.find("C") #start by finding first C
    while found != -1: 
        count +=1 
        found = count_Cs_here.find("C",found+1) #continue to find next C until you reach the STR

    #here, we find all the methylation calls that correspond to the C's upstream of the STR
    place = 0 
    num_of_Cs_before_STR = mc_arr[place] + 1
    while num_of_Cs_before_STR <= count and place < len(mc_arr)-1:
        place += 1
        num_of_Cs_before_STR = num_of_Cs_before_STR + mc_arr[place] + 1

    #now, we trim out everything upstream of the STR in the original sequence, mc_arr, and ml_arr
    seq_ltrim = seq[find:]
    mc_ltrim = mc_arr[place:]
    ml_ltrim = ml_arr[place:]
    mc_ltrim[0] = mc_ltrim[0]-(count-(sum(mc_arr[:place])+len(mc_arr[:place])))
    
    return (seq[find:],mc_ltrim,ml_ltrim)


def right_trim_read(seq_ltrim,str_seq,mc_ltrim,ml_ltrim,howmany):
    """discard all info downstream of the repeat"""

    last_instance = seq_ltrim.rfind(str_seq*howmany) #find last instance of 4 or 2 consecutive repeats
    seq_rtrim = seq_ltrim[:last_instance+len(str_seq*howmany)] #trim out all the sequences downstream of STR

    #here we count how many C's are in your STR
    count = 0
    found = seq_rtrim.find("C") #first first C in STR
    while found != -1:
        count +=1 
        found = seq_rtrim.find("C",found+1) #continue to find next C until you reach the end of the STR
        
    # now, read mc_ltrim until all C's within the STR have been accounted for
    place = 0 
    num_of_Cs_rec = mc_ltrim[place] + 1 #the first mC call 
    while num_of_Cs_rec <= count and place < len(mc_ltrim)-1:
        place += 1
        num_of_Cs_rec = num_of_Cs_rec + mc_ltrim[place] + 1
    
    #now, we trim out everything downstream of the STR in mc_arr and ml_arr
    mc_rtrim = mc_ltrim[:place]
    ml_rtrim = ml_ltrim[:place]
    return (seq_rtrim,mc_rtrim,ml_rtrim)


def plot_mC_and_interrupters(seq_rtrim,str_seq,mc_rtrim,ml_rtrim,output_name):
    """plots interrupting bases and 5mC methylation call values from trimmed reads
   also returns an array of methylation call values"""

    base = [None] * len(seq_rtrim)
    i = 0
    
    #for the entirety of the STR sequence
    while i < len(seq_rtrim):
        #if the next few bases matches the STR, keep reading 
        if seq_rtrim[i:i+len(str_seq)] == str_seq:
            i = i+len(str_seq)

        # if it doesn't match the STR, add a dot corresponding to the erroneous base (i.e. interrupting base)
        else:
            if seq_rtrim[i] == 'A':
                base[i] = 1.2
            elif seq_rtrim[i] == 'C':
                base[i] = 1.4
            elif seq_rtrim[i] == 'G':
                base[i] = 1.6
            else:
                base[i] = 1.8
            i += 1 #move on to the next base
    
    #make an array as long as the trimmed sequence
    mC = [None] * len(seq_rtrim)
    #mC_onlyCs will only contain methylation-likelihood values (no None)
    mC_onlyCs = [] 

    #in mC array, replace None with the methylation likelihood (ml) value for each C with a methylation call 
    #Note that the modified bam files store methylation likelihood as a value from 0 to 255. Here, we will plot from 0 to 1.
    i = 0
    mc_index = 0
    current_num_of_C = 0
    while i < len(seq_rtrim):
        target_num_of_C = mc_rtrim[mc_index]+1
        if seq_rtrim[i] == 'C':
            current_num_of_C += 1
            if current_num_of_C == target_num_of_C:
                mC[i] = ml_rtrim[mc_index]/255 #add ml value (store methylation likelihood from 0 to 1) 
                mC_onlyCs.append(ml_rtrim[mc_index]/255)
                current_num_of_C = 0
                mc_index += 1
        i += 1

        if mc_index >= len(mc_rtrim):
            break
    
    #determine how many subplots to make
    num_of_plots = len(seq_rtrim)//300+1

    f, axs = plt.subplots(num_of_plots, 1, figsize=(15,3*num_of_plots), constrained_layout=True)
    title_str = output_name.split('/')[-1] + '\n5mC calls in blue\nInterrupting bases in red (A=1.2, C=1.4, G=1.6, T=1.8)'
    f.suptitle(title_str, fontsize = 20)
    f.supxlabel("Base index in STR",fontsize = 20)
    f.supylabel("CpG methylation\nlikelihood",fontsize = 20)
    
    #print interrupting bases and methylation calls along the same x-axis 
    for x in range(num_of_plots):
        ax = plt.subplot(num_of_plots, 1, x+1)
        
        if x == num_of_plots-1:
            ax.plot(range(300*x,len(seq_rtrim)),base[300*x:len(seq_rtrim)],color = 'red',marker='o')
            ax.plot(range(300*x,len(seq_rtrim)),mC[300*x:len(seq_rtrim)],color = 'blue',marker='o')
            ax.set_xlim([300*x,len(seq_rtrim)])

        else:
            ax.plot(range(300*x,300*(x+1)),base[300*x:300*(x+1)],color = 'red',marker='o')
            ax.plot(range(300*x,300*(x+1)),mC[300*x:300*(x+1)],color = 'blue',marker='o')
            ax.set_xlim([300*x,300*(x+1)])

        ax.set_ylim([0, 2])
        ax.axhline(1,color='black',ls='--')
    
    plt.savefig(output_name)
    plt.close()
    
    return (mC_onlyCs) 


if __name__ == '__main__':
    main()
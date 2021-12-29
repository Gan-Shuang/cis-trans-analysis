import re
import sys
import argparse
import pymysql
import pysam
###EGFR cis-trans type list###
cis_trans_list=[(1, 'EGFR', 'p.T790M', 'EGFR', 'p.C797S'),
                (2, 'EGFR', 'p.T790M', 'EGFR', 'p.C797G'),
                (3, 'EGFR', 'p.T790M', 'EGFR', 'p.G796D'), 
                (4, 'EGFR', 'p.T790M', 'EGFR', 'p.L792H'), 
                (5, 'EGFR', 'p.T790M', 'EGFR', 'p.G796R'),
                (6, 'EGFR', 'p.T790M', 'EGFR', 'p.L792F'),
                (7, 'EGFR', 'p.T790M', 'EGFR', 'p.L792P'),
                (8, 'EGFR', 'p.T790M', 'EGFR', 'p.G796S'),
                (9, 'EGFR', 'p.T790M', 'EGFR', 'p.G796A'),
                (10, 'EGFR', 'p.T790M', 'EGFR', 'p.C797Y')]
def get_indel_info(vcf_file,cis_trans_list):
    result_vcf=open(vcf_file,"r")
    mut_genes={}
    pair_mut_phgvs={}
    for info in result_vcf.readlines():
        if not info.startswith("#"):
            gene=str(info).rstrip().split("\t")[7].split("|")[3]
            phgvs=str(info).rstrip().split("\t")[7].split("|")[10]
            chgvs=str(info).rstrip().split("\t")[7].split("|")[9]
            position=str(info).rstrip().split("\t")[0]+":"+str(info).rstrip().split("\t")[1]
            mutation=str(info).rstrip().split("\t")[4]
            mut_genes[position+">"+mutation+"|"+chgvs]=gene+":"+phgvs
    for gene_phgvs in cis_trans_list:
        gene_phgvs_pair1=gene_phgvs[1]+":"+gene_phgvs[2]
        gene_phgvs_pair2=gene_phgvs[3]+":"+gene_phgvs[4]
        for pos,mut in mut_genes.items():
            if gene_phgvs_pair2==mut:
                pair_mut_phgvs[list (mut_genes.keys()) [list (mut_genes.values()).index (gene_phgvs_pair1)]+"-"+pos]=gene_phgvs_pair1+"-"+mut
    return(mut_genes,pair_mut_phgvs)
def analyze_cis_trans(sample_bam,position_list):
    cis_trans_analysis_result={}
    for mut_info,cis_trans_mut in position_list.items():
        print("analysis on event:"+cis_trans_mut+"--"+mut_info)
        ##phgvsA
        mut_name1=cis_trans_mut.split("-")[0]
        mut_position1=int(mut_info.split("-")[0].split("|")[0].split(">")[0].split(":")[1])
        mut_type1=mut_info.split("-")[0].split("|")[0].split(">")[1]
        mut_length1=int(len(str(mut_type1)))
        chr_mut1=int(mut_info.split("-")[0].split("|")[0].split(":")[0])
        ##phgvsB
        mut_name2=cis_trans_mut.split("-")[1]
        mut_position2=int(mut_info.split("-")[1].split("|")[0].split(">")[0].split(":")[1])
        mut_type2=mut_info.split("-")[1].split("|")[0].split(">")[1]
        mut_length2=int(len(str(mut_type2)))
        chr_mut2=int(mut_info.split("-")[1].split("|")[0].split(":")[0])
        ##pysam read bam file
        samfile=pysam.AlignmentFile(sample_bam,"rb")
        phgvsA_list=[]
        phgvsB_list=[]
        for read_info in samfile.fetch(contig=str(chr_mut2), start=int(mut_position2-1), stop=int(mut_position2+mut_length2-1)):
            if int(str(read_info).split('\t')[3]) < mut_position1:
                if read_info.cigar[0][0]==4:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position2-read_start-1+read_info.cigar[0][1])
                        mut_position_end=int(mut_position2-read_start-1+mut_length2+read_info.cigar[0][1])
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if mut_allele==mut_type2 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                            phgvsB_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)==1:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position2-read_start-1)
                        mut_position_end=int(mut_position2-read_start-1+mut_length2)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if mut_allele==mut_type2 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                            phgvsB_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)>1 and read_info.cigar[1][0]==4:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position2-read_start-1)
                        mut_position_end=int(mut_position2-read_start-1+mut_length2)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if mut_allele==mut_type2 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                            phgvsB_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)>1 and read_info.cigar[1][0]==1:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position2-read_start-1)
                        mut_position_end=int(mut_position2-read_start-1+mut_length2)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if read_info.cigar[0][1]>mut_position_start:
                            mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                            if mut_allele==mut_type2 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsB_list.append(str(read_info))
                        if read_info.cigar[0][1]<mut_position_start:
                            insert_size=int(read_info.cigar[1][1])
                            mut_allele=str(read_info).split('\t')[9][mut_position_start+insert_size:mut_position_end+insert_size]
                            if mut_allele==mut_type2 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsB_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)>1 and read_info.cigar[1][0]==2:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position2-read_start-1)
                        mut_position_end=int(mut_position2-read_start-1+mut_length2)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if read_info.cigar[0][1]>mut_position_start:
                            mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                            if mut_allele==mut_type2 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsB_list.append(str(read_info))
                        if read_info.cigar[0][1]<mut_position_start:
                            del_size=int(read_info.cigar[1][1])
                            mut_allele=str(read_info).split('\t')[9][mut_position_start-del_size:mut_position_end-del_size]
                            if mut_allele==mut_type2 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsB_list.append(str(read_info))
                    except IndexError:
                        pass
        for read_info in samfile.fetch(contig=str(chr_mut1), start=int(mut_position1-1), stop=int(mut_position1+mut_length1-1)):
            if str(read_info) in phgvsB_list:
                if read_info.cigar[0][0]==4:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position1-read_start-1+read_info.cigar[0][1])
                        mut_position_end=int(mut_position1-read_start-1+mut_length1+read_info.cigar[0][1])
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if mut_allele==mut_type1 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                            phgvsA_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)==1:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position1-read_start-1)
                        mut_position_end=int(mut_position1-read_start-1+mut_length1)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if mut_allele==mut_type1 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                            phgvsA_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)>1 and read_info.cigar[1][0]==4:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position1-read_start-1)
                        mut_position_end=int(mut_position1-read_start-1+mut_length1)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if mut_allele==mut_type1 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                            phgvsA_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)>1 and read_info.cigar[1][0]==1:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position1-read_start-1)
                        mut_position_end=int(mut_position1-read_start-1+mut_length1)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if read_info.cigar[0][1]>mut_position_start:
                            mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                            if mut_allele==mut_type1 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsA_list.append(str(read_info))
                        if read_info.cigar[0][1]<mut_position_start:
                            insert_size=int(read_info.cigar[1][1])
                            mut_allele=str(read_info).split('\t')[9][mut_position_start+insert_size:mut_position_end+insert_size]
                            if mut_allele==mut_type1 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsA_list.append(str(read_info))
                    except IndexError:
                        pass
                if read_info.cigar[0][0]==0 and len(read_info.cigar)>1 and read_info.cigar[1][0]==2:
                    try:
                        read_start=int(str(read_info).split('\t')[3])
                        read_end=int(str(read_info).split('\t')[7])
                        mut_position_start=int(mut_position1-read_start-1)
                        mut_position_end=int(mut_position1-read_start-1+mut_length1)
                        mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                        if read_info.cigar[0][1]>mut_position_start:
                            mut_allele=str(read_info).split('\t')[9][mut_position_start:mut_position_end]
                            if mut_allele==mut_type1 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsA_list.append(str(read_info))
                        if read_info.cigar[0][1]<mut_position_start:
                            del_size=int(read_info.cigar[1][1])
                            mut_allele=str(read_info).split('\t')[9][mut_position_start-del_size:mut_position_end-del_size]
                            if mut_allele==mut_type1 and int(re.findall('\d+',str(read_info).split('\t')[10])[mut_position_start])>=30 and int(read_info.tags[0][1])<=4:
                                phgvsA_list.append(str(read_info))
                    except IndexError:
                        pass
        list_both=list(set(phgvsA_list))
        print("phgvsB coverage="+str(len(phgvsB_list)))
        print("both mutant coverage="+str(len(phgvsA_list)))
        percent=int(len(phgvsA_list))/int(len(phgvsB_list))
        print("percent_of_"+mut_name2+"="+str(percent))
        cis_trans_analysis_result[str(len(phgvsA_list))+"|"+str(len(phgvsB_list))+"|"+str(len(list_both))+"|"+str(percent)]=cis_trans_mut+"|"+mut_info
    return(cis_trans_analysis_result)
parser=argparse.ArgumentParser(description="Analyze cis-trans from vcf and bam")
parser.add_argument("-v","--vcf",dest="vcf_file",help="input vcf",required=True)
parser.add_argument("-b","--bam",dest="bam_file",help="input bam",required=True)
parser.add_argument("-o","--output",dest="result_file",help="cis-trans analysis result",required=True)
args = parser.parse_args()
sys.stdout.write("Vcf file: "+args.vcf_file+"\n")
sys.stdout.write("Bam file: "+args.bam_file+"\n")
sys.stdout.write("Output: "+args.result_file+"\n")
input_vcf_file=args.vcf_file
input_bam_file=args.bam_file
output_file=args.result_file

if __name__=="__main__":
    vcf_file=input_vcf_file
    bam_file=input_bam_file
    output_result=output_file
    sample_name=vcf_file.split("/")[-1].split(".")[0].split("-")[0]
    print("##########Start analysis on "+sample_name+"##########")
    mut_gene,sample_mut_genes=get_indel_info(vcf_file,cis_trans_list)
    print("Find "+str(len(mut_gene))+" mutations in vcf:"+str(mut_gene))
    sys.stdout.write("Find "+str(len(sample_mut_genes))+" cis-trans event in vcf:")
    print(sample_mut_genes)
    sys.stdout.write("--------------------"+"\n")
    if len(sample_mut_genes)==0:
        sys.stderr.write("Sample find none cis-trans mutation in file:"+vcf_file+"\n")
        sys.stdout.write("Analysis of trans-cis mutaion DONE!"+"\n")
#         sys.exit()
    result=analyze_cis_trans(bam_file,sample_mut_genes)
    sys.stdout.write("--------------------"+"\n")
    sys.stdout.write(sample_name+" result:")
    print(result)
    output=open(output_result,"w")
    output.write("sampleid\tgeneA\tphgvsA\tgeneB\tphgvsB\tcis_trans\tchgvsA\tchgvsB\tratio\n")
    for content,position in result.items():
        persentage=str(content).split("|")[3]
        geneA=str(position).split("|")[0].split("-")[0].split(":")[0]
        geneB=str(position).split("|")[0].split("-")[1].split(":")[0]
        phgvsA=str(position).split("|")[0].split("-")[0].split(":")[1]
        phgvsB=str(position).split("|")[0].split("-")[1].split(":")[1]
        chgvsA=str(position).split("|")[2].split("-")[0]
        chgvsB=str(position).split("|")[3]
        mut=str(position).split("|")[2].split("-")[1]
###persentage threshold set as 0.1###
        if float(persentage) <= 0.1:
            event_type="反式"
            output.write(sample_name+"\t"+geneA+"\t"+phgvsA+"\t"+geneB+"\t"+phgvsB+"\t"+event_type+"\t"+chgvsA+"\t"+chgvsB+"\t"+persentage+"\n")
        else:
            event_type="顺式"
            output.write(sample_name+"\t"+geneA+"\t"+phgvsA+"\t"+geneB+"\t"+phgvsB+"\t"+event_type+"\t"+chgvsA+"\t"+chgvsB+"\t"+persentage+"\n")
    output.close()
    print("##########Finish analysis on "+sample_name+"##########")

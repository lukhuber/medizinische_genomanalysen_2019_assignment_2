#! /usr/bin/env python3

import vcf

__author__ = 'Lukas HUBER'


class Assignment2:
    
    def __init__(self, input1_file, input2_file):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)

        self.vcf1_file = input1_file
        self.vcf2_file = input2_file
        

    def get_average_quality_of_file(self, file):
        '''
        Get the average PHRED quality of all variants
        :return:
        '''
        phred = 0

        for i in vcf.Reader(open(file, "r")):
            phred += i.QUAL

        if file == self.vcf1_file:
            average_phred = round((phred / self.total_number_of_variants1), 2)
        else:
            average_phred = round((phred / self.total_number_of_variants2), 2)

        return average_phred


        
    def get_total_number_of_variants_of_file(self, file):
        '''
        Get the total number of variants
        :return: total number of variants
        '''
        total_number_of_variants = 0

        for i in vcf.Reader(open(file, "r")):
            total_number_of_variants += 1

        return total_number_of_variants
    
    
    def get_variant_caller_of_vcf(self, file):
        '''
        Return the variant caller name
        :return: 
        '''
        variant_callers = set()

        for i in vcf.Reader(open(file, "r")):
            for j in i.INFO["callsetnames"]:
                variant_callers.add(j)

        variant_callers.remove('')

        return sorted(variant_callers)
        
        
    def get_human_reference_version(self, file):
        '''
        Return the genome reference version
        :return: 
        '''
        human_reference_versions = set()

        for i in vcf.Reader(open(file, "r")):
            if 'difficultregion' in i.INFO.keys():
                for human_reference_version in i.INFO['difficultregion']:
                        human_reference_versions.add(human_reference_version)

        return sorted(human_reference_versions)
        
        
    def get_number_of_indels(self, file):
        '''
        Return the number of identified INDELs
        :return:
        '''
        number_of_indels = 0

        for i in vcf.Reader(open(file, "r")):
            if i.is_indel:
                number_of_indels += 1

        return number_of_indels
        

    def get_number_of_snvs(self, file):
        '''
        Return the number of SNVs
        :return: 
        '''
        number_of_snvs = 0

        for i in vcf.Reader(open(file, "r")):
            if i.is_snp:
                number_of_snvs += 1

        return number_of_snvs
        


    def get_number_of_heterozygous_variants(self, file):
        '''
        Return the number of heterozygous variants
        :return: 
        '''
        number_of_heterozygous_variants = 0

        for i in vcf.Reader(open(file, "r")):
            number_of_heterozygous_variants += i.num_het

        return number_of_heterozygous_variants
        
    
    def merge_chrs_into_one_vcf(self, file1, file2):
        '''
        Creates one VCF containing all variants of chr21 and chr22
        :return:
        '''

        print("Merging chr21_new.vcf with chr22_new.vcf")

        vcf_file1 = vcf.Reader(open(file1), "r")
        vcf_file2 = vcf.Reader(open(file2), "r")

        vcf_writer = vcf.Writer(open("merged_file.vcf", "w"),   vcf_file1)

        for file in [vcf_file1, vcf_file2]:
            for i in file:
                vcf_writer.write_record(i)

        print("Merge successful. File 'merged_file.vcf' created.")
    
    def print_summary(self):
        print("Calculating total number of variants")
        self.total_number_of_variants2 = self.get_total_number_of_variants_of_file(self.vcf2_file)
        print("Calculating average qualitiy of file")
        self.average_quality_of_file2 = self.get_average_quality_of_file(self.vcf2_file)
        print("Getting variant caller")
        self.variant_caller2 = self.get_variant_caller_of_vcf(self.vcf2_file)
        print("Getting human reference version")
        self.human_reference_version2 = self.get_human_reference_version(self.vcf2_file)
        print("Calculating number of indels")
        self.number_of_indels2 = self.get_number_of_indels(self.vcf2_file)
        print("Calculating number of snvs")
        self.number_of_snvs2 = self.get_number_of_snvs(self.vcf2_file)
        print("Calculating number of heterozygous variants")
        self.number_of_heterozygous_variants2 = self.get_number_of_heterozygous_variants(self.vcf2_file)
        
        print("\n------------------- RESULTS -------------------")
        print("Total number of variants of chr22:       ", self.total_number_of_variants2)
        print("Average quality of chr22:                ", self.average_quality_of_file2)
        print("Number of indels in chr22:               ", self.number_of_indels2)
        print("Number of snvs in chr22:                 ", self.number_of_snvs2)
        print("Number of heterozygous variants in chr22:", self.number_of_heterozygous_variants2)
        print("Variant caller of chr22:                ")
        for caller in self.variant_caller2:
            print("                          ", caller)
        print("Human Reference Version of chr22:        ")
        for ref in self.human_reference_version2:
            print(" ", ref)
        print("-----------------------------------------------\n")

        self.merge_chrs_into_one_vcf(self.vcf1_file, self.vcf2_file)
    
    
def main():
    vcf1_file = "chr21_new.vcf"
    vcf2_file = "chr22_new.vcf"

    print("Assignment 2")
    assignment2 = Assignment2(vcf1_file, vcf2_file)
    assignment2.print_summary()
    print("Done with assignment 2")
        
        
if __name__ == '__main__':
    main()
   
    




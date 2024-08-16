import os
import subprocess
import tskit
import sys
from implantErrors import implant_error, error_imputer



class stdpopsim:
    def __init__(self, output_file_name: str, chromosome = "20", afr_pop = 500, eur_pop = 0, genetic_map_type = "HapMapII_GRCh38", seed = 1):
        self.chromosome = chromosome
        self.afr_pop = afr_pop
        self.eur_pop = eur_pop
        self.genetic_map_type = genetic_map_type
        self.output_file_name = output_file_name
        self.sim_run = False
        self.seed = seed

    def run(self):
        if os.path.splitext(self.output_file_name)[-1] != ".trees":
            self.output_file_name = self.output_file_name + ".trees"
        command = ["stdpopsim", "HomSap", "-c", f"{self.chromosome}", "-g", f"{self.genetic_map_type}", "-d", "OutOfAfrica_2T12", f"AFR:{self.afr_pop}", f"EUR:{self.eur_pop}",
                   "-o", self.output_file_name, "-s", f"{self.seed}"]
        subprocess.run(command)
        self.sim_run = True
        
                

    def toVCF(self):
        """
        extracts a vcf file from the tree sequence
        do not run if the simulation has not been run yet, as there will not be a tree sequence file to extract from
        """
        if self.sim_run == False:
            print("Simulation has not been run yet")
            return
        vcf_file = os.path.splitext(self.output_file_name)[0] + ".vcf"
        ts = tskit.load(self.output_file_name)
        f_vcf = open(vcf_file, "w+")
        ts.write_vcf(f_vcf, contig_id=self.chromosome)
        f_vcf.close()
        return vcf_file
    

def main():
    """
    Usage: 
    ts_file name: the name that you want to assign to your tree sequence
    sim_population_size: the number of individuals you want to simulate; note that the number of samples refers to the number of haplotypes in the cohort, so that number will be double this number
    maf_cutoff: the desired maximum minor allele frequency, given in decimal (not percent)
    genotyping_error_rate: the desired genotyping error to be implanted, given in decimal (not percent)
    """
    ts_file_name = sys.argv[1]
    if ts_file_name[len(ts_file_name) - 6:] != '.trees.':
        ts_file_name += '.trees'
    sim_population_size = sys.argv[2]
    maf_cutoff = float(sys.argv[3])
    genotyping_error_rate = float(sys.argv[4])
    random_seed = int(sys.argv[5])
    sim = stdpopsim(ts_file_name, afr_pop=sim_population_size, seed = random_seed)
    print("runnning simulation...")
    sim.run()
    print("extracting VCF from tree sequence...")
    vcf_file_path = sim.toVCF()
    print(f"Processing VCF File (removing multiallelic sites and filtering by minor allele frequency of {maf_cutoff * 100}%)...")
    with open("tmp_file", "w+") as file:
        subprocess.run(['bcftools', 'view', '-m2', '-M2', '-v', 'snps', '-q', f'{maf_cutoff}:minor', f"{vcf_file_path}"], stdout = file)
    vcf_file_path = os.path.splitext(vcf_file_path)[0] + f"_maf{maf_cutoff}" + ".vcf"
    subprocess.run(["mv", "tmp_file", vcf_file_path])
    vcf_file_path_e = os.path.splitext(vcf_file_path)[0] + f"_e{genotyping_error_rate}" + ".vcf"
    print(f"Implanting genotyping error of {genotyping_error_rate * 100}%...")
    implant_error(vcf_file_path, genotyping_error_rate, vcf_file_path_e)
    subprocess.run(["rm", vcf_file_path, "error_loci.txt", ts_file_name[0:len(ts_file_name) - 6] + ".vcf"])
    print(" ")
    print(vcf_file_path_e)

    
    

if __name__ == "__main__":
    main()

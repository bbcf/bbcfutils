package ch.epfl.bbcf.bbcfutils.access.genrep;

import java.io.IOException;
import java.util.List;


import ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo.Assembly;
import ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo.Chromosome;
import ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo.Genome;
import ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo.Organism;




public class Examples {

	/**
	 * Example method on how to use this library
	 */
	public static void main(String[] args) {

		try {
			//fetch the mus musculus nr assembly (mm9)
			//NR_Assembly nrAssembly = GenrepWrapper.getNrAssemblyFromId(70);
			Assembly assembly = GenrepWrapper.getAssemblyFromNrAssemblyId(158);
			//print different info on the assembly
			System.out.println("assembly "+assembly.getName()+" from "+assembly.getSource_name());
			//print more information on the genome associated with this assembly
			Genome genome = GenrepWrapper.getGenomeFromId(assembly.getGenome_id());
			System.out.println(" genome "+genome.getName()+" with tax_id "+genome.getTax_id());
			//get the chromosomes of the assembly
			List<Chromosome> chromosomes = assembly.getChromosomes();
			//print the chromosomes name and id
			for(Chromosome chromosome : chromosomes){
				System.out.println("\tchromosome "+chromosome.getName()+" with id "+chromosome.getId()+"  "+chromosome.getRefseq_version()+"  "+chromosome.getRefseq_locus());
			}
			//I have a BED file with two chromosomes in input : chrI and an NC number NC_000069.5
			//the first one is obviously the chromosome 1, but let's convert it automatically. And
			// I don't know about the second.
			System.out.println("guess chr1");
			Chromosome one = GenrepWrapper.guessChromosome("chr1",assembly.getId());
			if(one!=null){
				System.out.println("chromosome "+one.getName());
			} else {
				System.err.println("chromosome chr1 not found");
			}
			System.out.println("guess NC_000069.5");
			Chromosome iDontKnow = GenrepWrapper.guessChromosome("NC_000069.5",assembly.getId());
			if(iDontKnow!=null){
				System.out.println("chromosome : "+iDontKnow.getName());
			} else {
				System.err.println("chromosome NC_000069.5 not found");
			}
			//get all organisms
			List<Organism> organisms = GenrepWrapper.getOrganisms();
			for(Organism o : organisms){
				System.out.println(o.getSpecies());
			}

		} catch (IOException e) {
			e.printStackTrace();
		} catch (MethodNotFoundException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	}



}

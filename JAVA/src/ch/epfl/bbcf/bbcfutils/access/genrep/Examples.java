package ch.epfl.bbcf.bbcfutils.access.genrep;

import java.io.IOException;
import java.util.List;

import org.codehaus.jackson.type.TypeReference;

import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.FORMAT;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.KEY;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.METHOD;
import ch.epfl.bbcf.bbcfutils.access.genrep.pojo.Assembly;
import ch.epfl.bbcf.bbcfutils.access.genrep.pojo.Chromosome;
import ch.epfl.bbcf.bbcfutils.access.genrep.pojo.Genome;
import ch.epfl.bbcf.bbcfutils.access.genrep.pojo.GenrepObject;
import ch.epfl.bbcf.bbcfutils.access.genrep.pojo.NR_Assembly;
import ch.epfl.bbcf.bbcfutils.access.genrep.pojo.Organism;




public class Examples {

	/**
	 * Example method on how to use this library
	 */
	public static void main(String[] args) {

		try {
			//fetch the mus musculus nr assembly (mm9)
			NR_Assembly nrAssembly = getNrAssemblyFromId(70);
			Assembly assembly = getAssemblyFromNrAssemblyId(70);
			//print different info on the assembly
			System.out.println("assembly "+assembly.getName()+" from "+assembly.getSource_name());
			//print more information on the genome associated with this assembly
			Genome genome = getGenomeFromId(assembly.getGenome_id());
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
			Chromosome one = guessChromosome("chr1",assembly.getId());
			if(one!=null){
				System.out.println("chromosome "+one.getName());
			} else {
				System.err.println("chromosome chr1 not found");
			}
			System.out.println("guess NC_000069.5");
			Chromosome iDontKnow = guessChromosome("NC_000069.5",assembly.getId());
			if(iDontKnow!=null){
				System.out.println("chromosome : "+iDontKnow.getName());
			} else {
				System.err.println("chromosome NC_000069.5 not found");
			}
			//get all organisms
			List<Organism> organisms = getOrganisms();
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



	private static Chromosome guessChromosome(String guessable, int assemblyId) throws MethodNotFoundException, IOException {
		@SuppressWarnings("unchecked")
		List<Chromosome> chromosomes = (List<Chromosome>) GenRepAccess.doQueryList(
				Constants.URL, METHOD.INDEX, FORMAT.json, new TypeReference<List<GenrepObject>>() {},KEY.chromosomes,
				"assembly_id="+assemblyId+"&identifier="+guessable);
		if(!chromosomes.isEmpty()){
			return chromosomes.get(0);
		}
		return null;
	}



	private static Assembly getAssemblyFromNrAssemblyId(int id) throws MethodNotFoundException, IOException {
		@SuppressWarnings("unchecked")
		List<Assembly> assemblies = (List<Assembly>) GenRepAccess.doQueryList(
				Constants.URL, METHOD.ALL, FORMAT.json, new TypeReference<List<GenrepObject>>() {},KEY.assemblies,null);
		for(Assembly assembly : assemblies){
			if(assembly.getNr_assembly_id()==id){
				return assembly;
			}
		}
		return null;
	}



	private static NR_Assembly getNrAssemblyFromId(int id) throws IOException, MethodNotFoundException, ClassNotFoundException {
		return (NR_Assembly) GenRepAccess.doQuery(
				Constants.URL, METHOD.SHOW, FORMAT.json, NR_Assembly.class,KEY.nr_assemblies, id);
	}



	/**
	 * Get an Assembly from it's id
	 * @param id - the id of the assembly
	 * @throws IOException
	 * @throws MethodNotFoundException
	 * @throws ClassNotFoundException
	 * @return an Assembly
	 */
	public static Assembly getAssemblyFromId(int id) throws IOException, MethodNotFoundException, ClassNotFoundException{
		return (Assembly) GenRepAccess.doQuery(
				Constants.URL, METHOD.SHOW, FORMAT.json, Assembly.class,KEY.assemblies, id);
	}

	/**
	 * Get a genome from it's id
	 * @param id - the genome identifier
	 * @return a Genome
	 * @throws ClassNotFoundException 
	 * @throws MethodNotFoundException 
	 * @throws IOException 
	 */
	private static Genome getGenomeFromId(int id) throws IOException, MethodNotFoundException, ClassNotFoundException {
		return (Genome) GenRepAccess.doQuery(
				Constants.URL, METHOD.SHOW, FORMAT.json, Genome.class,KEY.genomes, id);
	}

	/**
	 * Get all organisms present in Genrep
	 * @return A List of Organism
	 * @throws IOException 
	 * @throws MethodNotFoundException 
	 */
	@SuppressWarnings("unchecked")
	private static List<Organism> getOrganisms() throws MethodNotFoundException, IOException {
		return (List<Organism>) GenRepAccess.doQueryList(
				Constants.URL, METHOD.ALL, FORMAT.json, new TypeReference<List<GenrepObject>>() {},KEY.organisms,null);
	}


}

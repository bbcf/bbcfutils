package ch.epfl.bbcf.access.genrep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.codehaus.jackson.JsonNode;
import org.codehaus.jackson.map.DeserializationConfig;
import org.codehaus.jackson.map.DeserializationConfig.Feature;
import org.codehaus.jackson.map.ObjectMapper;

import com.sun.tools.example.debug.bdi.MethodNotFoundException;

import ch.epfl.bbcf.access.InternetConnection;
import ch.epfl.bbcf.access.genrep.Constants.FORMAT;
import ch.epfl.bbcf.access.genrep.Constants.KEY;
import ch.epfl.bbcf.access.genrep.Constants.METHOD;
import ch.epfl.bbcf.access.genrep.pojo.Assembly;
import ch.epfl.bbcf.access.genrep.pojo.Chromosome;
import ch.epfl.bbcf.access.genrep.pojo.Genome;
import ch.epfl.bbcf.access.genrep.pojo.GenrepObject;
import ch.epfl.bbcf.access.genrep.pojo.Organism;
import ch.epfl.bbcf.access.genrep.pojo.Source;

public class Examples {

	/**
	 * Example method on how to use this library
	 */
	public static void main(String[] args) {

		try {
			//fetch the Assembly mm9
			Assembly assembly = getAssemblyFromId(7);
			//print different info on the assembly
			System.out.println("assembly "+assembly.getName()+" from "+assembly.getSource_name());
			//print more information on the genome associated with this assembly
			Genome genome = getGenomeFromId(assembly.getGenome_id());
			System.out.println(" genome "+genome.getName()+" with tax_id "+genome.getTax_id());
			//get the chromosomes of the assembly
			List<Chromosome> chromosomes = assembly.getChromosomes();
			//print the chromosomes name and id
			for(Chromosome chromosome : chromosomes){
				System.out.println("\tchromosome "+chromosome.getName()+" with id "+chromosome.getId());
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (MethodNotFoundException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
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
	



	
}

package ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo;

import org.codehaus.jackson.annotate.JsonProperty;

import ch.epfl.bbcf.bbcfutils.access.genrep.Constants;



public class GenrepObject {

	protected GenrepObject instance;
	
	/**
	 * return the inner class of this Object
	 * DO NOT CALL THIS METHOD
	 * @return
	 */
	public GenrepObject getInstance() {
		return this.instance;
	}

	private void setInstance(GenrepObject instance) {
		this.instance = instance;
	}

	@SuppressWarnings("unused")
	@JsonProperty(Constants.ASSEMBLY)
	private void setAssembly(Assembly instance){
		setInstance(instance);
	}

	@SuppressWarnings("unused")
	@JsonProperty(Constants.CHROMOSOME)
	private void setChromosome(Chromosome instance){
		setInstance(instance);
	}

	@SuppressWarnings("unused")
	@JsonProperty(Constants.ORGANISM)
	private void setOrganism(Organism instance){
		setInstance(instance);
	}

	@SuppressWarnings("unused")
	@JsonProperty(Constants.GENOME)
	private void setGenome(Genome instance){
		setInstance(instance);
	}

	@SuppressWarnings("unused")
	@JsonProperty(Constants.NR_ASSEMBLY)
	private void setNrAssembly(NR_Assembly instance){
		setInstance(instance);
	}


	@SuppressWarnings("unused")
	@JsonProperty(Constants.SOURCE)
	private void setSource(Source instance){
		setInstance(instance);
	}


	//	public static <T> T getInstance(T object,Class<T> clazz){
	//		GenrepObject genrepObject = (GenrepObject) object;
	//		return (T) genrepObject.getInstance();
	//	}
	////	public static <T> T<? extends GenrepObject> getInstance(T object,Class<T> clazz) throws ClassNotFoundException{
	//		GenrepObject genrepObject = (GenrepObject) object;
	//		if(Assembly.class.equals(clazz)){
	//			return  genrepObject.getAssembly();
	//		} else if(Chromosome.class.equals(clazz)){
	//			return (T) genrepObject.getChromosome();
	//		} else if(NR_Assembly.class.equals(clazz)){
	//			return (T) genrepObject.getNRAssembly();
	//		} else if(Genome.class.equals(clazz)){
	//			return (T) genrepObject.getGenome();
	//		} else if(Organism.class.equals(clazz)){
	//			return (T) genrepObject.getOrganism();
	//		} else if(Source.class.equals(clazz)){
	//			return (T) genrepObject.getSource();
	//		} else {
	//			throw new ClassNotFoundException(" The class "+clazz+" is not an GerepObject");
	//		}
	//
	//	}
	//	public static Chromosome getInstance(GenrepObject o, Class<Genome> class1) {
	//		if(Assembly.class.equals(clazz)){
	//			return (T) genrepObject.getAssembly();
	//		} else if(Chromosome.class.equals(clazz)){
	//			return (T) genrepObject.getChromosome();
	//		} else if(NR_Assembly.class.equals(clazz)){
	//			return (T) genrepObject.getNRAssembly();
	//		} else if(Genome.class.equals(clazz)){
	//			return (T) genrepObject.getGenome();
	//		} else if(Organism.class.equals(clazz)){
	//			return (T) genrepObject.getOrganism();
	//		} else if(Source.class.equals(clazz)){
	//			return (T) genrepObject.getSource();
	//		} else {
	//			throw new ClassNotFoundException(" The class "+clazz+" is not an GerepObject");
	//		}
	//	}

}

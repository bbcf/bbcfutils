package ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo;

public class Source extends GenrepObject{

	private String name;
	private int id;
	
	public void setName(String name) {
		this.name = name;
	}
	public String getName() {
		return name;
	}
	public void setId(int id) {
		this.id = id;
	}
	public int getId() {
		return id;
	}
	@Override
	public Source getInstance() {
		return (Source) this.instance;
	}
}

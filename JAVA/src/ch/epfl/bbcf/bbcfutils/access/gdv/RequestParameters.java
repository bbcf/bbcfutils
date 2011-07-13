package ch.epfl.bbcf.bbcfutils.access.gdv;



public class RequestParameters {

	public static enum CMD {new_project,new_track,status,assemblies};

	public static enum PARAM {id,mail,key,command,
		name,seq_id,Public,project_id,
		url,job_id,assemblies,gdv_url}

	private String url,projectId,
	sequenceId,name,isPublic,job_id,
	assemblies,id,gdv_url;

	private String mail,key;

	private CMD command;

	public void setProjectId(String projectId) {
		this.projectId = projectId;
	}
	public String getProjectId() {
		return projectId;
	}
	public void setSequenceId(String sequenceId) {
		this.sequenceId = sequenceId;
	}
	public String getSequenceId() {
		return sequenceId;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getName() {
		return name;
	}
	public void setUrl(String url) {
		this.url = url;
	}
	public String getUrl() {
		return url;
	}
	public void setMail(String mail) {
		this.mail = mail;
	}
	public String getMail() {
		return mail;
	}
	public void setKey(String key) {
		this.key = key;
	}
	public String getKey() {
		return key;
	}
	public void setCommand(String command) {
		try {
			this.command = CMD.valueOf(command);
		} catch(IllegalArgumentException e){}
	}
	public CMD getCommand() {
		return command;
	}
	public void setIsPublic(String isPublic) {
		this.isPublic = isPublic;
	}
	public String getIsPublic() {
		return isPublic;
	}
	public void setAssemblies(String assemblies) {
		this.assemblies = assemblies;
	}
	public String getAssemblies() {
		return assemblies;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getId() {
		return id;
	}
	public void setJob_id(String job_id) {
		this.job_id = job_id;
	}
	public String getJob_id() {
		return job_id;
	}
	public void setGdv_url(String gdv_url) {
		this.gdv_url = gdv_url;
	}
	public String getGdv_url() {
		return gdv_url;
	}
}

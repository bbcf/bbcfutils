package ch.epfl.bbcf.access.gdv;

import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class RequestParameters {

	public static final String NEW_PROJECT_COMMAND ="new_project";
	public static final String ADD_TRACK_COMMAND ="add_track";
	public static final String ADD_SQLITE_COMMAND ="add_sqlite";

	public static final String[] commands = 
	{NEW_PROJECT_COMMAND,ADD_TRACK_COMMAND,ADD_SQLITE_COMMAND};

	public static final String 
	COMMAND_PARAM = "command",
	MAIL_PARAM = "mail",
	KEY_PARAM = "key",
	PROJECT_ID_PARAM = "project_id",
	DATATYPE_PARAM="datatype",
	TYPE_PARAM="type",
	OBFUSCATED_PARAM="obfuscated",
	SEQUENCE_ID_PARAM="seq_id",
	NAME_PARAM="name",
	URL_PARAM="url",
	PUBLIC_PARAM="public";
	
	private static final Map<String,String[]> mapCommands = buildMapCommands();

	private String url,projectId,
	datatype,type,obfuscated,
	sequenceId,name,isPublic;
	
	private String mail,key,command;


	

	/**
	 * method to tell the user
	 * what a command need for arguments
	 * in order to work
	 * @return Map
	 */
	private static Map<String,String[]> buildMapCommands() {
		Map<String,String[]> map = new HashMap<String,String[]>();
		String required[] = {MAIL_PARAM,KEY_PARAM,COMMAND_PARAM};
		String newProject[] = {SEQUENCE_ID_PARAM,NAME_PARAM};
		String groupProject[] = {OBFUSCATED_PARAM};
		String publicProject[] = {PUBLIC_PARAM};
		String addTrack[] = {URL_PARAM,PROJECT_ID_PARAM};
		String addSqlite[] = {URL_PARAM,PROJECT_ID_PARAM,DATATYPE_PARAM};
		map.put("Required for login",required);
		map.put(NEW_PROJECT_COMMAND,newProject);
		map.put(NEW_PROJECT_COMMAND+" public (optionnal)",publicProject);
		map.put(ADD_TRACK_COMMAND,addTrack);
		map.put(ADD_SQLITE_COMMAND,addSqlite);
		return map;
	}



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



	public void setType(String type) {
		this.type = type;
	}



	public String getType() {
		return type;
	}



	public void setObfuscated(String obfuscated) {
		this.obfuscated = obfuscated;
	}



	public String getObfuscated() {
		return obfuscated;
	}



	public void setName(String name) {
		this.name = name;
	}



	public String getName() {
		return name;
	}



	public void setDatatype(String datatype) {
		this.datatype = datatype;
	}



	public String getDatatype() {
		return datatype;
	}



	public void setUrl(String url) {
		this.url = url;
	}



	public String getUrl() {
		return url;
	}

	public static Map<String,String[]> getMapcommands() {
		return mapCommands;
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
		this.command = command;
	}

	public String getCommand() {
		return command;
	}



	public void setIsPublic(String isPublic) {
		this.isPublic = isPublic;
	}



	public String getIsPublic() {
		return isPublic;
	}


	

}

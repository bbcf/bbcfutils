package ch.epfl.bbcf.access.gdv;

import java.io.IOException;

import ch.epfl.bbcf.access.InternetConnection;




public abstract class Command {


	private String command;
	private String mail;
	private String pass;

	protected Command(String command,String mail,String pass){
		this.command = command;
		this.mail = mail;
		this.pass = pass;
	}
	/**
	 * prepare the parameters to send to GDV
	 * @param id
	 * @param params
	 */
	protected abstract boolean doRequest(String id, RequestParameters params);

	/**
	 * send a post to GDV application
	 * @param body
	 * @return
	 */
	protected boolean send(String body2){
		String body = "mail="+mail+"&key="+pass+"&command="+command+"&"+body2;
		try {
			System.out.println(body);
			System.out.println(InternetConnection.sendPOSTConnection(PostToGDV.GDV_ADRESS, body,InternetConnection.MIME_TYPE_FORM_APPLICATION));
			return true;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return false;
	}



	/**
	 * Creation of a new project on GDV
	 * 
	 * @author Yohan Jarosz
	 *
	 */
	public static class Project extends Command{
		protected Project(String command,String mail,String pass){
			super(command, mail, pass);
		}
		@Override
		protected boolean doRequest(String id, RequestParameters params) {
			if(PostToGDV.checkParams(params.getSequenceId(),params.getName())){
				String body = "id="+id
				+"&"+RequestParameters.SEQUENCE_ID_PARAM+"="+params.getSequenceId()
				+"&"+RequestParameters.NAME_PARAM+"="+params.getName();
				if(null!=params.getIsPublic()){
					body+="&"+RequestParameters.PUBLIC_PARAM+"="+params.getIsPublic();
				}
				return send(body);
			}
			return false;
		}


	}

	/**
	 * Adding a track on GDV, already in SQLite format
	 * 
	 * @author Yohan Jarosz
	 *
	 */
	public static class SQLite extends Command{
		protected SQLite(String command,String mail,String pass){
			super(command, mail, pass);
		}

		@Override
		protected boolean doRequest(String id, RequestParameters params) {
			if(PostToGDV.checkParams(params.getUrl(),params.getProjectId(),params.getDatatype())){
				String body = "id="+id+"&"+RequestParameters.URL_PARAM+"="+params.getUrl()
				+"&"+RequestParameters.PROJECT_ID_PARAM+"="+params.getProjectId()
				+"&"+RequestParameters.DATATYPE_PARAM+"="+params.getDatatype();
				if(null!=params.getName()){
					body+="&"+RequestParameters.NAME_PARAM+"="+params.getName();
				}
				return send(body);
			}
			return false;

		}

	}





	/**
	 * Adding a new track on GDV
	 * 
	 * @author Yohan Jarosz
	 *
	 */
	public static class Track extends Command{
		protected Track(String command,String mail,String pass){
			super(command, mail, pass);
		}

		@Override
		protected boolean doRequest(String id, RequestParameters params) {
			if(PostToGDV.checkParams(params.getUrl(),params.getProjectId())){
				String body = "id="+id+"&"+RequestParameters.URL_PARAM+"="+params.getUrl()
				+"&"+RequestParameters.PROJECT_ID_PARAM+"="+params.getProjectId();
				if(null!=params.getName()){
					body+="&"+RequestParameters.NAME_PARAM+"="+params.getName();
				}
				return send(body);
			}
			return false;
		}
	}




	/**
	 * return the right method for the command given
	 * @param command - the command 
	 * @param pass - the password
	 * @param mail - the mail
	 * @return
	 */
	public static Command getControl(String command, String mail, String pass) {
		if(command.equalsIgnoreCase(RequestParameters.NEW_PROJECT_COMMAND)){
			return new Project(command,mail,pass);
		} else if(command.equalsIgnoreCase(RequestParameters.ADD_SQLITE_COMMAND)){
			return new SQLite(command,mail,pass);
		} else if(command.equalsIgnoreCase(RequestParameters.ADD_TRACK_COMMAND)){
			return new Track(command,mail,pass);
		} else {
			return null;
		}
	}



}


package ch.epfl.bbcf.bbcfutils.access.gdv;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import ch.epfl.bbcf.bbcfutils.access.gdv.RequestParameters.CMD;
import ch.epfl.bbcf.bbcfutils.access.gdv.RequestParameters.PARAM;


/**
 * 
 * @author Yohan Jarosz
 *
 */
public class PostToGDV {

	public static final String DEFAULT_GDV_ADRESS = "http://salt.epfl.ch/gdv/post";
	public static final String DEFAULT_GDV_ID = "gdv_post";

	/**
	 * Main entry of the program
	 * @param args
	 */
	public static void main(String[] args) {
		if(!parseArgs(args)){
			usage();
		}
	}




	/**
	 * parse the arguments provided by the user
	 * @param args - the args
	 * @param mail - the login
	 * @param pass - the password
	 * @return
	 */
	private static boolean parseArgs(String[] args) {
		/* help section */
		if(args.length>0){
			if(args[0].equalsIgnoreCase("h")||args[0].equalsIgnoreCase("help")){
				System.out.println(README);	
				return true;
			} else if(args[0].equalsIgnoreCase("gendoc")){
				File file = new File("README");
				try {
					FileWriter writer = new FileWriter(file);
					writer.write(README);
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
				return true;
			}
		}

		/* get all parameters */
		RequestParameters params = buildRequestParameters(args);
		/* check default values */
		if(null==params.getGdv_url()){
			params.setGdv_url(DEFAULT_GDV_ADRESS);
		}
		if(null==params.getId()){
			params.setId(DEFAULT_GDV_ID);
		}
		if(null!=params){
			boolean ok = true;
			ok = Command.checkParam(PARAM.mail,params.getMail());
			ok = Command.checkParam(PARAM.key,params.getKey());
			ok = Command.checkCommand(params.getCommand());
			if(ok){
				return Command.doRequest(params.getGdv_url(),params);
			}
		}

		/* display usage */
		return false;
	}





	/**
	 * build an RequestParamters object with all args provided 
	 * by user
	 * @param args - the arguments
	 * @return RequestParamters
	 */
	private static RequestParameters buildRequestParameters(String[] args) {
		RequestParameters params = new RequestParameters();
		for(int i=0;i<args.length;i++){
			if(i+1<args.length){
				if(args[i].startsWith("--")){
					String cmd = args[i].substring(2, args[i].length());
					PARAM p = PARAM.valueOf(cmd);
					switch(p){
					case project_id :
						params.setProjectId(args[i+1]); 
						break;
					case url : 
						params.setUrl(args[i+1]);
						break;
					case command :
						params.setCommand(args[i+1]);
						break;
					case key :
						params.setKey(args[i+1]);
						break;
					case mail :
						params.setMail(args[i+1]);
						break;
					case seq_id :
						params.setSequenceId(args[i+1]);
						break;
					case name :
						params.setName(args[i+1]);
						break;
					case Public :
						params.setIsPublic(args[i+1]);
						break;
					case assemblies:
						params.setAssemblies(args[i+1]);
						break;
					case id:
						params.setId(args[i+1]);
						break;
					case job_id:
						params.setJob_id(args[i+1]);
						break;
					case gdv_url :
						params.setGdv_url(args[i+1]);
						break;
					}
				} else {
					System.err.println("the parameter "+args[i]+" is not a parameter.");
				}
			}
			i++;
		} 
		return params;
	}







	/**
	 * check if a parameter required is given or not
	 * @param params
	 * @return
	 */
	public static boolean checkParams(String... params) {
		for(String p : params){
			if(null==p){
				System.err.println("\nmissing param(s)");
				return false;
			}
		}
		return true;
	}



	/**
	 * convenient method to write a table to the output
	 * @param values
	 * @return
	 */
	private static String writeList(String[] values) {
		String str = "";
		for(String v :values){
			str+=v+", ";
		}
		str=str.substring(0, str.length()-2);
		return str;
	}

	/**
	 * print help to the console output
	 */
	private static void usage() {
		System.out.println(README);
	}

	private static final String eol = System.getProperty("line.separator");
	private static final String README =
		eol+"################" +
		eol+"#### README ####" +
		eol+"################" +
		eol+"" +
		eol+"This utility will help send jobs to GDV" +
		eol+"" +
		eol+"First you need a mail and an user key." +
		eol+"Login in GDV and go to 'preferences' page you will have your personal user key." +
		eol+"" +
		eol+"### COMMANDS ###" +
		eol+"" +
		eol+"You can launch several commands : " +
		eol+"- new_project : create a new project on GDV " +
		eol+"    return {project_id:<>,public_url:<>}" +
		eol+"- new_track : add a track to a GDV project" +
		eol+"    return {job_id:<>}" +
		eol+"- assemblies : fetch all assemblies (and id) created on GDV" +
		eol+"    return [{id:<>,species:<>,assembly:<>},...]" +
		eol+"- status : get the status of a job" +
		eol+"    return {job_id:<>,status:<'running','error' or 'success'>}" +
		eol+"" +
		eol+"### PARAMETERS ###" +
		eol+""+
		eol+"## requested for all queries:" +
		eol+"\t- mail : login in tequila" +
		eol+"\t- key  : user key" +
		eol+"\t- command : one of [new_project, add_track,add_sqlite]" +
		eol+""+
		eol+"## command 'new_project'" +
		eol+"\t- seq_id : the nr_assembly id of the species in Genrep" +
		eol+"\t- name : name of the project" +
		eol+"\t- public (OPTIONNAL default=false): true to make the project public" +
		eol+""+
		eol+"## command 'add_track' :" +
		eol+"\t- url : url to fetch" +
		eol+"\t- project_id : the project id in GDV" +
		eol+"\t- name (OPTIONNAL default=name of the file): the name you want to give to the track" +
		eol+""+
		eol+"## command 'status' :" +
		eol+"\t- job_id : the identifier of the job" +
		eol+"";
}

